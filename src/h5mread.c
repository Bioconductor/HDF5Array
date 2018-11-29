/****************************************************************************
 *            Exploring alternate rhdf5::h5read() implementations           *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************
 *
 * Some useful links:
 * - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5S.html
 * - Documentation of H5Dread():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
 * - An H5Dread() example:
 *     https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
 */
#include "HDF5Array.h"
#include "hdf5.h"

#include <zlib.h>

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */
#include <limits.h>  /* for INT_MAX */

//#include <time.h>


typedef struct {
	hid_t dset_id, space_id, plist_id;
	int ndim;
	hsize_t *h5dim, *h5chunk_spacings;
	int *h5nchunk;
	hid_t dtype_id;
	H5T_class_t class;
	size_t size;
	SEXPTYPE ans_type;
	size_t ans_elt_size, chunk_data_buf_size;
	hid_t mem_type_id;  /* the memory type we'll use to load the data */
} DSetDesc;


/****************************************************************************
 * Low-level helpers
 */

/* See hdf5-1.10.3/src/H5Tpublic.h for the list of datatype classes. We only
   support H5T_INTEGER, H5T_FLOAT, and H5T_STRING for now. */
static int get_ans_type_from_class(H5T_class_t class, int as_int, size_t size,
				   SEXPTYPE *ans_type)
{
	const char *classname;

	switch (class) {
	    case H5T_INTEGER:
		*ans_type = (as_int || size <= sizeof(int)) ? INTSXP : REALSXP;
		return 0;
	    case H5T_FLOAT:
		*ans_type = as_int ? INTSXP : REALSXP;
		return 0;
	    case H5T_STRING:
		if (as_int)
			warning("'as.integer' is ignored when "
				"dataset class is H5T_STRING");
		*ans_type = STRSXP;
		return 0;
	    case H5T_TIME: classname = "H5T_TIME"; break;
	    case H5T_BITFIELD: classname = "H5T_BITFIELD"; break;
	    case H5T_OPAQUE: classname = "H5T_OPAQUE"; break;
	    case H5T_COMPOUND: classname = "H5T_COMPOUND"; break;
	    case H5T_REFERENCE: classname = "H5T_REFERENCE"; break;
	    case H5T_ENUM: classname = "H5T_ENUM"; break;
	    case H5T_VLEN: classname = "H5T_VLEN"; break;
	    case H5T_ARRAY: classname = "H5T_ARRAY"; break;
	    default:
		PRINT_TO_ERRMSG_BUF("unknown dataset class identifier: %d",
				    class);
	    return -1;
	}
	PRINT_TO_ERRMSG_BUF("unsupported dataset class: %s", classname);
	return -1;
}

static size_t get_ans_elt_size_from_ans_type(SEXPTYPE ans_type, size_t size)
{
	switch (ans_type) {
	    case INTSXP:  return sizeof(int);
	    case REALSXP: return sizeof(double);
	    case STRSXP:  return size;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(ans_type)));
	return 0;
}

static hid_t get_mem_type_id_from_ans_type(SEXPTYPE ans_type, hid_t dtype_id)
{
	switch (ans_type) {
	    case INTSXP:  return H5T_NATIVE_INT;
	    case REALSXP: return H5T_NATIVE_DOUBLE;
	    case STRSXP:  return dtype_id;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(ans_type)));
	return -1;
}

static hsize_t *alloc_hsize_t_buf(size_t buflength, const char *what)
{
	hsize_t *buf;

	buf = (hsize_t *) malloc(buflength * sizeof(hsize_t));
	if (buf == NULL)
		PRINT_TO_ERRMSG_BUF("failed to allocate memory for %s", what);
	return buf;
}

static void destroy_dset_desc(DSetDesc *dset_desc)
{
	free(dset_desc->h5nchunk);
	free(dset_desc->h5dim);
	H5Tclose(dset_desc->dtype_id);
	H5Pclose(dset_desc->plist_id);
	H5Sclose(dset_desc->space_id);
}

static int get_dset_desc(hid_t dset_id, int ndim, int as_int,
			 DSetDesc *dset_desc)
{
	hid_t space_id, plist_id, dtype_id, mem_type_id;
	int dset_ndim, *h5nchunk, h5along;
	hsize_t *h5dim, *h5chunk_spacings, d, spacing, nchunk;
	H5T_class_t class;
	size_t size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE ans_type;

	dset_desc->dset_id = dset_id;

	/* Set 'space_id'. */
	space_id = H5Dget_space(dset_id);
	if (space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		return -1;
	}
	dset_desc->space_id = space_id;

	/* Set 'ndim'. */
	dset_ndim = H5Sget_simple_extent_ndims(space_id);
	if (dset_ndim < 0) {
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF(
			"H5Sget_simple_extent_ndims() returned an error");
		return -1;
	}
	if (ndim != dset_ndim) {
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF(
			"Dataset has %d dimensions but 'starts' has %d list "
			"element%s.\n  'starts' must have one list element "
			"per dimension in the dataset.",
			dset_ndim, ndim, ndim > 1 ? "s" : "");
		return -1;
	}
	dset_desc->ndim = ndim;

	/* Set 'plist_id'. */
	plist_id = H5Dget_create_plist(dset_id);
	if (plist_id < 0) {
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		return -1;
	}
	dset_desc->plist_id = plist_id;

	/* Set 'dtype_id'. */
	dtype_id = H5Dget_type(dset_id);
	if (dtype_id < 0) {
		H5Pclose(plist_id);
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		return -1;
	}
	dset_desc->dtype_id = dtype_id;

	/* Allocate 'h5dim', 'h5chunk_spacings', and 'h5nchunk'. */
	h5dim = alloc_hsize_t_buf(2 * ndim, "'h5dim' and 'h5chunk_spacings'");
	if (h5dim == NULL) {
		H5Tclose(dtype_id);
		H5Pclose(plist_id);
		H5Sclose(space_id);
		return -1;
	}
	h5chunk_spacings = h5dim + ndim;
	h5nchunk = (int *) malloc(ndim * sizeof(int));
	if (h5nchunk == NULL) {
		free(h5dim);
		H5Tclose(dtype_id);
		H5Pclose(plist_id);
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'h5nchunk'");
		return -1;
	}

	/* Set 'h5dim'. */
	if (H5Sget_simple_extent_dims(space_id, h5dim, NULL) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		goto on_error;
	}
	dset_desc->h5dim = h5dim;

	/* Set 'h5chunk_spacings'. */
	if (H5Pget_chunk(plist_id, ndim, h5chunk_spacings) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
				    "an unexpected value");
		goto on_error;
	}
	dset_desc->h5chunk_spacings = h5chunk_spacings;

	/* Set 'h5nchunk'. */
	for (h5along = 0; h5along < ndim; h5along++) {
		d = h5dim[h5along];
		if (d == 0) {
			h5nchunk[h5along] = 0;
			continue;
		}
		spacing = h5chunk_spacings[h5along];
		nchunk = d / spacing;
		if (d % spacing != 0)
			nchunk++;
		if (nchunk > INT_MAX) {
			PRINT_TO_ERRMSG_BUF("datasets with more than "
				"INT_MAX chunks along any dimension\n  "
				"are not supported at the moment");
			goto on_error;
		}
		h5nchunk[h5along] = nchunk;
	}
	dset_desc->h5nchunk = h5nchunk;

	/* Set 'class'. */
	class = H5Tget_class(dtype_id);
	if (class == H5T_NO_CLASS) {
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		goto on_error;
	}
	dset_desc->class = class;

	/* Set 'size'. */
	size = H5Tget_size(dtype_id);
	if (size == 0) {
		PRINT_TO_ERRMSG_BUF("H5Tget_size() returned 0");
		goto on_error;
	}
	dset_desc->size = size;

	/* Set 'ans_type'. */
	if (get_ans_type_from_class(class, as_int, size, &ans_type) < 0)
		goto on_error;
	if (ans_type == STRSXP && H5Tis_variable_str(dtype_id) != 0) {
		PRINT_TO_ERRMSG_BUF("reading variable-length string data "
				    "is not supported at the moment");
		goto on_error;
	}
	dset_desc->ans_type = ans_type;

	/* Set 'ans_elt_size'. */
	ans_elt_size = get_ans_elt_size_from_ans_type(ans_type, size);
	if (ans_elt_size == 0)
		goto on_error;
	dset_desc->ans_elt_size = ans_elt_size;

	/* Set 'chunk_data_buf_size'. */
	chunk_data_buf_size = ans_elt_size;
	for (h5along = 0; h5along < ndim; h5along++)
		chunk_data_buf_size *= h5chunk_spacings[h5along];
	dset_desc->chunk_data_buf_size = chunk_data_buf_size;

	/* Set 'mem_type_id'. */
	mem_type_id = get_mem_type_id_from_ans_type(ans_type, dtype_id);
	if (mem_type_id < 0)
		goto on_error;
	dset_desc->mem_type_id = mem_type_id;
	return 0;

    on_error:
	destroy_dset_desc(dset_desc);
	return -1;
}

static hid_t get_mem_space(int ndim, const int *ans_dim)
{
	hsize_t *h5dim;
	int along, h5along;
	hid_t mem_space_id;

	/* Allocate and set 'h5dim'. */
	h5dim = alloc_hsize_t_buf(ndim, "'h5dim'");
	if (h5dim == NULL)
		return -1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		h5dim[h5along] = ans_dim[along];

	mem_space_id = H5Screate_simple(ndim, h5dim, NULL);
	if (mem_space_id < 0)
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
	free(h5dim);
	return mem_space_id;
}

static int check_selection_against_dset(const DSetDesc *dset_desc,
			SEXP starts, SEXP counts, int *ans_dim)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = dset_desc->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) dset_desc->h5dim[h5along];
	return _check_selection(starts, counts, dim_buf->elts, ans_dim);
}

static int check_ordered_selection_against_dset(const DSetDesc *dset_desc,
			SEXP starts, SEXP counts, int *ans_dim,
			int *nstart, int *nblock,
			long long int *last_block_start)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = dset_desc->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) dset_desc->h5dim[h5along];
	return _check_ordered_selection(starts, counts, dim_buf->elts, ans_dim,
					nstart, nblock, last_block_start);
}

static int map_starts_to_chunks(const DSetDesc *dset_desc,
				SEXP starts,
				int *ans_dim,
				IntAEAE *breakpoint_bufs,
				LLongAEAE *chunkidx_bufs)
{
	int ndim, along, h5along;
	LLongAE *dim_buf, *chunkdim_buf;

	ndim = dset_desc->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	chunkdim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		dim_buf->elts[along] =
			(long long int) dset_desc->h5dim[h5along];
		chunkdim_buf->elts[along] =
			(long long int) dset_desc->h5chunk_spacings[h5along];
	}
	return _map_starts_to_chunks(starts, dim_buf->elts, chunkdim_buf->elts,
				     ans_dim,
				     breakpoint_bufs, chunkidx_bufs);
}

static inline int next_midx(int ndim, const int *max_idx_plus_one,
			    int *midx_buf)
{
	int along, i;

	for (along = 0; along < ndim; along++) {
		i = midx_buf[along] + 1;
		if (i < max_idx_plus_one[along]) {
			midx_buf[along] = i;
			break;
		}
		midx_buf[along] = 0;
	}
	return along;
}


/****************************************************************************
 * read_data_1_2()
 *
 * A single call to H5Dread().
 *
 * More precisely:
 *   - First walk over the blocks described by 'starts' and 'counts' and
 *     add each block to the selection.
 *   - Then make a single call to H5Dread().
 */

static size_t set_nblock_buf(int ndim, SEXP starts,
			     int expand, const int *ans_dim, int *nblock_buf)
{
	size_t num_blocks;
	int along, nblock;
	SEXP start;

	num_blocks = 1;
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		if (start != R_NilValue) {
			nblock = LENGTH(start);
		} else {
			nblock = expand ? ans_dim[along] : 1;
		}
		num_blocks *= nblock_buf[along] = nblock;
	}
	return num_blocks;
}

static void init_h5blockoff_and_h5blockdim_bufs(const DSetDesc *dset_desc,
			SEXP starts,
			hsize_t *h5blockoff_buf, hsize_t *h5blockdim_buf)
{
	int ndim, along, h5along;

	ndim = dset_desc->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			h5blockoff_buf[h5along] = 0;
			h5blockdim_buf[h5along] = dset_desc->h5dim[h5along];
		} else {
			h5blockdim_buf[h5along] = 1;
		}
	}
	return;
}

static void update_h5blockoff_and_h5blockdim_bufs(int ndim,
			const int *midx, int moved_along,
			SEXP starts, SEXP counts,
			hsize_t *h5blockoff_buf, hsize_t *h5blockdim_buf)
{
	int along, h5along, i;
	SEXP start, count;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		h5blockoff_buf[h5along] = _get_trusted_elt(start, i) - 1;
		if (counts == R_NilValue)
			continue;
		count = VECTOR_ELT(counts, along);
		if (count == R_NilValue)
			continue;
		h5blockdim_buf[h5along] = _get_trusted_elt(count, i);
	}
	return;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_hyperslabs(const DSetDesc *dset_desc,
			SEXP starts, SEXP counts, const int *ans_dim,
			int *nblock_buf, int *midx_buf)
{
	int ret, ndim, moved_along;
	hsize_t *h5blockoff_buf, *h5blockdim_buf;
	long long int num_hyperslabs;

	ret = H5Sselect_none(dset_desc->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = dset_desc->ndim;
	set_nblock_buf(ndim, starts, 0, ans_dim, nblock_buf);

	/* Allocate 'h5blockoff_buf' and 'h5blockdim_buf'. */
	h5blockoff_buf = alloc_hsize_t_buf(2 * ndim, "'h5blockoff_buf' and "
						     "'h5blockdim_buf'");
	if (h5blockoff_buf == NULL)
		return -1;
	h5blockdim_buf = h5blockoff_buf + ndim;

	init_h5blockoff_and_h5blockdim_bufs(dset_desc, starts,
			h5blockoff_buf, h5blockdim_buf);

	/* Walk on the hyperslabs. */
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		update_h5blockoff_and_h5blockdim_bufs(ndim,
				midx_buf, moved_along,
				starts, counts,
				h5blockoff_buf, h5blockdim_buf);
		/* Add to current selection. */
		ret = H5Sselect_hyperslab(dset_desc->space_id, H5S_SELECT_OR,
				h5blockoff_buf, NULL, h5blockdim_buf, NULL);
		if (ret < 0) {
			PRINT_TO_ERRMSG_BUF(
				"H5Sselect_hyperslab() returned an error");
			break;
		}
		moved_along = next_midx(ndim, nblock_buf, midx_buf);
	} while (moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);

	free(h5blockoff_buf);
	return ret < 0 ? -1 : num_hyperslabs;
}

static inline hsize_t *add_element(int ndim, const int *outer_midx,
				   SEXP starts, hsize_t *coord_p)
{
	int h5along, i;
	SEXP start;
	long long int coord;

	for (h5along = ndim - 1; h5along >= 0; h5along--) {
		i = outer_midx[h5along];
		start = VECTOR_ELT(starts, h5along);
		if (start == R_NilValue) {
			coord = i;
		} else {
			coord = _get_trusted_elt(start, i) - 1;
		}
		*(coord_p++) = (hsize_t) coord;
	}
	return coord_p;
}

/* Return nb of selected elements (or -1 on error). */
static long long int select_elements(const DSetDesc *dset_desc,
			SEXP starts, const int *ans_dim,
			int *nblock_buf, int *midx_buf)
{
	int ndim, outer_moved_along, ret;
	size_t num_elements;
	hsize_t *coord_buf, *coord_p;

	ndim = dset_desc->ndim;
	num_elements = set_nblock_buf(ndim, starts, 1, ans_dim, nblock_buf);

	/* Allocate 'coord_buf'. */
	coord_buf = alloc_hsize_t_buf(num_elements * ndim, "'coord_buf'");
	if (coord_buf == NULL)
		return -1;

	/* Walk on the selected elements. */
	coord_p = coord_buf;
	outer_moved_along = ndim;
	do {
		coord_p = add_element(ndim, midx_buf, starts, coord_p);
		outer_moved_along = next_midx(ndim, nblock_buf, midx_buf);
	} while (outer_moved_along < ndim);

	ret = H5Sselect_elements(dset_desc->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	free(coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

static int set_selection(const DSetDesc *dset_desc, int method,
			 SEXP starts, SEXP counts, const int *ans_dim)
{
	int ndim;
	IntAE *nblock_buf, *midx_buf;
	long long int num_blocks;

	ndim = dset_desc->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	//clock_t t0 = clock();
	if (method == 2) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "'method' is set to 2");
			return -1;
		}
		num_blocks = select_elements(dset_desc,
					starts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	} else {
		num_blocks = select_hyperslabs(dset_desc,
					starts, counts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	}
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("time for setting selection: %e\n", dt);
	//printf("num_blocks: %lld, time per block: %e\n",
	//	num_blocks, dt / num_blocks);
	return num_blocks < 0 ? -1 : 0;
}

static int read_data_1_2(const DSetDesc *dset_desc, int method,
			 SEXP starts, SEXP counts, const int *ans_dim,
			 void *out, hid_t mem_space_id)
{
	int ret;

	ret = set_selection(dset_desc, method,
			    starts, counts, ans_dim);
	if (ret < 0)
		return -1;

	ret = H5Sselect_all(mem_space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_all() returned an error");
		return -1;
	}

	//clock_t t0 = clock();
	ret = H5Dread(dset_desc->dset_id,
		      dset_desc->mem_type_id, mem_space_id,
		      dset_desc->space_id, H5P_DEFAULT, out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	//printf("time for reading data from selection: %e\n",
	//	(1.0 * clock() - t0) / CLOCKS_PER_SEC);
	return ret;
}


/****************************************************************************
 * read_data_3()
 *
 * One call to H5Dread() per block in the selection.
 *
 * More precisely, walk over the blocks described by 'starts' and 'counts'
 * and make one call to H5Dread() per block.
 *
 * WARNING: Currently broken!

    library(HDF5Array)
    m0 <- matrix(1:60, ncol=6)
    M0 <- writeHDF5Array(m0, chunkdim=c(4, 5))
    h5mread(M0@seed@filepath, M0@seed@name, list(NULL, c(2:3, 5)), method=3L)
         [,1] [,2]      [,3]
    [1,]   11   41 209047208
    [2,]   12   42         0
    [3,]   13   43 209046632
    [4,]   14   44         0
    [5,]   15   45 209046568
    [6,]   16   46         0
    [7,]   17   47 209046440
    [8,]   18   48         0
    [9,]   19   49 209046312
   [10,]   20   50         0

 */

static int read_hyperslab(const DSetDesc *dset_desc,
		SEXP starts, SEXP counts,
		const int *midx, int moved_along,
		hsize_t *h5offset_in_buf, hsize_t *h5count_buf,
		hsize_t *h5offset_out_buf,
		void *out, hid_t mem_space_id)
{
	int ndim, along, h5along, i, ret;
	SEXP start, count;

	ndim = dset_desc->ndim;

	/* Update 'h5offset_in_buf', 'h5count_buf', and 'h5offset_out_buf'. */
	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		h5offset_in_buf[h5along] = _get_trusted_elt(start, i) - 1;
		if (counts != R_NilValue) {
			count = VECTOR_ELT(counts, along);
			if (count != R_NilValue)
				h5count_buf[h5along] =
					_get_trusted_elt(count, i);
		}
		if (along < moved_along) {
			h5offset_out_buf[h5along] = 0;
		} else {
			h5offset_out_buf[h5along] += h5count_buf[h5along];
		}
	}

	ret = H5Sselect_hyperslab(dset_desc->space_id, H5S_SELECT_SET,
				  h5offset_in_buf, NULL, h5count_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	ret = H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
				  h5offset_out_buf, NULL, h5count_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	ret = H5Dread(dset_desc->dset_id,
		      dset_desc->mem_type_id, mem_space_id,
		      dset_desc->space_id, H5P_DEFAULT, out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int read_data_3(const DSetDesc *dset_desc,
		       SEXP starts, SEXP counts, const int *ans_dim,
		       void *out, hid_t mem_space_id)
{
	int ndim, along, h5along, moved_along, ret;
	hsize_t *h5offset_in_buf, *h5count_buf, *h5offset_out_buf;
	IntAE *nblock_buf, *midx_buf;
	long long int num_hyperslabs;

	ndim = dset_desc->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	set_nblock_buf(ndim, starts, 0, ans_dim, nblock_buf->elts);

	/* Allocate 'h5offset_in_buf', 'h5count_buf', and 'h5offset_out_buf'. */
	h5offset_in_buf = alloc_hsize_t_buf(3 * ndim, "'h5offset_in_buf', "
				"'h5count_buf', and 'h5offset_out_buf'");
	if (h5offset_in_buf == NULL)
		return -1;
	h5count_buf = h5offset_in_buf + ndim;
	h5offset_out_buf = h5count_buf + ndim;

	/* Initialize 'h5offset_in_buf', 'h5count_buf',
	   and 'h5offset_out_buf'. */
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			h5offset_in_buf[h5along] =
				h5offset_out_buf[h5along] = 0;
			h5count_buf[h5along] = dset_desc->h5dim[h5along];
		} else {
			h5count_buf[h5along] = 1;
		}
	}

	/* Walk on the hyperslabs. */
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		ret = read_hyperslab(dset_desc,
				     starts, counts,
				     midx_buf->elts, moved_along,
				     h5offset_in_buf, h5count_buf,
				     h5offset_out_buf,
				     out, mem_space_id);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nblock_buf->elts, midx_buf->elts);
	} while (moved_along < ndim);

	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);
	free(h5offset_in_buf);
	return ret;
}


/****************************************************************************
 * read_data_4_5()
 *
 * One call to H5Dread() per chunk touched by the selection.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Make one call to H5Dread() to load the **entire** chunk data to an
 *     intermediate buffer.
 *   - Copy the selected data from the intermediate buffer to the output
 *     array.
 */

static void set_nchunk_buf(const DSetDesc *dset_desc,
			   const SEXP starts, const LLongAEAE *chunkidx_bufs,
			   IntAE *nchunk_buf)
{
	int ndim, along, h5along, nchunk;

	ndim = dset_desc->ndim;
	//printf("nchunk_buf:");
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) != R_NilValue) {
			nchunk = LLongAE_get_nelt(chunkidx_bufs->elts[along]);
		} else {
			nchunk = dset_desc->h5nchunk[h5along];
		}
		nchunk_buf->elts[along] = nchunk;
		//printf(" %d/%d", nchunk_buf->elts[along],
		//		 dset_desc->h5nchunk[h5along]);
	}
	//printf("\n");
	return;
}

static void update_h5chunkoff_and_h5chunkdim_bufs(const DSetDesc *dset_desc,
			const int *midx, int moved_along,
			SEXP starts, const LLongAEAE *chunkidx_bufs,
			hsize_t *h5chunkoff_buf, hsize_t *h5chunkdim_buf)
{
	int ndim, along, h5along, i;
	long long int chunkidx;
	hsize_t spacing;

	ndim = dset_desc->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = midx[along];
		if (VECTOR_ELT(starts, along) != R_NilValue) {
			chunkidx = chunkidx_bufs->elts[along]->elts[i];
		} else {
			chunkidx = i;
		}
		spacing = dset_desc->h5chunk_spacings[h5along];
		h5chunkoff_buf[h5along] = chunkidx * spacing;
		h5chunkdim_buf[h5along] = dset_desc->h5dim[h5along] -
					  h5chunkoff_buf[h5along];
		if (h5chunkdim_buf[h5along] > spacing)
			h5chunkdim_buf[h5along] = spacing;
	}
	//printf("# h5chunkoff_buf:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//	printf(" %llu", h5chunkoff_buf[h5along]);
	//printf("\n");
	//printf("# h5chunkdim_buf:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//	printf(" %llu", h5chunkdim_buf[h5along]);
	//printf("\n");
	return;
}

static void update_out_blockoff_and_out_blockdim_bufs(int ndim,
			const int *midx, int moved_along,
			SEXP starts, const IntAEAE *breakpoint_bufs,
			const hsize_t *h5chunkoff, const hsize_t *h5chunkdim,
			int *out_blockoff_buf, int *out_blockdim_buf)
{
	int along, h5along, i, off, d;
	const int *breakpoint;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = midx[along];
		if (VECTOR_ELT(starts, along) != R_NilValue ) {
			breakpoint = breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = h5chunkoff[h5along];
			d = h5chunkdim[h5along];
		}
		out_blockoff_buf[along] = off;
		out_blockdim_buf[along] = d;
	}
	//printf("# out_blockoff_buf:");
	//for (along = 0; along < ndim; along++)
	//	printf(" %d", out_blockoff_buf[along]);
	//printf("\n");
	//printf("# out_blockdim_buf:");
	//for (along = 0; along < ndim; along++)
	//	printf(" %d", out_blockdim_buf[along]);
	//printf("\n");
	return;
}

/* It takes about 218s on my laptop to load all the chunks from the EH1040
 * dataset (big 10x Genomics brain dataset in dense format, chunks of 100x100,
 * wrapped in the TENxBrainData package). That's 60 microseconds per chunk! */
static int load_chunk(const DSetDesc *dset_desc,
		      const hsize_t *h5chunkoff_buf,
		      const hsize_t *h5chunkdim_buf,
		      const hsize_t *h5zeros_buf,
		      void *chunk_data_out, hid_t mem_space_id)
{
	int ret;

	ret = H5Sselect_hyperslab(dset_desc->space_id, H5S_SELECT_SET,
				  h5chunkoff_buf, NULL, h5chunkdim_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	ret = H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
				  h5zeros_buf, NULL, h5chunkdim_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	ret = H5Dread(dset_desc->dset_id,
		      dset_desc->mem_type_id, mem_space_id,
		      dset_desc->space_id, H5P_DEFAULT, chunk_data_out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int uncompress_chunk_data(const void *compressed_chunk_data,
				 hsize_t compressed_size,
				 void *uncompressed_chunk_data,
				 size_t uncompressed_size)
{
	int ret;
	uLong destLen;

	destLen = uncompressed_size;
	ret = uncompress((Bytef *) uncompressed_chunk_data, &destLen,
			 compressed_chunk_data, (uLong) compressed_size);
	if (ret == Z_OK) {
		if (destLen == uncompressed_size)
			return 0;
		PRINT_TO_ERRMSG_BUF("error in uncompress_chunk_data(): "
				    "chunk data smaller than expected "
				    "after decompression");
		return -1;
	}
	switch (ret) {
	    case Z_MEM_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "not enough memory to uncompress chunk");
	    break;
	    case Z_BUF_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "not enough room in output buffer");
	    break;
	    case Z_DATA_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "chunk data corrupted or incomplete");
	    break;
	    default:
		PRINT_TO_ERRMSG_BUF("unknown error in uncompress()");
	}
	return -1;
}

/*
 * Unfortunately H5Dread_chunk() is NOT listed here:
 *   https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html
 * Header file for declaration:
 *   hdf5-1.10.3/src/H5Dpublic.h
 * See hdf5-1.10.3/test/direct_chunk.c for plenty of examples.
 *
 * Call stack for H5Dread_chunk()
 *   H5Dread_chunk                (H5Dio.c)
 *     H5D__chunk_direct_read     (H5Dchunk.c)
 *       H5F_block_read           (H5Fio.c)
 *         H5PB_read              (H5PB.c)
 *           H5F__accum_read      (H5Faccum.c)
 *             or
 *           H5FD_read            (H5FDint.c)
 *            ??
 *
 * Call stack for H5Dread()
 *   H5Dread                      (H5Dio.c)
 *     H5D__read                  (H5Dio.c)
 *       H5D__chunk_read          (H5Dchunk.c)
 *         H5D__select_read
 *           or
 *         H5D__scatgath_read     (H5Dscatgath.c)
 *           H5D__gather_file     (H5Dscatgath.c)
 *       call ser_read member of a H5D_layout_ops_t object
 *            ??
 */
#define COMPRESSION_OVERHEAD 8	// empirical (increase if necessary)

static int direct_load_chunk(const DSetDesc *dset_desc,
			     const hsize_t *h5chunkoff_buf,
			     const hsize_t *h5chunkdim_buf,
			     void *chunk_data_out,
			     void *compressed_chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(dset_desc->dset_id,
					h5chunkoff_buf, &chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > dset_desc->chunk_data_buf_size +
				 COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%lu + %d)",
				    chunk_storage_size,
				    dset_desc->chunk_data_buf_size,
				    COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(dset_desc->dset_id, H5P_DEFAULT,
			    h5chunkoff_buf, &filters,
			    compressed_chunk_data_buf);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dread_chunk() returned an error");
		return -1;
	}

	//printf("filters = %u\n", filters);

	//FIXME: This will error if chunk data is not compressed!
	//TODO: Decompress only if chunk data is compressed. There should be
	// a bit in the returned 'filters' that indicates this.
	return uncompress_chunk_data(compressed_chunk_data_buf,
				     chunk_storage_size,
				     chunk_data_out,
				     dset_desc->chunk_data_buf_size);
}

static void init_in_offset_and_out_offset(int ndim, SEXP starts,
			const int *outdim, const int *out_blockoff,
			const hsize_t *h5chunkoff,
			const hsize_t *h5chunk_spacings,
			size_t *in_offset, size_t *out_offset)
{
	size_t in_off, out_off;
	int along, h5along, i;
	SEXP start;

	in_off = out_off = 0;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunk_spacings[h5along];
		out_off *= outdim[along];
		i = out_blockoff[along];
		start = VECTOR_ELT(starts, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  h5chunkoff[h5along];
		out_off += i;
	}
	*in_offset = in_off;
	*out_offset = out_off;
	//printf("# in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static inline void update_in_offset_and_out_offset(int ndim,
			const int *inner_midx, int inner_moved_along,
			SEXP starts,
			const int *outdim, const int *out_blockoff,
			const int *out_blockdim,
			const hsize_t *h5chunk_spacings,
			size_t *in_offset, size_t *out_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc, out_off_inc;

	start = VECTOR_ELT(starts, inner_moved_along);
	if (start != R_NilValue) {
		i1 = out_blockoff[inner_moved_along] +
		     inner_midx[inner_moved_along];
		i0 = i1 - 1;
		in_off_inc = _get_trusted_elt(start, i1) -
			     _get_trusted_elt(start, i0);
	} else {
		in_off_inc = 1;
	}
	out_off_inc = 1;
	if (inner_moved_along >= 1) {
		along = inner_moved_along - 1;
		h5along = ndim - inner_moved_along;
		do {
			in_off_inc *= h5chunk_spacings[h5along];
			out_off_inc *= outdim[along];
			di = 1 - out_blockdim[along];
			start = VECTOR_ELT(starts, along);
			if (start != R_NilValue) {
				i1 = out_blockoff[along];
				i0 = i1 - di;
				in_off_inc += _get_trusted_elt(start, i1) -
					      _get_trusted_elt(start, i0);
			} else {
				in_off_inc += di;
			}
			out_off_inc += di;
			along--;
			h5along++;
		} while (along >= 0);
	}
	*in_offset += in_off_inc;
	*out_offset += out_off_inc;
	//printf("## in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static void gather_chunk_data(const DSetDesc *dset_desc,
			const int *outer_midx, int outer_moved_along,
			SEXP starts, const IntAEAE *breakpoint_bufs,
			const void *in,
			const hsize_t *h5chunkoff,
			SEXP ans, const int *outdim,
			const int *out_blockoff,
			const int *out_blockdim,
			int *inner_midx_buf,
			hid_t mem_type_id)
{
	int ndim, inner_moved_along;
	size_t in_offset, out_offset, s_len;
	long long int num_elts;
	const char *s;
	SEXP ans_elt;

	ndim = dset_desc->ndim;
	init_in_offset_and_out_offset(ndim, starts,
			outdim, out_blockoff,
			h5chunkoff, dset_desc->h5chunk_spacings,
			&in_offset, &out_offset);

	/* Walk on the selected elements in current chunk. */
	num_elts = 0;
	while (1) {
		num_elts++;
		if (dset_desc->ans_type == STRSXP) {
			s = (char *) in + in_offset * dset_desc->size;
			for (s_len = 0; s_len < dset_desc->size; s_len++)
				if (s[s_len] == 0)
					break;
			ans_elt = PROTECT(mkCharLen(s, s_len));
			SET_STRING_ELT(ans, out_offset, ans_elt);
			UNPROTECT(1);
		} else {
			if (dset_desc->ans_type == INTSXP) {
				INTEGER(ans)[out_offset] =
					((int *) in)[in_offset];
			} else {
				REAL(ans)[out_offset] =
					((double *) in)[in_offset];
			}
		}

		inner_moved_along = next_midx(ndim, out_blockdim,
					      inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset_and_out_offset(ndim,
				inner_midx_buf, inner_moved_along,
				starts,
				outdim, out_blockoff,
				out_blockdim, dset_desc->h5chunk_spacings,
				&in_offset, &out_offset);
	};
	//printf("# nb of selected elements in current chunk = %lld\n",
	//       num_elts);
	return;
}

/*
  WARNING: method 5 is not working properly on some datasets:
      library(HDF5Array)
      library(ExperimentHub)
      hub <- ExperimentHub()
      fname0 <- hub[["EH1039"]]
      h5mread(fname0, "mm10/barcodes", list(1), method=4L)
      # [1] "AAACCTGAGATAGGAG-1"
      h5mread(fname0, "mm10/barcodes", list(1), method=5L)
      # [1] "AAAAAAAAAAAAAAAAAAAA"
  Looks like the chunk data has been shuffled (transposed in that case)
  before being written to disk in order to improve compression.
  TODO: Investigate this further. I suspect we need to check whether a
  "Data shuffling filter" (H5Z_FILTER_SHUFFLE) was used at creation time.
  Check H5Pget_filter() for how to know whether this filter was used or not.
  There should be a way to retrieve information about how the data was
  shuffled.
 */
static int read_data_4_5(const DSetDesc *dset_desc, int method,
			 SEXP starts,
			 const IntAEAE *breakpoint_bufs,
			 const LLongAEAE *chunkidx_bufs,
			 const int *outdim, SEXP ans)
{
	int ndim, along, moved_along, ret;
	hid_t mem_space_id;
	void *chunk_data_buf, *compressed_chunk_data_buf = NULL;
	hsize_t *h5chunkoff_buf, *h5chunkdim_buf, *h5zeros_buf;
	IntAE *nchunk_buf, *outer_midx_buf,
	      *out_blockoff_buf, *out_blockdim_buf, *inner_midx_buf;
	long long int num_chunks, ndot;

	ndim = dset_desc->ndim;
	mem_space_id = H5Screate_simple(ndim, dset_desc->h5chunk_spacings,
					NULL);
	if (mem_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}

	/* Allocate 'h5chunkoff_buf', 'h5chunkdim_buf', and 'h5zeros_buf'. */
	h5chunkoff_buf = alloc_hsize_t_buf(3 * ndim, "'h5chunkoff_buf', "
				"'h5chunkdim_buf', and 'h5zeros_buf'");
	if (h5chunkoff_buf == NULL) {
		H5Sclose(mem_space_id);
		return -1;
	}
	h5chunkdim_buf = h5chunkoff_buf + ndim;
	h5zeros_buf = h5chunkdim_buf + ndim;
	for (along = 0; along < ndim; along++)
		h5zeros_buf[along] = 0;

	if (method == 4) {
		chunk_data_buf = malloc(dset_desc->chunk_data_buf_size);
	} else {
		warning("method 5 is still experimental, use at your own risk");
		chunk_data_buf = malloc(2 * dset_desc->chunk_data_buf_size +
					COMPRESSION_OVERHEAD);
	}
	if (chunk_data_buf == NULL) {
		free(h5chunkoff_buf);
		H5Sclose(mem_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}
	if (method != 4)
		compressed_chunk_data_buf = chunk_data_buf +
					    dset_desc->chunk_data_buf_size;

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(dset_desc, starts, chunkidx_bufs, nchunk_buf);
	outer_midx_buf = new_IntAE(ndim, ndim, 0);
	out_blockoff_buf = new_IntAE(ndim, ndim, 0);
	out_blockdim_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	/* Walk on the chunks. */
	num_chunks = ndot = 0;
	moved_along = ndim;
	//clock_t t_load_chunk = 0, t_gather_chunk_data = 0, t0;
	do {
		num_chunks++;
		//if (num_chunks % 20000 == 0) {
		//	printf(".");
		//	if (++ndot % 50 == 0)
		//		printf(" [%lld chunks processed]\n",
		//		       num_chunks);
		//	fflush(stdout);
		//}
		update_h5chunkoff_and_h5chunkdim_bufs(dset_desc,
			outer_midx_buf->elts, moved_along,
			starts, chunkidx_bufs,
			h5chunkoff_buf, h5chunkdim_buf);

		//t0 = clock();
		if (method == 4) {
			ret = load_chunk(dset_desc,
					 h5chunkoff_buf, h5chunkdim_buf,
					 h5zeros_buf,
					 chunk_data_buf, mem_space_id);
		} else {
			ret = direct_load_chunk(dset_desc,
					h5chunkoff_buf, h5chunkdim_buf,
					chunk_data_buf,
					compressed_chunk_data_buf);
		}
		if (ret < 0)
			break;
		//t_load_chunk += clock() - t0;

		//t0 = clock();
		update_out_blockoff_and_out_blockdim_bufs(ndim,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			h5chunkoff_buf, h5chunkdim_buf,
			out_blockoff_buf->elts, out_blockdim_buf->elts);

		gather_chunk_data(dset_desc,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			chunk_data_buf,
			h5chunkoff_buf,
			ans, outdim,
			out_blockoff_buf->elts, out_blockdim_buf->elts,
			inner_midx_buf->elts,
			dset_desc->mem_type_id);
		//t_gather_chunk_data += clock() - t0;

		moved_along = next_midx(ndim, nchunk_buf->elts,
					outer_midx_buf->elts);
	} while (moved_along < ndim);
	//printf("\n");
	//printf("total chunks processed = %lld\n", num_chunks);
	//double dt1 = 1.0 * t_load_chunk / CLOCKS_PER_SEC;
	//printf("total time spent loading chunks: %3.6f\n", dt1);
	//double dt2 = 1.0 * t_gather_chunk_data / CLOCKS_PER_SEC;
	//printf("total time spent gathering the chunk data: %3.6f\n", dt2);
	//printf("total time: %3.6f\n", dt1 + dt2);

	free(chunk_data_buf);
	free(h5chunkoff_buf);
	H5Sclose(mem_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to H5Dread() per chunk touched by the selection. No intermediate
 * buffer.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk
 * make one call to H5Dread() to load the selected data **directly** to the
 * output array.
 */

/*
static hsize_t *alloc_coord_buf(int ndim, const hsize_t *h5chunk_spacings)
{
	size_t max_inner_block_per_chunk;
	int along;

	max_inner_block_per_chunk = 1;
	for (along = 0; along < ndim; along++)
		max_inner_block_per_chunk *= (h5chunk_spacings[along] + 1) / 2;
	//printf("max_inner_block_per_chunk = %lu\n",
	//	 max_inner_block_per_chunk);
	return alloc_hsize_t_buf(max_inner_block_per_chunk * ndim,
				 "'coord_buf'");
}
*/

static void update_out_h5blockoff_and_out_h5blockdim_bufs(int ndim,
			const int *midx, int moved_along,
			SEXP starts, const IntAEAE *breakpoint_bufs,
			const hsize_t *h5chunkoff, const hsize_t *h5chunkdim,
			hsize_t *out_h5blockoff_buf,
			hsize_t *out_h5blockdim_buf)
{
	int along, h5along, i, off, d;
	const int *breakpoint;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = midx[along];
		if (VECTOR_ELT(starts, along) != R_NilValue ) {
			breakpoint = breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = h5chunkoff[h5along];
			d = h5chunkdim[h5along];
		}
		out_h5blockoff_buf[h5along] = off;
		out_h5blockdim_buf[h5along] = d;
	}
	return;
}

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP starts,
		const int *out_blockoff, const int *out_blockdim,
		IntAEAE *inner_breakpoint_bufs, int *inner_nblock_buf)
{
	int along, d, off, nblock, i;
	IntAE *inner_breakpoint_buf;
	SEXP start;
	long long int s0, s1;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		inner_breakpoint_buf = inner_breakpoint_bufs->elts[along];
		IntAE_set_nelt(inner_breakpoint_buf, 0);
		d = out_blockdim[along];
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue) {
			IntAE_insert_at(inner_breakpoint_buf, 0, d);
			inner_nblock_buf[along] = 1;
			continue;
		}
		off = out_blockoff[along];
		s1 = _get_trusted_elt(start, off);
		nblock = 0;
		for (i = 1; i < d; i++) {
			s0 = s1;
			s1 = _get_trusted_elt(start, off + i);
			if (s1 != s0 + 1)
				IntAE_insert_at(inner_breakpoint_buf,
						nblock++, i);
		}
		IntAE_insert_at(inner_breakpoint_buf, nblock++, d);
		inner_nblock_buf[along] = nblock;
	}
	return;
}

static void init_in_h5blockoff_and_in_h5blockdim_bufs(
			int ndim, SEXP starts,
			const hsize_t *h5chunkoff,
			const hsize_t *h5chunkdim,
			hsize_t *in_h5blockoff_buf,
			hsize_t *in_h5blockdim_buf)
{
	int along, h5along;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			in_h5blockoff_buf[h5along] = h5chunkoff[h5along];
			in_h5blockdim_buf[h5along] = h5chunkdim[h5along];
		} else {
			in_h5blockdim_buf[h5along] = 1;
		}
	}
	return;
}

static void update_in_h5blockoff_and_in_h5blockdim_bufs(
			int ndim,
			SEXP starts, const int *out_blockoff,
			const int *inner_midx, int inner_moved_along,
			const IntAEAE *inner_breakpoint_bufs,
			hsize_t *in_h5blockoff_buf, hsize_t *in_h5blockdim_buf)
{
	int along, h5along, idx, off, d, i;
	SEXP start;
	const int *inner_breakpoint;

	for (along = 0; along < ndim; along++) {
		if (along > inner_moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		inner_breakpoint = inner_breakpoint_bufs->elts[along]->elts;
		idx = inner_midx[along];
		off = idx == 0 ? 0 : inner_breakpoint[idx - 1];
		d = inner_breakpoint[idx] - off;
		i = out_blockoff[along] + off;
		h5along = ndim - 1 - along;
		in_h5blockoff_buf[h5along] = _get_trusted_elt(start, i) - 1;
		in_h5blockdim_buf[h5along] = d;
	}
	return;
}

/* Return nb of selected elements (or -1 on error). */
static long long int select_elements_from_chunk(const DSetDesc *dset_desc,
			SEXP starts,
			const hsize_t *out_h5blockoff, const int *out_blockdim,
			int *inner_midx_buf, hsize_t *coord_buf)
{
	int ret, ndim, inner_moved_along, along, h5along, i;
	size_t num_elements;
	hsize_t *coord_p;
	long long int coord;
	SEXP start;

	ret = H5Sselect_none(dset_desc->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = dset_desc->ndim;

	/* Walk on the selected elements from the current chunk. */
	num_elements = 0;
	coord_p = coord_buf;
	inner_moved_along = ndim;
	do {
		num_elements++;
		for (along = ndim - 1, h5along = 0;
		     along >= 0;
		     along--, h5along++)
		{
			i = out_h5blockoff[h5along] + inner_midx_buf[along];
			start = VECTOR_ELT(starts, along);
			if (start == R_NilValue) {
				coord = i;
			} else {
				coord = _get_trusted_elt(start, i) - 1;
			}
			*(coord_p++) = (hsize_t) coord;
		}
		inner_moved_along = next_midx(ndim, out_blockdim,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);

	ret = H5Sselect_elements(dset_desc->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_hyperslabs_from_chunk(const DSetDesc *dset_desc,
			SEXP starts, const int *out_blockoff,
			const hsize_t *h5chunkoff,
			const hsize_t *h5chunkdim,
			const IntAEAE *inner_breakpoint_bufs,
			const int *inner_nblock,
			int *inner_midx_buf,
			hsize_t *in_h5blockoff_buf,
			hsize_t *in_h5blockdim_buf)
{
	int ret, ndim, inner_moved_along;
	long long int num_hyperslabs;

	ret = H5Sselect_none(dset_desc->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = dset_desc->ndim;

	init_in_h5blockoff_and_in_h5blockdim_bufs(ndim, starts,
			h5chunkoff, h5chunkdim,
			in_h5blockoff_buf, in_h5blockdim_buf);

	/* Walk on the inner blocks. */
	num_hyperslabs = 0;
	inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_in_h5blockoff_and_in_h5blockdim_bufs(ndim,
				starts, out_blockoff,
				inner_midx_buf, inner_moved_along,
				inner_breakpoint_bufs,
				in_h5blockoff_buf, in_h5blockdim_buf);
		ret = H5Sselect_hyperslab(dset_desc->space_id, H5S_SELECT_OR,
				in_h5blockoff_buf, NULL,
				in_h5blockdim_buf, NULL);
		if (ret < 0) {
			PRINT_TO_ERRMSG_BUF(
				"H5Sselect_hyperslab() returned an error");
			return -1;
		}
		inner_moved_along = next_midx(ndim, inner_nblock,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(inner_nblock)

	return num_hyperslabs;
}

static int read_selection(const DSetDesc *dset_desc,
			  const hsize_t *out_h5blockoff,
			  const hsize_t *out_h5blockdim,
			  void *out, hid_t mem_space_id)
{
	int ret;

	ret = H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
				  out_h5blockoff, NULL, out_h5blockdim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	ret = H5Dread(dset_desc->dset_id,
		      dset_desc->mem_type_id, mem_space_id,
		      dset_desc->space_id, H5P_DEFAULT, out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int read_data_6(const DSetDesc *dset_desc,
		       SEXP starts,
		       const IntAEAE *breakpoint_bufs,
		       const LLongAEAE *chunkidx_bufs,
		       const int *outdim, SEXP ans)
{
	void *out;
	int ndim, moved_along, ret;
	hid_t mem_space_id;
	hsize_t *h5chunkoff_buf, *h5chunkdim_buf,
		*out_h5blockoff_buf, *out_h5blockdim_buf,
		*in_h5blockoff_buf, *in_h5blockdim_buf;
		//*coord_buf;
	IntAE *nchunk_buf, *outer_midx_buf,
	      *out_blockoff_buf, *out_blockdim_buf,
	      *inner_nblock_buf, *inner_midx_buf;
	IntAEAE *inner_breakpoint_bufs;
	long long int num_chunks, ndot;

	if (dset_desc->ans_type == INTSXP) {
		out = INTEGER(ans);
	} else {
		out = REAL(ans);
	}
	ndim = dset_desc->ndim;
	mem_space_id = get_mem_space(ndim, outdim);
	if (mem_space_id < 0)
		return -1;

	/* Allocate 'h5chunkoff_buf', 'h5chunkdim_buf',
	   'out_h5blockoff_buf', 'out_h5blockdim_buf',
	   'in_h5blockoff_buf', and 'in_h5blockdim_buf'. */
	h5chunkoff_buf = alloc_hsize_t_buf(6 * ndim,
				"'h5chunkoff_buf', 'h5chunkdim_buf', "
				"'out_h5blockoff_buf', 'out_h5blockdim_buf', "
				"'in_h5blockoff_buf', and 'in_h5blockdim_buf'");
	if (h5chunkoff_buf == NULL) {
		H5Sclose(mem_space_id);
		return -1;
	}
	h5chunkdim_buf = h5chunkoff_buf + ndim;
	out_h5blockoff_buf = h5chunkdim_buf + ndim;
	out_h5blockdim_buf = out_h5blockoff_buf + ndim;
	in_h5blockoff_buf = out_h5blockdim_buf + ndim;
	in_h5blockdim_buf = in_h5blockoff_buf + ndim;

	/* Allocate 'coord_buf'. */
	//coord_buf = alloc_coord_buf(ndim, dset_desc->h5chunk_spacings);
	//if (coord_buf == NULL) {
	//	free(h5chunkoff_buf);
	//	H5Sclose(mem_space_id);
	//	return -1;
	//}

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(dset_desc, starts, chunkidx_bufs, nchunk_buf);
	outer_midx_buf = new_IntAE(ndim, ndim, 0);
	out_blockoff_buf = new_IntAE(ndim, ndim, 0);
	out_blockdim_buf = new_IntAE(ndim, ndim, 0);
	inner_nblock_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	/* Walk on the chunks. */
	num_chunks = ndot = 0;
	moved_along = ndim;
	//clock_t t_select_elements = 0, t_read_selection = 0, t0;
	do {
		num_chunks++;
		//if (num_chunks % 20000 == 0) {
		//	printf(".");
		//	if (++ndot % 50 == 0)
		//		printf(" [%lld chunks processed]\n",
		//		       num_chunks);
		//	fflush(stdout);
		//}
		update_h5chunkoff_and_h5chunkdim_bufs(dset_desc,
			outer_midx_buf->elts, moved_along,
			starts, chunkidx_bufs,
			h5chunkoff_buf, h5chunkdim_buf);

		update_out_blockoff_and_out_blockdim_bufs(ndim,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			h5chunkoff_buf, h5chunkdim_buf,
			out_blockoff_buf->elts, out_blockdim_buf->elts);

		update_out_h5blockoff_and_out_h5blockdim_bufs(ndim,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			h5chunkoff_buf, h5chunkdim_buf,
			out_h5blockoff_buf, out_h5blockdim_buf);

		update_inner_breakpoints(ndim, moved_along,
			starts,
			out_blockoff_buf->elts, out_blockdim_buf->elts,
			inner_breakpoint_bufs, inner_nblock_buf->elts);

		//t0 = clock();
		/* Having 'inner_nblock_buf->elts' identical to
		   'out_blockdim_buf->elts' means that all the inner blocks
		   are single elements so we use select_elements_from_chunk()
		   which is faster than select_hyperslabs_from_chunk() in that
		   case. NO IT'S NOT FASTER! */
		//ret = memcmp(inner_nblock_buf->elts, out_blockdim_buf->elts,
		//	     ndim * sizeof(int));
		//printf("# chunk %lld: %s\n", num_chunks,
		//       ret == 0 ? "select_elements" : "select_hyperslabs");
		//if (ret == 0) {
		//	ret = select_elements_from_chunk(dset_desc,
		//		starts,
		//		out_h5blockoff_buf, out_blockdim_buf->elts,
		//		inner_midx_buf->elts, coord_buf);
		//} else {
			ret = select_hyperslabs_from_chunk(dset_desc,
				starts, out_blockoff_buf->elts,
				h5chunkoff_buf, h5chunkdim_buf,
				inner_breakpoint_bufs, inner_nblock_buf->elts,
				inner_midx_buf->elts,
				in_h5blockoff_buf, in_h5blockdim_buf);
		//}
		if (ret < 0)
			break;
		//t_select_elements += clock() - t0;

		//t0 = clock();
		ret = read_selection(dset_desc,
			out_h5blockoff_buf, out_h5blockdim_buf,
			out, mem_space_id);
		if (ret < 0)
			break;
		//t_read_selection += clock() - t0;

		moved_along = next_midx(ndim, nchunk_buf->elts,
					outer_midx_buf->elts);
	} while (moved_along < ndim);
	//printf("\n");
	//printf("total chunks processed = %lld\n", num_chunks);
	//double dt1 = 1.0 * t_select_elements / CLOCKS_PER_SEC;
	//printf("total time spent selecting elements: %3.6f\n", dt1);
	//double dt2 = 1.0 * t_read_selection / CLOCKS_PER_SEC;
	//printf("total time spent reading elements: %3.6f\n", dt2);
	//printf("total time: %3.6f\n", dt1 + dt2);

	//free(coord_buf);
	free(h5chunkoff_buf);
	H5Sclose(mem_space_id);
	return ret;
}


/****************************************************************************
 * C_h5mread()
 */

/* Return R__NilValue on error. */
static SEXP h5mread_1_2_3(const DSetDesc *dset_desc,
			  SEXP starts, SEXP counts, int noreduce,
			  int method, int *ans_dim)
{
	int ndim, ret, along;
	IntAE *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	int *nstart, *nblock;
	long long int *last_block_start;
	R_xlen_t ans_len;
	SEXP ans, reduced;
	void *out;
	hid_t mem_space_id;
	int nprotect = 0;

	ndim = dset_desc->ndim;

	if (noreduce || method == 2) {
		/* This call will populate 'ans_dim'. */
		ret = check_selection_against_dset(dset_desc,
					starts, counts, ans_dim);
		if (ret < 0)
			return R_NilValue;
	} else {
		nstart_buf = new_IntAE(ndim, ndim, 0);
		nblock_buf = new_IntAE(ndim, ndim, 0);
		last_block_start_buf = new_LLongAE(ndim, ndim, 0);
		nstart = nstart_buf->elts;
		nblock = nblock_buf->elts;
		last_block_start = last_block_start_buf->elts;
		/* This call will populate 'ans_dim', 'nblock',
		   and 'last_block_start'. */
		ret = check_ordered_selection_against_dset(dset_desc,
					starts, counts, ans_dim,
					nstart, nblock, last_block_start);
		if (ret < 0)
			return R_NilValue;
		if (_selection_can_be_reduced(ndim, nstart, nblock)) {
			reduced = PROTECT(_reduce_selection(starts, counts,
							    ans_dim,
							    nblock,
							    last_block_start));
			nprotect++;
			starts = VECTOR_ELT(reduced, 0);
			counts = VECTOR_ELT(reduced, 1);
		}
	}

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(dset_desc->ans_type, ans_len));
	nprotect++;

	if (ans_len != 0) {
		if (dset_desc->ans_type == INTSXP) {
			out = INTEGER(ans);
		} else {
			out = REAL(ans);
		}
		mem_space_id = get_mem_space(ndim, ans_dim);
		if (mem_space_id < 0)
			goto on_error;
		if (method <= 2) {
			ret = read_data_1_2(dset_desc, method,
				starts, counts, ans_dim,
				out, mem_space_id);
		} else {
			ret = read_data_3(dset_desc,
				starts, counts, ans_dim,
				out, mem_space_id);
		}
		H5Sclose(mem_space_id);
	}
	if (ret < 0)
		goto on_error;

	UNPROTECT(nprotect);
	return ans;

    on_error:
	UNPROTECT(nprotect);
	return R_NilValue;
}

/* Return R__NilValue on error. */
static SEXP h5mread_4_5_6(const DSetDesc *dset_desc,
			  SEXP starts, SEXP counts,
			  int method, int *ans_dim)
{
	int ndim, ret, along;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *chunkidx_bufs;
	R_xlen_t ans_len;
	SEXP ans;

	ndim = dset_desc->ndim;
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunkidx_bufs = new_LLongAEAE(ndim, ndim);

	/* This call will populate 'ans_dim', 'breakpoint_bufs',
	   and 'chunkidx_bufs'. */
	ret = map_starts_to_chunks(dset_desc, starts,
				   ans_dim, breakpoint_bufs, chunkidx_bufs);
	if (ret < 0)
		return R_NilValue;

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(dset_desc->ans_type, ans_len));

	if (ans_len != 0) {
		if (method <= 5) {
			/* methods 4 and 5 */
			ret = read_data_4_5(dset_desc, method,
				starts, breakpoint_bufs, chunkidx_bufs,
				ans_dim, ans);
		} else {
			/* method 6 */
			ret = read_data_6(dset_desc,
				starts, breakpoint_bufs, chunkidx_bufs,
				ans_dim, ans);
		}
		if (ret < 0)
			goto on_error;
	}

	UNPROTECT(1);
	return ans;

    on_error:
	UNPROTECT(1);
	return R_NilValue;
}

/* Return R__NilValue on error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce,
		    int as_int, int method)
{
	int ndim;
	DSetDesc dset_desc;
	SEXP ans, ans_dim;

	if (method < 0 || method > 6) {
		PRINT_TO_ERRMSG_BUF("'method' must be >= 0 and <= 6");
		return R_NilValue;
	}

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0)
		return R_NilValue;

	if (get_dset_desc(dset_id, ndim, as_int, &dset_desc) < 0)
		return R_NilValue;

	ans = R_NilValue;

	if (dset_desc.ans_type == STRSXP) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "reading string data");
			goto on_error;
		}
		if (method == 0) {
			method = 4;
		} else if (method != 4 && method != 5) {
			PRINT_TO_ERRMSG_BUF("reading string data is only "
					    "supported by methods 4 and 5");
			goto on_error;
		}
	} else if (method == 0) {
		method = counts == R_NilValue ? 6 : 1;
	}

	ans_dim = PROTECT(NEW_INTEGER(ndim));

	if (method <= 3) {
		ans = h5mread_1_2_3(&dset_desc, starts, counts, noreduce,
				    method, INTEGER(ans_dim));
	} else if (counts != R_NilValue) {
		PRINT_TO_ERRMSG_BUF("'counts' must be NULL for "
				    "methods 4, 5, and 6");
	} else {
		ans = h5mread_4_5_6(&dset_desc, starts, counts,
				    method, INTEGER(ans_dim));
	}
	if (ans != R_NilValue) {
		PROTECT(ans);
		SET_DIM(ans, ans_dim);
		UNPROTECT(1);
	}

	UNPROTECT(1);

    on_error:
	destroy_dset_desc(&dset_desc);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce,
	       SEXP as_integer, SEXP method)
{
	SEXP filepath0, name0, ans;
	int noreduce0, as_int, method0;
	hid_t file_id, dset_id;

	/* Check 'filepath'. */
	if (!(IS_CHARACTER(filepath) && LENGTH(filepath) == 1))
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' cannot be NA");

	/* Check 'name'. */
	if (!(IS_CHARACTER(name) && LENGTH(name) == 1))
		error("'name' must be a single string");
	name0 = STRING_ELT(name, 0);
	if (name0 == NA_STRING)
		error("'name' cannot be NA");

	/* Check 'noreduce'. */
	if (!(IS_LOGICAL(noreduce) && LENGTH(noreduce) == 1))
		error("'noreduce' must be TRUE or FALSE");
	noreduce0 = LOGICAL(noreduce)[0];

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	/* Check 'method'. */
	if (!(IS_INTEGER(method) && LENGTH(method) == 1))
		error("'method' must be a single integer");
	method0 = INTEGER(method)[0];

	file_id = H5Fopen(CHAR(filepath0), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		error("failed to open file %s", CHAR(filepath0));
	dset_id = H5Dopen(file_id, CHAR(name0), H5P_DEFAULT);
	if (dset_id < 0) {
		H5Fclose(file_id);
		error("failed to open dataset %s from file %s",
		      CHAR(name0), CHAR(filepath0));
	}
	ans = PROTECT(h5mread(dset_id, starts, counts, noreduce0,
			      as_int, method0));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ans == R_NilValue) {
		UNPROTECT(1);
		error(_HDF5Array_errmsg_buf);
	}
	UNPROTECT(1);
	return ans;
}

