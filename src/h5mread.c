/****************************************************************************
 *       Experimenting with alternate rhdf5::h5read() implementations       *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"
#include "hdf5.h"

#include <zlib.h>

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */

//#include <time.h>


typedef struct {
	hid_t dset_id, space_id, plist_id;
	int ndim;
	hsize_t *h5dim, *h5chunkdim;
} DSetDesc;


/****************************************************************************
 * Low-level helpers
 */

static int get_ans_type(hid_t dset_id, int as_int, SEXPTYPE *ans_type)
{
	hid_t type_id;
	H5T_class_t class;
	size_t size;
	const char *classname;

	if (as_int) {
		*ans_type = INTSXP;
		return 0;
	}
	type_id = H5Dget_type(dset_id);
	if (type_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		return -1;
	}
	class = H5Tget_class(type_id);
	size = H5Tget_size(type_id);
	H5Tclose(type_id);
	if (class == H5T_NO_CLASS) {
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		return -1;
	}
	if (size == 0) {
		PRINT_TO_ERRMSG_BUF("H5Tget_size() returned 0");
		return -1;
	}

	if (class == H5T_INTEGER) {
		*ans_type = size <= sizeof(int) ? INTSXP : REALSXP;
		return 0;
	}
	if (class == H5T_FLOAT) {
		*ans_type = REALSXP;
		return 0;
	}
	switch (class) {
		case H5T_TIME: classname = "H5T_TIME"; break;
		case H5T_STRING: classname = "H5T_STRING"; break;
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

static int get_dset_desc(hid_t dset_id, int ndim, DSetDesc *dset_desc)
{
	hid_t space_id, plist_id;
	int dset_ndim;
	hsize_t *h5dim, *h5chunkdim;

	space_id = H5Dget_space(dset_id);
	if (space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		return -1;
	}

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

	plist_id = H5Dget_create_plist(dset_id);
	if (plist_id < 0) {
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		return -1;
	}

	/* Allocate 'h5dim' and 'h5chunkdim'. */
	h5dim = (hsize_t *) malloc(2 * ndim * sizeof(hsize_t));
	if (h5dim == NULL) {
		H5Pclose(plist_id);
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'h5dim' and 'h5chunkdim'");
		return -1;
	}
	h5chunkdim = h5dim + ndim;

	/* Set 'h5dim'. */
	if (H5Sget_simple_extent_dims(space_id, h5dim, NULL) != ndim) {
		free(h5dim);
		H5Pclose(plist_id);
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		return -1;
	}

	/* Set 'h5chunkdim'. */
	if (H5Pget_chunk(plist_id, ndim, h5chunkdim) != ndim) {
		free(h5dim);
		H5Pclose(plist_id);
		H5Sclose(space_id);
		PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
				    "an unexpected value");
		return -1;
	}

	dset_desc->dset_id = dset_id;
	dset_desc->space_id = space_id;
	dset_desc->plist_id = plist_id;
	dset_desc->ndim = ndim;
	dset_desc->h5dim = h5dim;
	dset_desc->h5chunkdim = h5chunkdim;
	return 0;
}

static void destroy_dset_desc(DSetDesc *dset_desc)
{
	free(dset_desc->h5dim);
	H5Pclose(dset_desc->plist_id);
	H5Sclose(dset_desc->space_id);
}

static int check_selection(const DSetDesc *dset_desc,
			   SEXP starts, SEXP counts,
			   int *nstart, int *ans_dim,
			   int *nblock, long long int *last_block_start)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = dset_desc->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) dset_desc->h5dim[h5along];
	return _deep_check_selection(starts, counts, dim_buf->elts,
				     nstart, ans_dim,
				     nblock, last_block_start);
}

static int map_starts_to_chunks(const DSetDesc *dset_desc,
				SEXP starts,
				int *ans_dim,
				IntAEAE *breakpoint_bufs,
				IntAEAE *chunkidx_bufs)
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
			(long long int) dset_desc->h5chunkdim[h5along];
	}
	return _map_starts_to_chunks(starts, dim_buf->elts, chunkdim_buf->elts,
				     ans_dim,
				     breakpoint_bufs, chunkidx_bufs);
}

static hid_t get_mem_space(int ndim, const int *ans_dim)
{
	hsize_t *h5dim;
	int along, h5along;
	hid_t mem_space_id;

	/* Allocate and set 'h5dim'. */
	h5dim = (hsize_t *) malloc(ndim * sizeof(hsize_t));
	if (h5dim == NULL) {
		PRINT_TO_ERRMSG_BUF("failed to allocate memory for 'h5dim'");
		return -1;
	}
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		h5dim[h5along] = ans_dim[along];

	mem_space_id = H5Screate_simple(ndim, h5dim, NULL);
	if (mem_space_id < 0)
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
	free(h5dim);
	return mem_space_id;
}

static int next_midx(int ndim, int *midx, const int *nstart)
{
	int along, i;

	for (along = 0; along < ndim; along++) {
		i = midx[along] + 1;
		if (i < nstart[along]) {
			midx[along] = i;
			break;
		}
		midx[along] = 0;
	}
	return along;
}


/****************************************************************************
 * read_data1()
 *
 * Some useful links:
 * - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5S.html
 * - Documentation of H5Dread():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
 * - An H5Dread() example:
 *     https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
 */

static int select_hyperslab(hid_t file_space_id, int ndim,
			    SEXP starts, SEXP counts,
			    const int *midx, int moved_along,
			    hsize_t *h5offset_buf, hsize_t *h5count_buf)
{
	int along, h5along, i, ret;
	SEXP start, count;

	/* Set 'h5offset_buf' and 'h5count_buf'. */
	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		h5offset_buf[h5along] = _get_trusted_elt(start, i) - 1;
		if (counts != R_NilValue) {
			count = VECTOR_ELT(counts, along);
			if (count != R_NilValue)
				h5count_buf[h5along] =
					_get_trusted_elt(count, i);
		}
	}

	/* Add to current selection. */
	ret = H5Sselect_hyperslab(file_space_id, H5S_SELECT_OR,
				  h5offset_buf, NULL, h5count_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_hyperslabs(const DSetDesc *dset_desc,
			SEXP starts, SEXP counts, const int *nstart)
{
	int ndim, along, h5along, moved_along, ret;
	long long int num_hyperslabs;
	hsize_t *h5offset_buf, *h5count_buf;  // hyperslab offsets and dims
	IntAE *midx_buf;

	ndim = dset_desc->ndim;
	for (along = 0; along < ndim; along++)
		if (nstart[along] == 0)
			return 0;  // empty region (no hyperslab)

	ret = H5Sselect_none(dset_desc->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	/* Allocate 'h5offset_buf' and 'h5count_buf'. */
	h5offset_buf = (hsize_t *) malloc(2 * ndim * sizeof(hsize_t));
	if (h5offset_buf == NULL) {
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'h5offset_buf' and 'h5count_buf'");
		return -1;
	}
	h5count_buf = h5offset_buf + ndim;

	/* Initialize 'h5offset_buf' and 'h5count_buf'. */
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			h5offset_buf[h5along] = 0;
			h5count_buf[h5along] = dset_desc->h5dim[h5along];
		} else {
			h5count_buf[h5along] = 1;
		}
	}

	/* Walk on the blocks. */
	midx_buf = new_IntAE(ndim, ndim, 0);
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		ret = select_hyperslab(dset_desc->space_id, ndim,
				       starts, counts,
				       midx_buf->elts, moved_along,
				       h5offset_buf, h5count_buf);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, midx_buf->elts, nstart);
	} while (moved_along < ndim);

	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(nstart)
	free(h5offset_buf);
	return num_hyperslabs;
}

/* Return nb of elements (or -1 on error). */
static long long int select_elements(const DSetDesc *dset_desc,
			SEXP starts, const int *ans_dim)
{
	int ndim, along, moved_along, i, ret;
	size_t num_elements;
	hsize_t *coord_buf, *coord_p;
	IntAE *midx_buf;
	SEXP start;
	long long int coord;

	ndim = dset_desc->ndim;
	num_elements = 1;
	for (along = 0; along < ndim; along++)
		num_elements *= ans_dim[along];
	//printf("nb of elements = %lu\n", num_elements);  // = length(ans)
	if (num_elements == 0)
		return 0;  // no elements

	/* Allocate 'coord_buf'. */
	coord_buf = (hsize_t *) malloc(num_elements * ndim * sizeof(hsize_t));
	if (coord_buf == NULL) {
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'coord_buf'");
		return -1;
	}

	/* Walk on the elements. */
	midx_buf = new_IntAE(ndim, ndim, 0);
	coord_p = coord_buf;
	moved_along = ndim;
	do {
		for (along = ndim - 1; along >= 0; along--) {
			i = midx_buf->elts[along];
			start = VECTOR_ELT(starts, along);
			if (start == R_NilValue) {
				coord = i;
			} else {
				coord = _get_trusted_elt(start, i) - 1;
			}
			*(coord_p++) = (hsize_t) coord;
		}
		moved_along = next_midx(ndim, midx_buf->elts, ans_dim);
	} while (moved_along < ndim);

	ret = H5Sselect_elements(dset_desc->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	free(coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

static int set_selection(const DSetDesc *dset_desc,
			 SEXP starts, SEXP counts, int noreduce,
			 const int *nstart,
			 const int *ans_dim,
			 const int *nblock,
			 const long long int *last_block_start)
{
	int ndim;
	SEXP reduced;
	long long int nselection;

	ndim = dset_desc->ndim;
	if (!noreduce && _selection_can_be_reduced(ndim, nstart, nblock)) {
		reduced = PROTECT(_reduce_selection(starts, counts,
						    ans_dim,
						    nblock, last_block_start));
		starts = VECTOR_ELT(reduced, 0);
		counts = VECTOR_ELT(reduced, 1);
		nstart = nblock;
	}

	//clock_t t0 = clock();
	nselection = counts != R_NilValue ?
		     select_hyperslabs(dset_desc, starts, counts, nstart) :
		     select_elements(dset_desc, starts, ans_dim);
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("time for setting selection: %e\n", dt);
	//printf("nselection: %lld, time per selection unit: %e\n",
	//	nselection, dt / nselection);

	if (nstart == nblock)
		UNPROTECT(1);  // unprotect 'reduced'

	return nselection < 0 ? -1 : 0;
}

static int read_data1(const DSetDesc *dset_desc,
		      SEXP starts, SEXP counts, int noreduce,
		      const int *nstart,
		      const int *ans_dim,
		      const int *nblock,
		      const long long int *last_block_start,
		      void *out, hid_t mem_space_id, hid_t mem_type_id)
{
	int ret;

	ret = set_selection(dset_desc, starts, counts, noreduce,
			    nstart, ans_dim,
			    nblock, last_block_start);
	if (ret < 0)
		return -1;

	ret = H5Sselect_all(mem_space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_all() returned an error");
		return -1;
	}

	//clock_t t0 = clock();
	ret = H5Dread(dset_desc->dset_id, mem_type_id,
		      mem_space_id, dset_desc->space_id,
		      H5P_DEFAULT, out);
	//printf("time for reading data from selection: %e\n",
	//	(1.0 * clock() - t0) / CLOCKS_PER_SEC);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}


/****************************************************************************
 * read_data2()
 */

static int read_selection_unit(const DSetDesc *dset_desc,
		SEXP starts, SEXP counts,
		const int *midx, int moved_along,
		hsize_t *h5offset_in_buf, hsize_t *h5count_buf,
		hsize_t *h5offset_out_buf,
		void *out, hid_t mem_space_id, hid_t mem_type_id)
{
	int ndim, along, h5along, i, ret;
	SEXP start, count;

	ndim = dset_desc->ndim;

	/* Set 'h5offset_in_buf', 'h5count_buf', and 'h5offset_out_buf'. */
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

	return H5Dread(dset_desc->dset_id, mem_type_id,
		       mem_space_id, dset_desc->space_id,
		       H5P_DEFAULT, out);
}

static int read_data2(const DSetDesc *dset_desc,
		      SEXP starts, SEXP counts, int noreduce,
		      const int *nstart,
		      const int *ans_dim,
		      const int *nblock,
		      const long long int *last_block_start,
		      void *out, hid_t mem_space_id, hid_t mem_type_id)
{
	int ndim, along, h5along, moved_along, ret;
	hsize_t *h5offset_in_buf, *h5count_buf, *h5offset_out_buf;
	IntAE *midx_buf;
	long long int num_selection_units;

	ndim = dset_desc->ndim;

	/* Allocate 'h5offset_in_buf', 'h5count_buf', and 'h5offset_out_buf'. */
	h5offset_in_buf = (hsize_t *) malloc(3 * ndim * sizeof(hsize_t));
	if (h5offset_in_buf == NULL) {
		PRINT_TO_ERRMSG_BUF(
			"failed to allocate memory "
			"for 'h5offset_in_buf', 'h5count_buf', "
			"and 'h5offset_out_buf'");
		return -1;
	}
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

	/* Walk on the selection units. */
	midx_buf = new_IntAE(ndim, ndim, 0);
	num_selection_units = 0;
	moved_along = ndim;
	do {
		num_selection_units++;
		ret = read_selection_unit(dset_desc,
					  starts, counts,
					  midx_buf->elts, moved_along,
					  h5offset_in_buf, h5count_buf,
					  h5offset_out_buf,
					  out, mem_space_id, mem_type_id);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, midx_buf->elts, nstart);
	} while (moved_along < ndim);

	//printf("nb of selection units = %lld\n", num_selection_units); // = prod(nstart)
	free(h5offset_in_buf);
	return ret;
}


/****************************************************************************
 * read_data3()
 */

static void set_offset_and_count_bufs(const DSetDesc *dset_desc,
				      int moved_along,
				      hsize_t *h5offset_buf,
				      hsize_t *h5count_buf)
{
	int ndim, along, h5along;
	hsize_t h5chunk_d;

	ndim = dset_desc->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		h5chunk_d = dset_desc->h5chunkdim[h5along];
		if (along == moved_along) {
			h5offset_buf[h5along] += h5chunk_d;
		} else {
			h5offset_buf[h5along] = 0;
		}
		h5count_buf[h5along] = dset_desc->h5dim[h5along] -
				       h5offset_buf[h5along];
		if (h5count_buf[h5along] > h5chunk_d)
			h5count_buf[h5along] = h5chunk_d;
	}
	return;
}

/* It takes about 218s on my laptop to load all the chunks from the EH1040
 * dataset (big 10x Genomics brain dataset in dense format, chunks of 100x100,
 * wrapped in the TENxBrainData package). That's 60 microseconds per chunk! */
static int read_chunk(const DSetDesc *dset_desc,
		      const hsize_t *h5offset_buf, const hsize_t *h5count_buf,
		      const hsize_t *h5zeros_buf,
		      void *chunkout, hid_t mem_space_id, hid_t mem_type_id)
{
	int ret;

	ret = H5Sselect_hyperslab(dset_desc->space_id, H5S_SELECT_SET,
				  h5offset_buf, NULL, h5count_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	ret = H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET,
				  h5zeros_buf, NULL, h5count_buf, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}

	return H5Dread(dset_desc->dset_id, mem_type_id,
		       mem_space_id, dset_desc->space_id,
		       H5P_DEFAULT, chunkout);
}

/*
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

static void test_zlib_compress2_uncompress_cycle(int level)
{
	const int chunk_len = 100000;
	int chunk_data[chunk_len], buf[chunk_len], buf2[chunk_len];
	uLongf chunk_size = sizeof(chunk_data), compressed_size;
	int i, ret;

	for (i = 0; i < chunk_len; i++)
		chunk_data[i] = 99999999 - 777 * i;
	compressed_size = sizeof(buf);
	printf("ok1\n");
	ret = compress2((Bytef *) buf, &compressed_size,
			(Bytef *) chunk_data, chunk_size, level);
	printf("ok2\n");
	if (ret != 0)
		error("compress2() returned an error");
	printf("size of compressed data = %lu\n", compressed_size);
	ret = uncompress_chunk_data(buf, (hsize_t) compressed_size,
				    buf2, chunk_size);
	if (ret != 0)
		error(_HDF5Array_errmsg_buf);
	if (memcmp(chunk_data, buf2, chunk_size) != 0)
		error("data differs");
	return;
}
*/

/*
 * For some unclear reason, H5Dread_chunk() is NOT listed here:
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
 *
static int direct_read_chunk(const DSetDesc *dset_desc,
			     const hsize_t *h5offset_buf,
			     void *chunkout, void *compressed_chunk_data,
			     size_t chunk_size)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(dset_desc->dset_id,
					h5offset_buf, &chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > chunk_size) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() returned "
				    "an unexpected chunk storage size (%llu)",
				    chunk_storage_size);
		return -1;
	}
	if (chunk_storage_size == chunk_size)
		compressed_chunk_data = chunkout;
	ret = H5Dread_chunk(dset_desc->dset_id, H5P_DEFAULT,
			    h5offset_buf, &filters, compressed_chunk_data);
	if (ret < 0)
		return -1;
	if (chunk_storage_size == chunk_size)
		return 0;
	return uncompress_chunk_data(compressed_chunk_data, chunk_storage_size,
				     chunkout, chunk_size);
}
*/

static int read_data3(const DSetDesc *dset_desc,
		      SEXP starts,
		      const int *ans_dim,
		      const IntAEAE *breakpoint_bufs,
		      const IntAEAE *chunkidx_bufs,
		      void *out, hid_t mem_space_id, hid_t mem_type_id)
{
	int ndim, along, h5along, moved_along, ret;
	size_t chunk_len, chunk_size;
	void *chunkout; // *chunkout2, *compressed_chunk_data;
	hsize_t *h5offset_buf, *h5count_buf, *h5zeros_buf;
	IntAE *midx_buf, *nchunk_buf;
	long long int num_chunks, ndot;

	ndim = dset_desc->ndim;
	mem_space_id = H5Screate_simple(ndim, dset_desc->h5chunkdim, NULL);
	if (mem_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}

	chunk_len = 1;
	for (along = 0; along < ndim; along++)
		chunk_len *= dset_desc->h5chunkdim[along];
	printf("chunk_len = %lu\n", chunk_len);

	if (mem_type_id == H5T_NATIVE_INT) {
		chunk_size = chunk_len * sizeof(int);
        } else {
		chunk_size = chunk_len * sizeof(double);
	}
	chunkout = malloc(3 * chunk_size);
	if (chunkout == NULL)
		return -1;
	//chunkout2 = chunkout + chunk_len;
	//compressed_chunk_data = chunkout2 + chunk_len;

        /* Allocate 'h5offset_buf', 'h5count_buf', and 'h5zeros_buf'. */
        h5offset_buf = (hsize_t *) malloc(3 * ndim * sizeof(hsize_t));
        if (h5offset_buf == NULL) {
		free(chunkout);
                PRINT_TO_ERRMSG_BUF("failed to allocate memory "
                                    "for 'h5offset_buf', 'h5count_buf', "
				    "and 'h5zeros_buf'");
                return -1;
        }
        h5count_buf = h5offset_buf + ndim;
        h5zeros_buf = h5count_buf + ndim;

	/* Walk on the chunks. */
	midx_buf = new_IntAE(ndim, ndim, 0);
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	printf("nchunk_buf:");
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		h5zeros_buf[h5along] = 0;
		nchunk_buf->elts[along] =
		    dset_desc->h5dim[h5along] / dset_desc->h5chunkdim[h5along];
		if (dset_desc->h5dim[h5along] % dset_desc->h5chunkdim[h5along])
			nchunk_buf->elts[along]++;
		printf(" %d", nchunk_buf->elts[along]);
	}
	printf("\n");

	num_chunks = ndot = 0;
	moved_along = ndim;
	do {
		num_chunks++;
		if (num_chunks % 20000 == 0) {
			printf(".");
			if (++ndot % 50 == 0)
				printf(" [%lld chunks processed]\n",
				       num_chunks);
			fflush(stdout);
		}
		set_offset_and_count_bufs(dset_desc,
				 moved_along, h5offset_buf, h5count_buf);

		ret = read_chunk(dset_desc,
				 h5offset_buf, h5count_buf, h5zeros_buf,
				 chunkout, mem_space_id, mem_type_id);
		if (ret < 0)
			break;
/*
		// Doesn't work yet! (For some obscure reason, trying to
		// decompress the chunk data fails with uncompress() returns
		// Z_DATA_ERROR (data corrupted or incomplete).
		ret = direct_read_chunk(dset_desc, h5offset_buf,
					chunkout2, compressed_chunk_data,
					chunk_size);
		if (ret < 0)
			break;

		if (memcmp(chunkout, chunkout2, chunk_len) != 0) {
			printf("chunkout:\n");
			for (size_t i = 0; i < 10; i++) {
			  for (size_t j = 0; j < 10; j++)
			    printf(" %4d", ((int *) chunkout)[i + 100 * j]);
			  printf("\n");
			}
			printf("chunkout2:\n");
			for (size_t i = 0; i < 10; i++) {
			  for (size_t j = 0; j < 10; j++)
			    printf(" %4d", ((int *) chunkout2)[i + 100 * j]);
			  printf("\n");
			}
			PRINT_TO_ERRMSG_BUF("chunkout != chunkout2");
			ret = -1;
			break;
		}
*/
		moved_along = next_midx(ndim, midx_buf->elts, nchunk_buf->elts);
	} while (moved_along < ndim);
	printf("\n");

	printf("nb of chunks = %lld\n", num_chunks);

	free(h5offset_buf);
	free(chunkout);
	return ret;
}


/****************************************************************************
 * C_h5mread()
 */

/* Return R_NilValue on error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce,
		    int method, int as_int)
{
	int ret, ndim, along;
	SEXPTYPE ans_type;
	DSetDesc dset_desc;
	SEXP ans_dim, ans;

	IntAE *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	int *nstart, *nblock;
	long long int *last_block_start;

	IntAEAE *breakpoint_bufs, *chunkidx_bufs;

	R_xlen_t ans_len;
	void *out;
	hid_t mem_space_id, mem_type_id;

	ret = get_ans_type(dset_id, as_int, &ans_type);
	if (ret < 0)
		return R_NilValue;

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0)
		return R_NilValue;

	ret = get_dset_desc(dset_id, ndim, &dset_desc);
	if (ret < 0)
		return R_NilValue;

	ans_dim = PROTECT(NEW_INTEGER(ndim));

	if (method != 3) {
		nstart_buf = new_IntAE(ndim, ndim, 0);
		nblock_buf = new_IntAE(ndim, ndim, 0);
		last_block_start_buf = new_LLongAE(ndim, ndim, 0);

		nstart = nstart_buf->elts;
		nblock = nblock_buf->elts;
		last_block_start = last_block_start_buf->elts;

		/* This call will populate 'nstart', 'ans_dim', 'nblock',
		   and 'last_block_start'. */
		ret = check_selection(&dset_desc, starts, counts,
				      nstart, INTEGER(ans_dim),
				      nblock, last_block_start);
	} else {
		breakpoint_bufs = new_IntAEAE(ndim, ndim);
		chunkidx_bufs = new_IntAEAE(ndim, ndim);
		/* This call will populate 'ans_dim', 'breakpoint_bufs',
		   and 'chunkidx_bufs'. */
		ret = map_starts_to_chunks(&dset_desc, starts,
				      INTEGER(ans_dim),
				      breakpoint_bufs, chunkidx_bufs);
	}
	if (ret < 0)
		goto on_error1;

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= INTEGER(ans_dim)[along];
	ans = PROTECT(allocVector(ans_type, ans_len));
	SET_DIM(ans, ans_dim);

	if (ans_type == INTSXP) {
		out = INTEGER(ans);
		mem_type_id = H5T_NATIVE_INT;
	} else {
		out = REAL(ans);
		mem_type_id = H5T_NATIVE_DOUBLE;
	}

	mem_space_id = get_mem_space(dset_desc.ndim, INTEGER(ans_dim));
	if (mem_space_id < 0)
		goto on_error2;

	if (ans_len != 0) {
		if (method == 1) {
			ret = read_data1(&dset_desc,
				starts, counts, noreduce,
				nstart, INTEGER(ans_dim),
				nblock, last_block_start,
				out, mem_space_id, mem_type_id);
		} else if (method == 2) {
			ret = read_data2(&dset_desc,
				starts, counts, noreduce,
				nstart, INTEGER(ans_dim),
				nblock, last_block_start,
				out, mem_space_id, mem_type_id);
		} else if (method == 3) {
			ret = read_data3(&dset_desc,
				starts, INTEGER(ans_dim),
				breakpoint_bufs, chunkidx_bufs,
				out, mem_space_id, mem_type_id);
		} else if (method != 0) {
			PRINT_TO_ERRMSG_BUF("'method' must be 0, 1, 2, or 3");
			ret = -1;
		}
		if (ret < 0)
			goto on_error3;
	}

	H5Sclose(mem_space_id);
	destroy_dset_desc(&dset_desc);
	UNPROTECT(2);
	return ans;

    on_error3:
	H5Sclose(mem_space_id);
    on_error2:
	UNPROTECT(1);
    on_error1:
	UNPROTECT(1);
	destroy_dset_desc(&dset_desc);
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce,
	       SEXP method, SEXP as_integer)
{
	SEXP filepath0, name0, ans;
	int noreduce0, method0, as_int;
	hid_t file_id, dset_id;

	//test_zlib_compress2_uncompress_cycle(9);

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

	/* Check 'method'. */
	if (!(IS_INTEGER(method) && LENGTH(method) == 1))
		error("'method' must be a single integer");
	method0 = INTEGER(method)[0];

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

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
			      method0, as_int));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ans == R_NilValue) {
		UNPROTECT(1);
		error(_HDF5Array_errmsg_buf);
	}
	UNPROTECT(1);
	return ans;
}

