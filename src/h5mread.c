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

#include <zlib.h>

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */

//#include <time.h>

typedef struct {
	hsize_t *h5off, *h5dim;
	int *off, *dim; // same as h5off and h5dim but stored as int
} Viewport;


/****************************************************************************
 * Low-level helpers
 */

static void *get_dataptr(SEXP x)
{
	SEXPTYPE Rtype;

	Rtype = TYPEOF(x);
	switch (Rtype) {
	    case LGLSXP:  return LOGICAL(x);
	    case INTSXP:  return INTEGER(x);
	    case REALSXP: return REAL(x);
	    case RAWSXP:  return RAW(x);
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return NULL;
}

/* 'mode' controls what fields we should allocate memory for:
 *     mode == 0: 'h5off' and 'h5dim' fields only
 *     mode == 1: 'off' and 'dim' fields only
 *     mode == 2: for all the fields
 */
static int alloc_Viewport(Viewport *vp, int ndim, int mode)
{
	vp->h5off = NULL;
	vp->off = NULL;
	if (mode != 1) {
		/* Allocate memory for the 'h5off' and 'h5dim' fields. */
		vp->h5off = _alloc_hsize_t_buf(2 * ndim, 0, "Viewport fields");
		if (vp->h5off == NULL)
			return -1;
		vp->h5dim = vp->h5off + ndim;
	}
	if (mode != 0) {
		/* Allocate memory for the 'off' and 'dim' fields. */
		vp->off = (int *) malloc(2 * ndim * sizeof(int));
		if (vp->off == NULL) {
			if (mode != 1)
				free(vp->h5off);
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for Viewport fields");
			return -1;
		}
		vp->dim = vp->off + ndim;
	}
	return 0;
}

static void free_Viewport(Viewport *vp)
{
	if (vp->h5off != NULL)
		free(vp->h5off);
	if (vp->off != NULL)
		free(vp->off);
	return;
}

/* Used in read_data_4_5() and read_data_7(). */
static int alloc_chunkvp_middlevp_destvp_bufs(int ndim,
			Viewport *chunkvp_buf,
			Viewport *middlevp_buf,
			Viewport *destvp_buf, int destvp_mode)
{
	if (alloc_Viewport(chunkvp_buf, ndim, 0) < 0)
		return -1;
	middlevp_buf->h5off =
		_alloc_hsize_t_buf(ndim, 1, "'middlevp_buf->h5off'");
	middlevp_buf->h5dim = chunkvp_buf->h5dim;
	if (middlevp_buf->h5off == NULL) {
		free_Viewport(chunkvp_buf);
		return -1;
	}
	if (alloc_Viewport(destvp_buf, ndim, destvp_mode) < 0) {
		free(middlevp_buf->h5off);
		free_Viewport(chunkvp_buf);
		return -1;
	}
	return 0;
}

static void free_chunkvp_middlevp_destvp_bufs(
			Viewport *chunkvp_buf,
			Viewport *middlevp_buf,
			Viewport *destvp_buf)
{
	free_Viewport(destvp_buf);
	free(middlevp_buf->h5off);
	free_Viewport(chunkvp_buf);
	return;
}

/* Used in read_data_6(). */
static int alloc_chunkvp_innervp_destvp_bufs(int ndim,
			Viewport *chunkvp_buf,
			Viewport *innervp_buf,
			Viewport *destvp_buf)
{
	if (alloc_Viewport(chunkvp_buf, ndim, 0) < 0)
		return -1;
	if (alloc_Viewport(innervp_buf, ndim, 0) < 0) {
		free_Viewport(chunkvp_buf);
		return -1;
	}
	if (alloc_Viewport(destvp_buf, ndim, 2) < 0) {
		free_Viewport(innervp_buf);
		free_Viewport(chunkvp_buf);
		return -1;
	}
	return 0;
}

static void free_chunkvp_innervp_destvp_bufs(
			Viewport *chunkvp_buf,
			Viewport *innervp_buf,
			Viewport *destvp_buf)
{
	free_Viewport(destvp_buf);
	free_Viewport(innervp_buf);
	free_Viewport(chunkvp_buf);
	return;
}

static int set_hyperslab(hid_t space_id, const Viewport *vp)
{
	int ret;

	ret = H5Sselect_hyperslab(space_id, H5S_SELECT_SET,
				  vp->h5off, NULL, vp->h5dim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

static int add_hyperslab(hid_t space_id, const Viewport *vp)
{
	int ret;

	ret = H5Sselect_hyperslab(space_id, H5S_SELECT_OR,
				  vp->h5off, NULL, vp->h5dim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

static hid_t get_mem_space(int ndim, const int *ans_dim)
{
	hsize_t *h5dim;
	int along, h5along;
	hid_t mem_space_id;

	/* Allocate and set 'h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
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

static int read_Viewport(const H5DSetDescriptor *h5dset,
			 const Viewport *dsetvp,
			 const Viewport *memvp,
			 void *mem, hid_t mem_space_id)
{
	int ret;

	ret = set_hyperslab(h5dset->space_id, dsetvp);
	if (ret < 0)
		return -1;
	ret = set_hyperslab(mem_space_id, memvp);
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static long long int check_selection_against_h5dset(
		const H5DSetDescriptor *h5dset,
		SEXP starts, SEXP counts, int *selection_dim_buf)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = h5dset->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) h5dset->h5dim[h5along];
	return _check_selection(ndim, dim_buf->elts, starts, counts,
				selection_dim_buf);
}

static long long int check_ordered_selection_against_h5dset(
		const H5DSetDescriptor *h5dset,
		SEXP starts, SEXP counts, int *selection_dim_buf,
		int *nstart_buf, int *nblock_buf,
		long long int *last_block_start_buf)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = h5dset->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) h5dset->h5dim[h5along];
	return _check_ordered_selection(ndim, dim_buf->elts, starts, counts,
					selection_dim_buf,
					nstart_buf, nblock_buf,
					last_block_start_buf);
}

static int map_starts_to_chunks(const H5DSetDescriptor *h5dset,
				SEXP starts,
				int *nstart_buf,
				IntAEAE *breakpoint_bufs,
				LLongAEAE *chunkidx_bufs)
{
	int ndim, along, h5along;
	LLongAE *dim_buf, *chunkdim_buf;

	ndim = h5dset->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	chunkdim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		dim_buf->elts[along] =
			(long long int) h5dset->h5dim[h5along];
		chunkdim_buf->elts[along] =
			(long long int) h5dset->h5chunkdim[h5along];
	}
	return _map_starts_to_chunks(ndim, dim_buf->elts, chunkdim_buf->elts,
				     starts,
				     nstart_buf,
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
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			nblock = LENGTH(start);
		} else {
			nblock = expand ? ans_dim[along] : 1;
		}
		num_blocks *= nblock_buf[along] = nblock;
	}
	return num_blocks;
}

static void init_srcvp_buf(const H5DSetDescriptor *h5dset,
			   SEXP starts, Viewport *srcvp_buf)
{
	int ndim, along, h5along;
	SEXP start;
	hsize_t d;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
		if (start == R_NilValue) {
			srcvp_buf->h5off[h5along] = 0;
			d = h5dset->h5dim[h5along];
		} else {
			d = 1;
		}
		srcvp_buf->h5dim[h5along] = d;
	}
	return;
}

static void update_srcvp_buf(int ndim,
			const int *midx, int moved_along,
			SEXP starts, SEXP counts,
			Viewport *srcvp_buf)
{
	int along, h5along, i;
	SEXP start, count;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		if (starts == R_NilValue)
			continue;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		srcvp_buf->h5off[h5along] = _get_trusted_elt(start, i) - 1;
		if (counts == R_NilValue)
			continue;
		count = VECTOR_ELT(counts, along);
		if (count == R_NilValue)
			continue;
		srcvp_buf->h5dim[h5along] = _get_trusted_elt(count, i);
	}
	return;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_hyperslabs(const H5DSetDescriptor *h5dset,
			SEXP starts, SEXP counts, const int *ans_dim,
			int *nblock_buf, int *midx_buf)
{
	int ret, ndim, moved_along;
	Viewport srcvp_buf;
	long long int num_hyperslabs;

	ret = H5Sselect_none(h5dset->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = h5dset->ndim;
	set_nblock_buf(ndim, starts, 0, ans_dim, nblock_buf);

	/* Allocate 'srcvp_buf'. */
	if (alloc_Viewport(&srcvp_buf, ndim, 0) < 0)
		return -1;

	init_srcvp_buf(h5dset, starts, &srcvp_buf);

	/* Walk on the hyperslabs. */
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		update_srcvp_buf(ndim, midx_buf, moved_along,
				 starts, counts, &srcvp_buf);
		/* Add to current selection. */
		ret = add_hyperslab(h5dset->space_id, &srcvp_buf);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nblock_buf, midx_buf);
	} while (moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);

	free_Viewport(&srcvp_buf);
	return ret < 0 ? -1 : num_hyperslabs;
}

static inline hsize_t *add_element(int ndim, const int *midx,
				   SEXP starts, hsize_t *coord_p)
{
	int h5along, i;
	SEXP start;
	long long int coord;

	for (h5along = ndim - 1; h5along >= 0; h5along--) {
		i = midx[h5along];
		start = GET_LIST_ELT(starts, h5along);
		if (start != R_NilValue) {
			coord = _get_trusted_elt(start, i) - 1;
		} else {
			coord = i;
		}
		*(coord_p++) = (hsize_t) coord;
	}
	return coord_p;
}

/* Return nb of selected elements (or -1 on error). */
static long long int select_elements(const H5DSetDescriptor *h5dset,
			SEXP starts, const int *ans_dim,
			int *nblock_buf, int *midx_buf)
{
	int ndim, outer_moved_along, ret;
	size_t num_elements;
	hsize_t *coord_buf, *coord_p;

	ndim = h5dset->ndim;
	num_elements = set_nblock_buf(ndim, starts, 1, ans_dim, nblock_buf);

	/* Allocate 'coord_buf'. */
	coord_buf = _alloc_hsize_t_buf(num_elements * ndim, 0, "'coord_buf'");
	if (coord_buf == NULL)
		return -1;

	/* Walk on the selected elements. */
	coord_p = coord_buf;
	outer_moved_along = ndim;
	do {
		coord_p = add_element(ndim, midx_buf, starts, coord_p);
		outer_moved_along = next_midx(ndim, nblock_buf, midx_buf);
	} while (outer_moved_along < ndim);

	ret = H5Sselect_elements(h5dset->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	free(coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

static int set_selection(const H5DSetDescriptor *h5dset, int method,
			 SEXP starts, SEXP counts, const int *ans_dim)
{
	int ndim;
	IntAE *nblock_buf, *midx_buf;
	long long int num_blocks;

	ndim = h5dset->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	//clock_t t0 = clock();
	if (method == 1) {
		num_blocks = select_hyperslabs(h5dset,
					starts, counts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	} else {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "'method' is set to 2");
			return -1;
		}
		num_blocks = select_elements(h5dset,
					starts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	}
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("time for setting selection: %e\n", dt);
	//printf("num_blocks: %lld, time per block: %e\n",
	//	num_blocks, dt / num_blocks);
	return num_blocks < 0 ? -1 : 0;
}

static int read_data_1_2(const H5DSetDescriptor *h5dset, int method,
			 SEXP starts, SEXP counts, const int *ans_dim,
			 void *mem, hid_t mem_space_id)
{
	int ret;

	ret = set_selection(h5dset, method,
			    starts, counts, ans_dim);
	if (ret < 0)
		return -1;

	ret = H5Sselect_all(mem_space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_all() returned an error");
		return -1;
	}

	//clock_t t0 = clock();
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
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
 */

static int read_hyperslab(const H5DSetDescriptor *h5dset,
		SEXP starts, SEXP counts,
		const int *midx, int moved_along,
		Viewport *srcvp_buf, Viewport *destvp_buf,
		void *mem, hid_t mem_space_id)
{
	int ndim, along, h5along, i, ret;
	SEXP start;

	ndim = h5dset->ndim;

	/* Update 'destvp_buf' and 'srcvp_buf' IN THAT ORDER! */
	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		start = GET_LIST_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		if (i == 0) {
			destvp_buf->h5off[h5along] = 0;
		} else {
			destvp_buf->h5off[h5along] += srcvp_buf->h5dim[h5along];
		}
	}
	update_srcvp_buf(ndim, midx, moved_along, starts, counts, srcvp_buf);

	ret = set_hyperslab(h5dset->space_id, srcvp_buf);
	if (ret < 0)
		return -1;
	ret = set_hyperslab(mem_space_id, destvp_buf);
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int read_data_3(const H5DSetDescriptor *h5dset,
		       SEXP starts, SEXP counts, const int *ans_dim,
		       void *mem, hid_t mem_space_id)
{
	int ndim, moved_along, ret;
	Viewport srcvp_buf, destvp_buf;
	IntAE *nblock_buf, *midx_buf;
	long long int num_hyperslabs;

	ndim = h5dset->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	set_nblock_buf(ndim, starts, 0, ans_dim, nblock_buf->elts);

	/* Allocate 'srcvp_buf' and 'destvp_buf'. */
	if (alloc_Viewport(&srcvp_buf, ndim, 0) < 0)
		return -1;
	destvp_buf.h5off = _alloc_hsize_t_buf(ndim, 1, "'destvp_buf.h5off'");
	destvp_buf.h5dim = srcvp_buf.h5dim;
	if (destvp_buf.h5off == NULL) {
		free_Viewport(&srcvp_buf);
		return -1;
	}

	/* Initialize 'srcvp_buf' (this also initializes 'destvp_buf.h5dim'). */
	init_srcvp_buf(h5dset, starts, &srcvp_buf);

	/* Walk on the hyperslabs. */
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		ret = read_hyperslab(h5dset, starts, counts,
				     midx_buf->elts, moved_along,
				     &srcvp_buf, &destvp_buf,
				     mem, mem_space_id);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nblock_buf->elts, midx_buf->elts);
	} while (moved_along < ndim);

	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);
	free(destvp_buf.h5off);
	free_Viewport(&srcvp_buf);
	return ret;
}


/****************************************************************************
 * read_data_4_5()
 *
 * One call to H5Dread() per chunk touched by the user selection.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Make one call to H5Dread() to load the **entire** chunk data to an
 *     intermediate buffer.
 *   - Copy the user-selected data from the intermediate buffer to 'ans'.
 *
 * Assumes that 'h5dset->h5chunkdim' and 'h5dset->h5nchunk' are NOT
 * NULL. This is NOT checked!
 */

static void set_nchunk_buf(const H5DSetDescriptor *h5dset,
			   const SEXP starts, const LLongAEAE *chunkidx_bufs,
			   IntAE *nchunk_buf)
{
	int ndim, along, h5along, nchunk;
	SEXP start;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			nchunk = LLongAE_get_nelt(chunkidx_bufs->elts[along]);
		} else {
			nchunk = h5dset->h5nchunk[h5along];
		}
		nchunk_buf->elts[along] = nchunk;
	}
	return;
}

static void update_chunkvp_buf(const H5DSetDescriptor *h5dset,
			const int *chunk_midx, int moved_along,
			SEXP starts, const LLongAEAE *chunkidx_bufs,
			Viewport *chunkvp_buf)
{
	int ndim, along, h5along, i;
	SEXP start;
	long long int chunkidx;
	hsize_t chunkd, off, d;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = chunk_midx[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			chunkidx = chunkidx_bufs->elts[along]->elts[i];
		} else {
			chunkidx = i;
		}
		chunkd = h5dset->h5chunkdim[h5along];
		off = chunkidx * chunkd;
		d = h5dset->h5dim[h5along] - off;
		if (d > chunkd)
			d = chunkd;
		chunkvp_buf->h5off[h5along] = off;
		chunkvp_buf->h5dim[h5along] = d;
	}
	//printf("# chunkvp_buf->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//	printf(" %llu", chunkvp_buf->h5off[h5along]);
	//printf("\n");
	//printf("# chunkvp_buf->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//	printf(" %llu", chunkvp_buf->h5dim[h5along]);
	//printf("\n");
	return;
}

static void update_destvp_buf(const H5DSetDescriptor *h5dset,
			const int *chunk_midx, int moved_along,
			SEXP starts, const IntAEAE *breakpoint_bufs,
			const Viewport *chunkvp, Viewport *destvp_buf)
{
	int ndim, along, h5along, i, off, d;
	SEXP start;
	const int *breakpoint;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = chunk_midx[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue ) {
			breakpoint = breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = chunkvp->h5off[h5along];
			d = chunkvp->h5dim[h5along];
		}
		if (destvp_buf->h5off != NULL) {
			destvp_buf->h5off[h5along] = off;
			destvp_buf->h5dim[h5along] = d;
		}
		destvp_buf->off[along] = off;
		destvp_buf->dim[along] = d;
	}
	//printf("# destvp_buf (offsets):");
	//for (along = 0; along < ndim; along++)
	//	printf(" %d", destvp_buf->off[along]);
	//printf("\n");
	//printf("# destvp_buf (dims):");
	//for (along = 0; along < ndim; along++)
	//	printf(" %d", destvp_buf->dim[along]);
	//printf("\n");
	return;
}

static void update_chunkvp_destvp_bufs(const H5DSetDescriptor *h5dset,
			const int *chunk_midx, int moved_along,
			SEXP starts,
			const IntAEAE *breakpoint_bufs,
			const LLongAEAE *chunkidx_bufs,
			Viewport *chunkvp_buf, Viewport *destvp_buf)
{
	update_chunkvp_buf(h5dset,
			chunk_midx, moved_along,
			starts, chunkidx_bufs,
			chunkvp_buf);
	update_destvp_buf(h5dset,
			chunk_midx, moved_along,
			starts, breakpoint_bufs,
			chunkvp_buf, destvp_buf);
	return;
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

static int load_chunk(const H5DSetDescriptor *h5dset,
		      const Viewport *chunkvp,
		      void *chunk_data_out,
		      void *compressed_chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(h5dset->dset_id,
					chunkvp->h5off,
					&chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > h5dset->chunk_data_buf_size +
				 COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%lu + %d)",
				    chunk_storage_size,
				    h5dset->chunk_data_buf_size,
				    COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(h5dset->dset_id, H5P_DEFAULT,
			    chunkvp->h5off, &filters,
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
				     h5dset->chunk_data_buf_size);
}

static void init_in_offset_and_out_offset(int ndim, SEXP starts,
			const int *out_dim, const Viewport *destvp,
			const Viewport *chunkvp,
			const hsize_t *h5chunkdim,
			size_t *in_offset, size_t *out_offset)
{
	size_t in_off, out_off;
	int along, h5along, i;
	SEXP start;

	in_off = out_off = 0;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunkdim[h5along];
		out_off *= out_dim[along];
		i = destvp->off[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  chunkvp->h5off[h5along];
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
			const int *out_dim, const Viewport *destvp,
			const hsize_t *h5chunkdim,
			size_t *in_offset, size_t *out_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc, out_off_inc;

	start = GET_LIST_ELT(starts, inner_moved_along);
	if (start != R_NilValue) {
		i1 = destvp->off[inner_moved_along] +
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
			in_off_inc *= h5chunkdim[h5along];
			out_off_inc *= out_dim[along];
			di = 1 - destvp->dim[along];
			start = GET_LIST_ELT(starts, along);
			if (start != R_NilValue) {
				i1 = destvp->off[along];
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

static int gather_chunk_data(const H5DSetDescriptor *h5dset,
			SEXP starts,
			const void *in,
			const Viewport *chunkvp,
			SEXP ans, const int *out_dim,
			const Viewport *destvp,
			int *inner_midx_buf)
{
	int ndim, is_na, inner_moved_along;
	size_t in_offset, out_offset, s_len;
	long long int num_elts;
	const char *s;
	SEXP ans_elt;

	ndim = h5dset->ndim;
	init_in_offset_and_out_offset(ndim, starts,
			out_dim, destvp,
			chunkvp, h5dset->h5chunkdim,
			&in_offset, &out_offset);

	/* Walk on the selected elements in current chunk. */
	num_elts = 0;
	while (1) {
		num_elts++;
		switch (h5dset->Rtype) {
		    case LGLSXP:
			LOGICAL(ans)[out_offset] = ((int *) in)[in_offset];
		    break;
		    case INTSXP:
			INTEGER(ans)[out_offset] = ((int *) in)[in_offset];
		    break;
		    case REALSXP:
			REAL(ans)[out_offset] = ((double *) in)[in_offset];
		    break;
		    case STRSXP:
			s = (char *) in + in_offset * h5dset->H5size;
			for (s_len = 0; s_len < h5dset->H5size; s_len++)
				if (s[s_len] == 0)
					break;
			is_na = h5dset->as_na_attr &&
				s_len == 2 && s[0] == 'N' && s[1] == 'A';
			if (is_na) {
				SET_STRING_ELT(ans, out_offset, NA_STRING);
			} else {
				ans_elt = PROTECT(mkCharLen(s, s_len));
				SET_STRING_ELT(ans, out_offset, ans_elt);
				UNPROTECT(1);
			}
		    break;
		    case RAWSXP:
			RAW(ans)[out_offset] = ((char *) in)[in_offset];
		    break;
		    default:
			PRINT_TO_ERRMSG_BUF("unsupported type: %s",
					    CHAR(type2str(h5dset->Rtype)));
			return -1;
		}
		inner_moved_along = next_midx(ndim, destvp->dim,
					      inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset_and_out_offset(ndim,
				inner_midx_buf, inner_moved_along,
				starts,
				out_dim, destvp,
				h5dset->h5chunkdim,
				&in_offset, &out_offset);
	};
	//printf("# nb of selected elements in current chunk = %lld\n",
	//       num_elts);
	return 0;
}

static int read_data_from_chunk_4_5(const H5DSetDescriptor *h5dset, int method,
			const int *chunk_midx, int moved_along,
			SEXP starts,
			const IntAEAE *breakpoint_bufs,
			const LLongAEAE *chunkidx_bufs,
			SEXP ans, const int *ans_dim,
			int *inner_midx_buf,
			Viewport *chunkvp_buf,
			const Viewport *middlevp,
			Viewport *destvp_buf,
			void *chunk_data_buf, hid_t chunk_space_id,
			void *compressed_chunk_data_buf)
{
	int ret;

	if (method == 4) {
		/* It takes about 218s on my laptop to load all the chunks
		   from the EH1040 dataset (big 10x Genomics brain dataset
		   in dense format, chunks of 100x100, wrapped in the
		   TENxBrainData package). That's 60 microseconds per chunk! */
		ret = read_Viewport(h5dset,
				chunkvp_buf, middlevp,
				chunk_data_buf, chunk_space_id);
	} else {
		ret = load_chunk(h5dset,
				chunkvp_buf,
				chunk_data_buf,
				compressed_chunk_data_buf);
	}
	if (ret < 0)
		return ret;
	ret = gather_chunk_data(h5dset,
			starts,
			chunk_data_buf,
			chunkvp_buf,
			ans, ans_dim, destvp_buf,
			inner_midx_buf);
	return ret;
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
static int read_data_4_5(const H5DSetDescriptor *h5dset, int method,
			 SEXP starts,
			 const IntAEAE *breakpoint_bufs,
			 const LLongAEAE *chunkidx_bufs,
			 SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ret;
	hid_t chunk_space_id;
	void *chunk_data_buf, *compressed_chunk_data_buf = NULL;
	Viewport chunkvp_buf, middlevp_buf, destvp_buf;
	IntAE *nchunk_buf, *chunk_midx_buf, *inner_midx_buf;
	long long int num_chunks, ndot;

	ndim = h5dset->ndim;
	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}

	/* Allocate 'chunkvp_buf', 'middlevp_buf', and 'destvp_buf'.
	   We set 'destvp_mode' to 1 because in the context of read_data_4_5()
	   we won't use 'destvp_buf.h5off' or 'destvp_buf.h5dim'. */
	if (alloc_chunkvp_middlevp_destvp_bufs(ndim,
		&chunkvp_buf, &middlevp_buf, &destvp_buf, 1) < 0)
	{
		H5Sclose(chunk_space_id);
		return -1;
	}

	if (method == 4) {
		chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	} else {
		warning("method 5 is still experimental, use at your own risk");
		chunk_data_buf = malloc(2 * h5dset->chunk_data_buf_size +
					COMPRESSION_OVERHEAD);
	}
	if (chunk_data_buf == NULL) {
		free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
						  &middlevp_buf,
						  &destvp_buf);
		H5Sclose(chunk_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}
	if (method != 4)
		compressed_chunk_data_buf = chunk_data_buf +
					    h5dset->chunk_data_buf_size;

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(h5dset, starts, chunkidx_bufs, nchunk_buf);
	chunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	/* Walk over the chunks touched by the user selection. */
	num_chunks = ndot = 0;
	moved_along = ndim;
	do {
		num_chunks++;
		//if (num_chunks % 20000 == 0) {
		//	printf(".");
		//	if (++ndot % 50 == 0)
		//		printf(" [%lld chunks processed]\n",
		//		       num_chunks);
		//	fflush(stdout);
		//}
		update_chunkvp_destvp_bufs(h5dset,
			chunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, chunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ret = read_data_from_chunk_4_5(h5dset, method,
			chunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, chunkidx_bufs,
			ans, ans_dim,
			inner_midx_buf->elts,
			&chunkvp_buf, &middlevp_buf, &destvp_buf,
			chunk_data_buf, chunk_space_id,
			compressed_chunk_data_buf);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nchunk_buf->elts,
					chunk_midx_buf->elts);
	} while (moved_along < ndim);
	free(chunk_data_buf);
	free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
					  &middlevp_buf,
					  &destvp_buf);
	H5Sclose(chunk_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to H5Dread() per chunk touched by the user selection.
 * No intermediate buffer.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Select the hyperslabs obtained by intersecting the selection with
 *     the current chunk.
 *   - Call H5Dread(). This loads the selected data **directly** to the
 *     output array.
 *
 * Assumes that 'h5dset->h5chunkdim' and 'h5dset->h5nchunk' are NOT
 * NULL. This is NOT checked!
 */

/*
static hsize_t *alloc_coord_buf(int ndim, const hsize_t *h5chunkdim)
{
	size_t max_inner_block_per_chunk;
	int along;

	max_inner_block_per_chunk = 1;
	for (along = 0; along < ndim; along++)
		max_inner_block_per_chunk *= (h5chunkdim[along] + 1) / 2;
	//printf("max_inner_block_per_chunk = %lu\n",
	//       max_inner_block_per_chunk);
	return _alloc_hsize_t_buf(max_inner_block_per_chunk * ndim,
				  "'coord_buf'");
}
*/

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP starts,
		const Viewport *destvp,
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
		d = destvp->dim[along];
		start = GET_LIST_ELT(starts, along);
		if (start == R_NilValue) {
			IntAE_insert_at(inner_breakpoint_buf, 0, d);
			inner_nblock_buf[along] = 1;
			continue;
		}
		off = destvp->off[along];
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

static void init_innervp_buf(int ndim,
			SEXP starts, const Viewport *chunkvp,
			Viewport *innervp_buf)
{
	int along, h5along;
	SEXP start;
	hsize_t d;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
		if (start == R_NilValue) {
			innervp_buf->h5off[h5along] =
					chunkvp->h5off[h5along];
			d = chunkvp->h5dim[h5along];
		} else {
			d = 1;
		}
		innervp_buf->h5dim[h5along] = d;
	}
	return;
}

static void update_innervp_buf(int ndim,
			SEXP starts, const Viewport *destvp,
			const int *inner_midx, int inner_moved_along,
			const IntAEAE *inner_breakpoint_bufs,
			Viewport *innervp_buf)
{
	int along, h5along, idx, off, d, i;
	SEXP start;
	const int *inner_breakpoint;

	for (along = 0; along < ndim; along++) {
		if (along > inner_moved_along)
			break;
		if (starts == R_NilValue)
			continue;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		inner_breakpoint = inner_breakpoint_bufs->elts[along]->elts;
		idx = inner_midx[along];
		off = idx == 0 ? 0 : inner_breakpoint[idx - 1];
		d = inner_breakpoint[idx] - off;
		i = destvp->off[along] + off;
		h5along = ndim - 1 - along;
		innervp_buf->h5off[h5along] = _get_trusted_elt(start, i) - 1;
		innervp_buf->h5dim[h5along] = d;
	}
	return;
}

/* Return nb of selected elements (or -1 on error). */
static long long int select_elements_from_chunk(
		const H5DSetDescriptor *h5dset,
		SEXP starts,
		const Viewport *destvp,
		int *inner_midx_buf, hsize_t *coord_buf)
{
	int ret, ndim, inner_moved_along, along, h5along, i;
	size_t num_elements;
	hsize_t *coord_p;
	long long int coord;
	SEXP start;

	ret = H5Sselect_none(h5dset->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = h5dset->ndim;

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
			i = destvp->h5off[h5along] + inner_midx_buf[along];
			start = GET_LIST_ELT(starts, along);
			if (start != R_NilValue) {
				coord = _get_trusted_elt(start, i) - 1;
			} else {
				coord = i;
			}
			*(coord_p++) = (hsize_t) coord;
		}
		inner_moved_along = next_midx(ndim, destvp->dim,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);

	ret = H5Sselect_elements(h5dset->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_intersection_of_blocks_with_chunk(
		const H5DSetDescriptor *h5dset,
		SEXP starts,
		const Viewport *destvp, const Viewport *chunkvp,
		const IntAEAE *inner_breakpoint_bufs,
		const int *inner_nblock,
		int *inner_midx_buf,
		Viewport *innervp_buf)
{
	int ret, ndim, inner_moved_along;
	long long int num_hyperslabs;

	ret = H5Sselect_none(h5dset->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = h5dset->ndim;

	init_innervp_buf(ndim, starts, chunkvp, innervp_buf);

	/* Walk on the block intersections with the currrent chunk. */
	num_hyperslabs = 0;
	inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_innervp_buf(ndim, starts, destvp,
				   inner_midx_buf, inner_moved_along,
				   inner_breakpoint_bufs,
				   innervp_buf);
		ret = add_hyperslab(h5dset->space_id, innervp_buf);
		if (ret < 0)
			return -1;
		inner_moved_along = next_midx(ndim, inner_nblock,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(inner_nblock)

	return num_hyperslabs;
}

static int read_selection(const H5DSetDescriptor *h5dset,
			  const Viewport *destvp_buf,
			  void *mem, hid_t mem_space_id)
{
	int ret;

	ret = set_hyperslab(mem_space_id, destvp_buf);
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int read_data_from_chunk_6(const H5DSetDescriptor *h5dset,
			const int *chunk_midx, int moved_along,
			SEXP starts,
			const IntAEAE *breakpoint_bufs,
			const LLongAEAE *chunkidx_bufs,
			void *mem, hid_t mem_space_id,
			int *inner_midx_buf,
			Viewport *chunkvp_buf,
			Viewport *innervp_buf,
			Viewport *destvp_buf,
			IntAEAE *inner_breakpoint_bufs,
			const IntAE *inner_nblock_buf)
			// hsize_t *coord_buf)
{
	int ret;

	update_inner_breakpoints(h5dset->ndim, moved_along,
			starts, destvp_buf,
			inner_breakpoint_bufs, inner_nblock_buf->elts);
	//t0 = clock();
	/* Having 'inner_nblock_buf->elts' identical to 'destvp_buf->dim'
	   means that all the inner blocks are single elements so we use
	   select_elements_from_chunk() which could be faster than
	   select_intersection_of_blocks_with_chunk() in that case.
	   NO IT'S NOT FASTER! */
	//ret = memcmp(inner_nblock_buf->elts, destvp_buf->dim,
	//	     ndim * sizeof(int));
	//printf("# chunk %lld: %s\n", num_chunks,
	//       ret == 0 ? "select_elements" : "select_hyperslabs");
	//if (ret == 0) {
	//	ret = select_elements_from_chunk(h5dset,
	//		starts, destvp_buf,
	//		inner_midx_buf, coord_buf);
	//} else {
		ret = select_intersection_of_blocks_with_chunk(
			h5dset, starts, destvp_buf, chunkvp_buf,
			inner_breakpoint_bufs, inner_nblock_buf->elts,
			inner_midx_buf, innervp_buf);
	//}
	if (ret < 0)
		return ret;
	//t_select_elements += clock() - t0;

	//t0 = clock();
	ret = read_selection(h5dset, destvp_buf,
			mem, mem_space_id);
	//t_read_selection += clock() - t0;
	return ret;
}

static int read_data_6(const H5DSetDescriptor *h5dset,
		       SEXP starts,
		       const IntAEAE *breakpoint_bufs,
		       const LLongAEAE *chunkidx_bufs,
		       SEXP ans, const int *ans_dim)
{
	void *mem;
	int ndim, moved_along, ret;
	hid_t mem_space_id;
	Viewport chunkvp_buf, innervp_buf, destvp_buf;
	//hsize_t *coord_buf;
	IntAE *nchunk_buf, *chunk_midx_buf,
	      *inner_nblock_buf, *inner_midx_buf;
	IntAEAE *inner_breakpoint_bufs;
	long long int num_chunks, ndot;

	mem = get_dataptr(ans);
	if (mem == NULL)
		return -1;
	ndim = h5dset->ndim;
	mem_space_id = get_mem_space(ndim, ans_dim);
	if (mem_space_id < 0)
		return -1;

	/* Allocate 'chunkvp_buf', 'innervp_buf', and 'destvp_buf'. */
	if (alloc_chunkvp_innervp_destvp_bufs(ndim,
		&chunkvp_buf, &innervp_buf, &destvp_buf) < 0)
	{
		H5Sclose(mem_space_id);
		return -1;
	}

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(h5dset, starts, chunkidx_bufs, nchunk_buf);
	chunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_nblock_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	/* Walk over the chunks touched by the user selection. */
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
		update_chunkvp_destvp_bufs(h5dset,
			chunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, chunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ret = read_data_from_chunk_6(h5dset,
			chunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, chunkidx_bufs,
			mem, mem_space_id,
			inner_midx_buf->elts,
			&chunkvp_buf, &innervp_buf, &destvp_buf,
			inner_breakpoint_bufs, inner_nblock_buf);
			//coord_buf);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nchunk_buf->elts,
					chunk_midx_buf->elts);
	} while (moved_along < ndim);
	//free(coord_buf);
	free_chunkvp_innervp_destvp_bufs(&chunkvp_buf,
					 &innervp_buf,
					 &destvp_buf);
	H5Sclose(mem_space_id);
	return ret;
}


/****************************************************************************
 * read_data_7()
 *
 * Like read_data_4_5() but bypasses the intermediate buffer if a chunk is
 * fully selected.
 *
 * Assumes that 'h5dset->h5chunkdim' and 'h5dset->h5nchunk' are NOT
 * NULL. This is NOT checked!
 */

static int chunk_is_fully_selected(int ndim,
				   const Viewport *chunkvp,
				   const Viewport *destvp)
{
	int ok, along, h5along;

	ok = 1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		//printf("chunkvp[%d]=%llu destvp[%d]=%llu\n",
		//       along, chunkvp->h5dim[h5along],
		//       along, destvp->h5dim[h5along]);
		if (chunkvp->h5dim[h5along] != destvp->h5dim[h5along]) {
			ok = 0;
			break;
		}
	}
	//printf("ok=%d\n", ok);
	return ok;
}

static int read_data_7(const H5DSetDescriptor *h5dset,
		       SEXP starts,
		       const IntAEAE *breakpoint_bufs,
		       const LLongAEAE *chunkidx_bufs,
		       SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ok, ret;
	hid_t chunk_space_id, mem_space_id;
	void *mem, *chunk_data_buf;
	Viewport chunkvp_buf, middlevp_buf, destvp_buf;
	IntAE *nchunk_buf, *chunk_midx_buf, *inner_midx_buf;
	long long int num_chunks;

	ndim = h5dset->ndim;
	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}
	mem = get_dataptr(ans);
	if (mem == NULL) {
		H5Sclose(chunk_space_id);
		return -1;
	}
	mem_space_id = get_mem_space(ndim, ans_dim);
	if (mem_space_id < 0) {
		H5Sclose(chunk_space_id);
		return -1;
	}

	/* Allocate 'chunkvp_buf', 'middlevp_buf', and 'destvp_buf'.
	   Unlike in read_data_4_5(), here we set 'destvp_mode' to 2 because
	   in the context of read_data_7() we will use 'destvp_buf.h5off'
	   and 'destvp_buf.h5dim'. */
	if (alloc_chunkvp_middlevp_destvp_bufs(ndim,
		&chunkvp_buf, &middlevp_buf, &destvp_buf, 2) < 0)
	{
		H5Sclose(mem_space_id);
		H5Sclose(chunk_space_id);
		return -1;
	}

	chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	if (chunk_data_buf == NULL) {
		free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
						  &middlevp_buf,
						  &destvp_buf);
		H5Sclose(mem_space_id);
		H5Sclose(chunk_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(h5dset, starts, chunkidx_bufs, nchunk_buf);
	chunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	/* Walk over the chunks touched by the user selection. */
	num_chunks = 0;
	moved_along = ndim;
	do {
		num_chunks++;
		update_chunkvp_destvp_bufs(h5dset,
			chunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, chunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ok = chunk_is_fully_selected(h5dset->ndim,
					     &chunkvp_buf, &destvp_buf);
		if (ok) {
			/* Load the chunk **directly** into 'ans' (no
			   intermediate buffer). */
			ret = read_Viewport(h5dset,
				&chunkvp_buf, &destvp_buf,
				mem, mem_space_id);
		} else {
			/* Load the **entire** chunk to an intermediate
			   buffer then copy the user-selected data from
			   the intermediate buffer to 'ans'. */
			ret = read_data_from_chunk_4_5(h5dset, 4,
				chunk_midx_buf->elts, moved_along,
				starts, breakpoint_bufs, chunkidx_bufs,
				ans, ans_dim,
				inner_midx_buf->elts,
				&chunkvp_buf, &middlevp_buf, &destvp_buf,
				chunk_data_buf, chunk_space_id, NULL);
		}
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nchunk_buf->elts,
					chunk_midx_buf->elts);
	} while (moved_along < ndim);
	free(chunk_data_buf);
	free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
					  &middlevp_buf,
					  &destvp_buf);
	H5Sclose(mem_space_id);
	H5Sclose(chunk_space_id);
	return ret;
}


/****************************************************************************
 * C_h5mread()
 */

/* Return R_NilValue on error. */
static SEXP h5mread_1_2_3(const H5DSetDescriptor *h5dset,
			  SEXP starts, SEXP counts, int noreduce,
			  int method, int *ans_dim)
{
	int ndim, ret;
	long long int ans_len;
	IntAE *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	SEXP ans, reduced;
	void *mem;
	hid_t mem_space_id;
	int nprotect = 0;

	ndim = h5dset->ndim;
	if (noreduce || method == 2) {
		/* This call will populate 'ans_dim'. */
		ans_len = check_selection_against_h5dset(h5dset,
					starts, counts, ans_dim);
		if (ans_len < 0)
			return R_NilValue;
	} else {
		nstart_buf = new_IntAE(ndim, ndim, 0);
		nblock_buf = new_IntAE(ndim, ndim, 0);
		last_block_start_buf = new_LLongAE(ndim, ndim, 0);
		/* This call will populate 'ans_dim', 'nblock_buf',
		   and 'last_block_start_buf'. */
		ans_len = check_ordered_selection_against_h5dset(h5dset,
					starts, counts, ans_dim,
					nstart_buf->elts, nblock_buf->elts,
					last_block_start_buf->elts);
		if (ans_len < 0)
			return R_NilValue;
		if (_selection_can_be_reduced(ndim, nstart_buf->elts,
						    nblock_buf->elts)) {
			reduced = PROTECT(_reduce_selection(
						ndim, starts, counts, ans_dim,
						nblock_buf->elts,
						last_block_start_buf->elts));
			nprotect++;
			starts = VECTOR_ELT(reduced, 0);
			counts = VECTOR_ELT(reduced, 1);
		}
	}

	ans = PROTECT(allocVector(h5dset->Rtype, (R_xlen_t) ans_len));
	nprotect++;

	if (ans_len != 0) {
		mem = get_dataptr(ans);
		if (mem == NULL)
			goto on_error;
		mem_space_id = get_mem_space(ndim, ans_dim);
		if (mem_space_id < 0)
			goto on_error;
		if (method <= 2) {
			ret = read_data_1_2(h5dset, method,
					starts, counts, ans_dim,
					mem, mem_space_id);
		} else {
			ret = read_data_3(h5dset,
					starts, counts, ans_dim,
					mem, mem_space_id);
		}
		H5Sclose(mem_space_id);
		if (ret < 0)
			goto on_error;
	}

	UNPROTECT(nprotect);
	return ans;

    on_error:
	UNPROTECT(nprotect);
	return R_NilValue;
}

/* Return R_NilValue on error. */
static SEXP h5mread_4_5_6_7(const H5DSetDescriptor *h5dset, SEXP starts,
			    int method, int *ans_dim)
{
	int ndim, ret, along;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *chunkidx_bufs;
	R_xlen_t ans_len;
	SEXP ans;

	ndim = h5dset->ndim;
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunkidx_bufs = new_LLongAEAE(ndim, ndim);

	/* This call will populate 'ans_dim', 'breakpoint_bufs',
	   and 'chunkidx_bufs'. */
	ret = map_starts_to_chunks(h5dset, starts,
				   ans_dim, breakpoint_bufs, chunkidx_bufs);
	if (ret < 0)
		return R_NilValue;

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(h5dset->Rtype, ans_len));

	if (ans_len != 0) {
		if (method <= 5) {
			/* methods 4 and 5 */
			ret = read_data_4_5(h5dset, method,
					starts, breakpoint_bufs, chunkidx_bufs,
					ans, ans_dim);
		} else if (method == 6) {
			/* method 6 */
			ret = read_data_6(h5dset,
					starts, breakpoint_bufs, chunkidx_bufs,
					ans, ans_dim);
		} else {
			/* method 7 */
			ret = read_data_7(h5dset,
					starts, breakpoint_bufs, chunkidx_bufs,
					ans, ans_dim);
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

/* If the H5 datatype that was used to store the logical data is an 8-bit
   or 16-bit signed integer type (e.g. H5T_STD_I8LE or H5T_STD_I16BE) then
   NA values got loaded as negative values that are not equal to NA_LOGICAL
   (e.g. as -128 for H5T_STD_I8LE and -2^16 for H5T_STD_I16BE).
   These values must be replaced with NA_LOGICAL. */
static void fix_logical_NAs(SEXP x)
{
	R_xlen_t x_len, i;
	int *x_p;

	x_len = XLENGTH(x);
	for (i = 0, x_p = LOGICAL(x); i < x_len; i++, x_p++) {
		if (*x_p < 0)
			*x_p = NA_LOGICAL;
	}
	return;
}

/* Return R_NilValue on error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce,
		    int as_int, int method)
{
	SEXP ans, ans_dim;
	H5DSetDescriptor h5dset;
	int ret, along;

	ans = R_NilValue;

	if (_init_H5DSetDescriptor(&h5dset, dset_id, as_int, 0) < 0)
		return ans;

	ret = _shallow_check_selection(h5dset.ndim, starts, counts);
	if (ret < 0)
		goto on_error;

	if (method < 0 || method > 7) {
		PRINT_TO_ERRMSG_BUF("'method' must be >= 0 and <= 7");
		goto on_error;
	}
	if (h5dset.Rtype == STRSXP) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "reading string data");
			goto on_error;
		}
		if (method == 0) {
			method = 4;
		} else if (method != 4 && method != 5) {
			PRINT_TO_ERRMSG_BUF("only methods 4 and 5 support "
					    "reading string data");
			goto on_error;
		}
	} else if (method == 0) {
		method = 1;
		/* March 27, 2019: My early testing (from Nov 2018) seemed
		   to indicate that method 6 was a better choice over method 4
		   when the layout is chunked and 'counts' is NULL. Turns out
		   that doing more testing today seems to indicate the opposite
		   i.e. method 4 now seems to perform better than method 6 on
		   all the datasets I've tested so far, including those used
		   by Pete Hickey here:
		     https://github.com/Bioconductor/DelayedArray/issues/13
		   and those used in the examples in man/h5mread.Rd.
		   Note sure what happened between Nov 2018 and today. Did I
		   do something stupid in my early testing? Did something
		   change in Rhdf5lib?
		   Anyway thanks to Pete for providing such a useful report.

		   Nov 26, 2019: I added method 7. Is like method 4 but
		   bypasses the intermediate buffer if a chunk is fully
		   selected. This is now preferred over methods 4 or 6. */
		if (h5dset.h5chunkdim != NULL &&
		    counts == R_NilValue &&
		    starts != R_NilValue)
		{
			for (along = 0; along < h5dset.ndim; along++) {
				if (VECTOR_ELT(starts, along) != R_NilValue) {
					//method = 6;
					//method = 4;
					method = 7;
					break;
				}
			}
		}
	}

	ans_dim = PROTECT(NEW_INTEGER(h5dset.ndim));

	if (method <= 3) {
		ans = h5mread_1_2_3(&h5dset, starts, counts, noreduce,
				    method, INTEGER(ans_dim));
	} else if (h5dset.h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("methods 4, 5, 6, and 7 cannot be used "
				    "on a contiguous dataset (unless\n"
				    "  it contains string data in which "
				    "case methods 4 and 5 can be used)");
	} else if (counts != R_NilValue) {
		PRINT_TO_ERRMSG_BUF("methods 4, 5, 6, and 7 can only be used "
				    "when 'counts' is NULL");
	} else {
		ans = h5mread_4_5_6_7(&h5dset, starts,
				      method, INTEGER(ans_dim));
	}
	if (ans != R_NilValue) {
		PROTECT(ans);
		if (TYPEOF(ans) == LGLSXP)
			fix_logical_NAs(ans);
		SET_DIM(ans, ans_dim);
		UNPROTECT(1);
	}

	UNPROTECT(1);

    on_error:
	_destroy_H5DSetDescriptor(&h5dset);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce,
	       SEXP as_integer, SEXP method)
{
	int noreduce0, as_int, method0;
	hid_t file_id, dset_id;
	SEXP ans;

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

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	ans = PROTECT(h5mread(dset_id, starts, counts, noreduce0,
			      as_int, method0));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	UNPROTECT(1);
	if (ans == R_NilValue)
		error(_HDF5Array_errmsg_buf());
	return ans;
}

