/****************************************************************************
 *                Workhorses behind h5mread methods 1, 2, 3                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_startscounts.h"

/*
 * Some useful links:
 * - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5S.html
 * - Documentation of H5Dread():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
 * - An H5Dread() example:
 *     https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
 */

#include "global_errmsg_buf.h"
#include "array_selection.h"
#include "h5mread_helpers.h"

#include <stdlib.h>  /* for malloc, free */
//#include <time.h>


/****************************************************************************
 * Low-level helpers
 */

static size_t set_nblocks(int ndim, SEXP starts,
			  const int *ans_dim, int expand, int *nblocks)
{
	size_t total_num_blocks;
	int along, nblock;
	SEXP start;

	total_num_blocks = 1;
	for (along = 0; along < ndim; along++) {
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			nblock = LENGTH(start);
		} else {
			nblock = expand ? ans_dim[along] : 1;
		}
		total_num_blocks *= nblocks[along] = nblock;
	}
	return total_num_blocks;
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

static void init_srcvp_buf(const H5DSetDescriptor *h5dset,
			   SEXP starts, H5Viewport *srcvp_buf)
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
			H5Viewport *srcvp_buf)
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
			int *nblocks, int *midx_buf)
{
	int ret, ndim, moved_along;
	H5Viewport srcvp_buf;
	long long int num_hyperslabs;

	ret = H5Sselect_none(h5dset->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = h5dset->ndim;
	set_nblocks(ndim, starts, ans_dim, 0, nblocks);

	/* Allocate 'srcvp_buf'. */
	if (_alloc_H5Viewport(&srcvp_buf, ndim, 0) < 0)
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
		ret = _add_H5Viewport_to_selection(h5dset->space_id,
						   &srcvp_buf);
		if (ret < 0)
			break;
		moved_along = _next_midx(ndim, nblocks, midx_buf);
	} while (moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);

	_free_H5Viewport(&srcvp_buf);
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
			int *nblocks, int *midx_buf)
{
	int ndim, outer_moved_along, ret;
	size_t num_elements;
	hsize_t *coord_buf, *coord_p;

	ndim = h5dset->ndim;
	num_elements = set_nblocks(ndim, starts, ans_dim, 1, nblocks);

	/* Allocate 'coord_buf'. */
	coord_buf = _alloc_hsize_t_buf(num_elements * ndim, 0, "'coord_buf'");
	if (coord_buf == NULL)
		return -1;

	/* Walk on the selected elements. */
	coord_p = coord_buf;
	outer_moved_along = ndim;
	do {
		coord_p = add_element(ndim, midx_buf, starts, coord_p);
		outer_moved_along = _next_midx(ndim, nblocks, midx_buf);
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
	long long int total_num_blocks;

	ndim = h5dset->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	//clock_t t0 = clock();
	if (method == 1) {
		total_num_blocks = select_hyperslabs(h5dset,
					starts, counts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	} else {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "'method' is set to 2");
			return -1;
		}
		total_num_blocks = select_elements(h5dset,
					starts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	}
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("time for setting selection: %e\n", dt);
	//printf("total_num_blocks: %lld, time per block: %e\n",
	//	total_num_blocks, dt / total_num_blocks);
	return total_num_blocks < 0 ? -1 : 0;
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
		H5Viewport *srcvp_buf, H5Viewport *destvp_buf,
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

	ret = _select_H5Viewport(h5dset->space_id, srcvp_buf);
	if (ret < 0)
		return -1;
	ret = _select_H5Viewport(mem_space_id, destvp_buf);
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
	H5Viewport srcvp_buf, destvp_buf;
	IntAE *nblock_buf, *midx_buf;
	long long int num_hyperslabs;

	ndim = h5dset->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	set_nblocks(ndim, starts, ans_dim, 0, nblock_buf->elts);

	/* Allocate 'srcvp_buf' and 'destvp_buf'. */
	if (_alloc_H5Viewport(&srcvp_buf, ndim, 0) < 0)
		return -1;
	destvp_buf.h5off = _alloc_hsize_t_buf(ndim, 1, "'destvp_buf.h5off'");
	destvp_buf.h5dim = srcvp_buf.h5dim;
	if (destvp_buf.h5off == NULL) {
		_free_H5Viewport(&srcvp_buf);
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
		moved_along = _next_midx(ndim, nblock_buf->elts,
					 midx_buf->elts);
	} while (moved_along < ndim);

	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);
	free(destvp_buf.h5off);
	_free_H5Viewport(&srcvp_buf);
	return ret;
}


/****************************************************************************
 * _h5mread_startscounts()
 */

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

/* Implements methods 1 to 3.
   Return an ordinary array or R_NilValue if an error occured. */
SEXP _h5mread_startscounts(const H5DSetDescriptor *h5dset,
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
		mem = DATAPTR(ans);
		if (mem == NULL)
			goto on_error;
		mem_space_id = _create_mem_space(ndim, ans_dim);
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

