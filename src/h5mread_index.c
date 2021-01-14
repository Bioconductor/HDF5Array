/****************************************************************************
 *                 Workhorses behind h5mread methods 4, 5, 6                *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_index.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "h5mread_helpers.h"
#include "ChunkIterator.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */
//#include <time.h>


/****************************************************************************
 * read_data_4_5()
 *
 * method 4: One call to _read_H5Viewport() or _read_h5chunk() (wrappers for
 * H5Dread() or H5Dread_chunk(), respectively) per chunk touched by the
 * user-supplied array selection.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Make one call to _read_H5Viewport() or _read_h5chunk() to load the
 *     **entire** chunk data to an intermediate buffer.
 *   - Copy the user-selected data from the intermediate buffer to 'ans'.
 *
 * method 5: Like read_data_4_5() but bypasses the intermediate buffer if a
 * chunk is fully selected.
 */

static void init_in_offset_and_out_offset(int ndim, SEXP index,
			const int *out_dim, const H5Viewport *dest_vp,
			const H5Viewport *tchunk_vp,
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
		i = dest_vp->off[along];
		start = GET_LIST_ELT(index, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  tchunk_vp->h5off[h5along];
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
			SEXP index,
			const int *out_dim, const H5Viewport *dest_vp,
			const hsize_t *h5chunkdim,
			size_t *in_offset, size_t *out_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc, out_off_inc;

	start = GET_LIST_ELT(index, inner_moved_along);
	if (start != R_NilValue) {
		i1 = dest_vp->off[inner_moved_along] +
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
			di = 1 - dest_vp->dim[along];
			start = GET_LIST_ELT(index, along);
			if (start != R_NilValue) {
				i1 = dest_vp->off[along];
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

static int load_val_to_array(const H5DSetDescriptor *h5dset,
		const void *in, size_t in_offset,
		void *out, size_t out_offset)
{
	switch (h5dset->Rtype) {
	    case LGLSXP: {
		int val = ((int *) in)[in_offset];
		LOGICAL(out)[out_offset] = val;
	    } break;
	    case INTSXP: {
		int val = ((int *) in)[in_offset];
		INTEGER(out)[out_offset] = val;
	    } break;
	    case REALSXP: {
		double val = ((double *) in)[in_offset];
		REAL(out)[out_offset] = val;
	    } break;
	    case STRSXP: {
		const char *val = ((char *) in) + in_offset * h5dset->H5size;
		int val_len;
		for (val_len = 0; val_len < h5dset->H5size; val_len++)
			if (val[val_len] == 0)
				break;
		int is_na = h5dset->as_na_attr &&
			    val_len == 2 && val[0] == 'N' && val[1] == 'A';
		if (is_na) {
			SET_STRING_ELT(out, out_offset, NA_STRING);
		} else {
			SEXP out_elt = PROTECT(mkCharLen(val, val_len));
			SET_STRING_ELT(out, out_offset, out_elt);
			UNPROTECT(1);
		}
	    } break;
	    case RAWSXP: {
		char val = ((char *) in)[in_offset];
		RAW(out)[out_offset] = val;
	    } break;
	    default:
		PRINT_TO_ERRMSG_BUF("unsupported type: %s",
				    CHAR(type2str(h5dset->Rtype)));
		return -1;
	}
	return 0;
}

static int copy_selected_chunk_data_to_ans(
		const H5DSetDescriptor *h5dset, SEXP index,
		const H5Viewport *tchunk_vp,
		const void *in,
		const H5Viewport *dest_vp, int *inner_midx_buf,
		const int *ans_dim, SEXP ans)
{
	int ndim, inner_moved_along, ret;
	size_t in_offset, out_offset;
	long long int num_elts;

	ndim = h5dset->ndim;
	init_in_offset_and_out_offset(ndim, index,
			ans_dim, dest_vp,
			tchunk_vp, h5dset->h5chunkdim,
			&in_offset, &out_offset);
	/* Walk on the selected elements in current chunk. */
	num_elts = 0;
	while (1) {
		ret = load_val_to_array(h5dset, in, in_offset, ans, out_offset);
		if (ret < 0)
			return ret;
		num_elts++;
		inner_moved_along = _next_midx(ndim, dest_vp->dim,
					       inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset_and_out_offset(ndim,
				inner_midx_buf, inner_moved_along,
				index,
				ans_dim, dest_vp,
				h5dset->h5chunkdim,
				&in_offset, &out_offset);
	};
	//printf("# nb of selected elements in current chunk = %lld\n",
	//       num_elts);
	return 0;
}

/*
  WARNING: 'use.H5Dread_chunk=TRUE' is not working properly on some datasets:
      library(HDF5Array)
      library(ExperimentHub)
      hub <- ExperimentHub()
      fname0 <- hub[["EH1039"]]
      h5mread(fname0, "mm10/barcodes", list(1), method=4L)
      # [1] "AAACCTGAGATAGGAG-1"
      h5mread(fname0, "mm10/barcodes", list(1),
              method=4L, use.H5Dread_chunk=TRUE)
      # [1] "AAAAAAAAAAAAAAAAAAAA"
  Looks like the chunk data has been shuffled (transposed in that case)
  before being written to disk in order to improve compression.
  TODO: Investigate this further. I suspect we need to check whether a
  "Data shuffling filter" (H5Z_FILTER_SHUFFLE) was used at creation time.
  Check H5Pget_filter() for how to know whether this filter was used or not.
  There should be a way to retrieve information about how the data was
  shuffled.
 */
static int read_data_4_5(ChunkIterator *chunk_iter,
		const int *ans_dim, SEXP ans,
		int method, int use_H5Dread_chunk)
{
	const H5DSetDescriptor *h5dset;
	int ndim, ret, direct_load;
	IntAE *inner_midx_buf;
	void *dest;
	hid_t dest_space_id;
	ChunkDataBuffer chunk_data_buf;

	if (use_H5Dread_chunk)
		warning("using 'use.H5Dread_chunk=TRUE' is still "
			"experimental, use at your own risk");

	h5dset = chunk_iter->h5dset;
	ndim = h5dset->ndim;
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	dest = DATAPTR(ans);
	if (dest == NULL)
		return -1;

	dest_space_id = _create_mem_space(ndim, ans_dim);
	if (dest_space_id < 0)
		return -1;

	ret = _init_ChunkDataBuffer(&chunk_data_buf, h5dset);
	if (ret < 0) {
		H5Sclose(dest_space_id);
		return ret;
	}
	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		//_print_tchunk_info(chunk_iter);
		direct_load = method == 5 &&
			      _tchunk_is_fully_selected(ndim,
					&chunk_iter->tchunk_vp,
					&chunk_iter->dest_vp);
		if (direct_load) {
			/* Load the chunk **directly** into 'ans' (no
			   intermediate buffer). */
			ret = _read_H5Viewport(h5dset,
				&chunk_iter->tchunk_vp,
				&chunk_iter->dest_vp,
				dest, dest_space_id);
		} else {
			/* Load the **entire** chunk to an intermediate
			   buffer then copy the user-selected chunk data
			   from the intermediate buffer to 'ans'. */
			ret = _load_chunk(chunk_iter, &chunk_data_buf,
					  use_H5Dread_chunk);
			if (ret < 0)
				break;
			ret = copy_selected_chunk_data_to_ans(
					chunk_iter->h5dset,
					chunk_iter->index,
					&chunk_iter->tchunk_vp,
					chunk_data_buf.data,
					&chunk_iter->dest_vp,
					inner_midx_buf->elts,
					ans_dim, ans);
		}
		if (ret < 0)
			break;
	}
	_destroy_ChunkDataBuffer(&chunk_data_buf);
	H5Sclose(dest_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to _read_h5selection() (wrapper for H5Dread()) per chunk touched
 * by the user-supplied array selection. No intermediate buffer.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Select the hyperslabs obtained by intersecting the user-supplied
 *     array selection with the current chunk.
 *   - Call _read_h5selection(). This loads the selected data **directly**
 *     to the final destination array (ans).
 */

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP index,
		const H5Viewport *dest_vp,
		IntAEAE *inner_breakpoint_bufs, int *inner_nchip_buf)
{
	int along, d, off, nchip, i;
	IntAE *inner_breakpoint_buf;
	SEXP start;
	long long int s0, s1;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		inner_breakpoint_buf = inner_breakpoint_bufs->elts[along];
		IntAE_set_nelt(inner_breakpoint_buf, 0);
		d = dest_vp->dim[along];
		start = GET_LIST_ELT(index, along);
		if (start == R_NilValue) {
			IntAE_insert_at(inner_breakpoint_buf, 0, d);
			inner_nchip_buf[along] = 1;
			continue;
		}
		off = dest_vp->off[along];
		s1 = _get_trusted_elt(start, off);
		nchip = 0;
		for (i = 1; i < d; i++) {
			s0 = s1;
			s1 = _get_trusted_elt(start, off + i);
			if (s1 != s0 + 1)
				IntAE_insert_at(inner_breakpoint_buf,
						nchip++, i);
		}
		IntAE_insert_at(inner_breakpoint_buf, nchip++, d);
		inner_nchip_buf[along] = nchip;
	}
	return;
}

static void init_inner_vp(int ndim, SEXP index,
		const H5Viewport *tchunk_vp,
		H5Viewport *inner_vp)
{
	int along, h5along;
	SEXP start;
	hsize_t d;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(index, along);
		if (start == R_NilValue) {
			inner_vp->h5off[h5along] =
					tchunk_vp->h5off[h5along];
			d = tchunk_vp->h5dim[h5along];
		} else {
			d = 1;
		}
		inner_vp->h5dim[h5along] = d;
	}
	return;
}

static void update_inner_vp(int ndim,
		SEXP index, const H5Viewport *dest_vp,
		const int *inner_midx, int inner_moved_along,
		const IntAEAE *inner_breakpoint_bufs,
		H5Viewport *inner_vp)
{
	int along, h5along, idx, off, d, i;
	SEXP start;
	const int *inner_breakpoint;

	for (along = 0; along < ndim; along++) {
		if (along > inner_moved_along)
			break;
		if (index == R_NilValue)
			continue;
		start = VECTOR_ELT(index, along);
		if (start == R_NilValue)
			continue;
		inner_breakpoint = inner_breakpoint_bufs->elts[along]->elts;
		idx = inner_midx[along];
		off = idx == 0 ? 0 : inner_breakpoint[idx - 1];
		d = inner_breakpoint[idx] - off;
		i = dest_vp->off[along] + off;
		h5along = ndim - 1 - along;
		inner_vp->h5off[h5along] = _get_trusted_elt(start, i) - 1;
		inner_vp->h5dim[h5along] = d;
	}
	return;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_intersection_of_chips_with_chunk(
		const H5DSetDescriptor *h5dset, SEXP index,
		const H5Viewport *dest_vp, const H5Viewport *tchunk_vp,
		const IntAEAE *inner_breakpoint_bufs,
		const int *inner_nchip,
		int *inner_midx_buf,
		H5Viewport *inner_vp)
{
	int ret, ndim, inner_moved_along;
	long long int num_hyperslabs;

	ret = H5Sselect_none(h5dset->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = h5dset->ndim;

	init_inner_vp(ndim, index, tchunk_vp, inner_vp);

	/* Walk on the "inner chips" i.e. on the intersections between
	   the "chips" in the user-supplied array selection and the currrent
	   chunk. */
	num_hyperslabs = 0;
	inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_inner_vp(ndim, index, dest_vp,
				inner_midx_buf, inner_moved_along,
				inner_breakpoint_bufs,
				inner_vp);
		ret = _add_H5Viewport_to_h5selection(h5dset->space_id,
						     inner_vp);
		if (ret < 0)
			return -1;
		inner_moved_along = _next_midx(ndim, inner_nchip,
					       inner_midx_buf);
	} while (inner_moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(inner_nchip)

	return num_hyperslabs;
}

static int read_data_from_chunk_6(const ChunkIterator *chunk_iter,
		void *dest, hid_t dest_space_id,
		int *inner_midx_buf,
		H5Viewport *inner_vp,
		IntAEAE *inner_breakpoint_bufs,
		const IntAE *inner_nchip_buf)
{
	const H5DSetDescriptor *h5dset;
	int ret;

	h5dset = chunk_iter->h5dset;
	update_inner_breakpoints(h5dset->ndim, chunk_iter->moved_along,
			chunk_iter->index, &chunk_iter->dest_vp,
			inner_breakpoint_bufs, inner_nchip_buf->elts);
	ret = select_intersection_of_chips_with_chunk(
			h5dset, chunk_iter->index,
			&chunk_iter->dest_vp, &chunk_iter->tchunk_vp,
			inner_breakpoint_bufs, inner_nchip_buf->elts,
			inner_midx_buf, inner_vp);
	if (ret < 0)
		return ret;
	ret = _read_h5selection(h5dset, &chunk_iter->dest_vp,
				dest, dest_space_id);
	return ret;
}

static int read_data_6(ChunkIterator *chunk_iter, const int *ans_dim, SEXP ans)
{
	const H5DSetDescriptor *h5dset;
	int ndim, ret;
	IntAE *inner_midx_buf, *inner_nchip_buf;
	IntAEAE *inner_breakpoint_bufs;
	void *dest;
	hid_t dest_space_id;
	H5Viewport inner_vp;
	ChunkDataBuffer chunk_data_buf;

	h5dset = chunk_iter->h5dset;
	ndim = h5dset->ndim;
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_nchip_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	dest = DATAPTR(ans);
	if (dest == NULL)
		return -1;
	dest_space_id = _create_mem_space(ndim, ans_dim);
	if (dest_space_id < 0)
		return -1;

	ret = _alloc_H5Viewport(&inner_vp, ndim, ALLOC_H5OFF_AND_H5DIM);
	if (ret < 0) {
		H5Sclose(dest_space_id);
		return ret;
	}

	ret = _init_ChunkDataBuffer(&chunk_data_buf, h5dset);
	if (ret < 0) {
		_free_H5Viewport(&inner_vp);
		H5Sclose(dest_space_id);
		return ret;
	}
	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		ret = read_data_from_chunk_6(chunk_iter,
			dest, dest_space_id,
			inner_midx_buf->elts,
			&inner_vp,
			inner_breakpoint_bufs, inner_nchip_buf);
		if (ret < 0)
			break;
	}
	_destroy_ChunkDataBuffer(&chunk_data_buf);
	_free_H5Viewport(&inner_vp);
	H5Sclose(dest_space_id);
	return ret;
}


/****************************************************************************
 * _h5mread_index()
 *
 * Implements methods 4 to 6.
 * Return an ordinary array or R_NilValue if an error occured.
 */

SEXP _h5mread_index(const H5DSetDescriptor *h5dset, SEXP index,
		    int method, int use_H5Dread_chunk, int *ans_dim)
{
	ChunkIterator chunk_iter;
	int ret, ndim, along;
	R_xlen_t ans_len;
	SEXP ans;

	/* In the context of methods 5 & 6, 'chunk_iter.dest_vp.h5off'
	   and 'chunk_iter.dest_vp.h5dim' will be used, not just
	   'chunk_iter.dest_vp.off' and 'chunk_iter.dest_vp.dim',
	   so we set 'alloc_full_dest_vp' (last arg) to 1. */
	ret = _init_ChunkIterator(&chunk_iter, h5dset, index,
				  ans_dim, method == 5 || method == 6);
	if (ret < 0)
		return R_NilValue;

	ndim = h5dset->ndim;
	for (along = 0, ans_len = 1; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(h5dset->Rtype, ans_len));

	if (method <= 5) {
		/* methods 4 and 5 */
		ret = read_data_4_5(&chunk_iter, ans_dim, ans,
				    method, use_H5Dread_chunk);
	} else {
		/* method 6 */
		ret = read_data_6(&chunk_iter, ans_dim, ans);
	}
	_destroy_ChunkIterator(&chunk_iter);
	UNPROTECT(1);
	return ret < 0 ? R_NilValue : ans;
}

