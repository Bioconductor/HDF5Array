/****************************************************************************
 *               Workhorses behind h5mread methods 4, 5, 6, 7               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_starts.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "h5mread_helpers.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */
//#include <time.h>


/****************************************************************************
 * read_data_4_5()
 *
 * One call to _read_H5Viewport() or _read_h5chunk() (wrappers for H5Dread()
 * or H5Dread_chunk(), respectively) per chunk touched by the user-supplied
 * array selection.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Make one call to _read_H5Viewport() or _read_h5chunk() to load the
 *     **entire** chunk data to an intermediate buffer.
 *   - Copy the user-selected data from the intermediate buffer to 'ans'.
 *
 * Assumes that 'h5dset->h5chunkdim' and 'h5dset->h5nchunk' are NOT
 * NULL. This is NOT checked!
 */

static void init_in_offset_and_out_offset(int ndim, SEXP starts,
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
		start = GET_LIST_ELT(starts, along);
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
			SEXP starts,
			const int *out_dim, const H5Viewport *dest_vp,
			const hsize_t *h5chunkdim,
			size_t *in_offset, size_t *out_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc, out_off_inc;

	start = GET_LIST_ELT(starts, inner_moved_along);
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
			start = GET_LIST_ELT(starts, along);
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

static int gather_selected_chunk_data(
		const H5DSetDescriptor *h5dset,
		SEXP starts, const void *in, const H5Viewport *tchunk_vp,
		SEXP ans, const int *out_dim,
		const H5Viewport *dest_vp, int *inner_midx_buf)
{
	int ndim, inner_moved_along, ret;
	size_t in_offset, out_offset;
	long long int num_elts;

	ndim = h5dset->ndim;
	init_in_offset_and_out_offset(ndim, starts,
			out_dim, dest_vp,
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
				starts,
				out_dim, dest_vp,
				h5dset->h5chunkdim,
				&in_offset, &out_offset);
	};
	//printf("# nb of selected elements in current chunk = %lld\n",
	//       num_elts);
	return 0;
}

static int read_data_from_chunk_4_5(const H5DSetDescriptor *h5dset, int method,
		SEXP starts,
		SEXP ans, const int *ans_dim,
		int *inner_midx_buf,
		const H5Viewport *tchunk_vp,
		const H5Viewport *middle_vp,
		const H5Viewport *dest_vp,
		void *chunk_data_buf, hid_t chunk_space_id,
		void *compressed_chunk_data_buf)
{
	int ret;

	if (method == 4) {
		/* It takes about 218s on my laptop to load all the chunks
		   from the EH1040 dataset (big 10x Genomics brain dataset
		   in dense format, chunks of 100x100, wrapped in the
		   TENxBrainData package). That's 60 microseconds per chunk! */
		ret = _read_H5Viewport(h5dset,
				tchunk_vp, middle_vp,
				chunk_data_buf, chunk_space_id);
	} else {
		ret = _read_h5chunk(h5dset,
				tchunk_vp,
				chunk_data_buf,
				compressed_chunk_data_buf);
	}
	if (ret < 0)
		return ret;
	ret = gather_selected_chunk_data(
			h5dset,
			starts, chunk_data_buf, tchunk_vp,
			ans, ans_dim,
			dest_vp, inner_midx_buf);
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
		const LLongAEAE *tchunkidx_bufs,
		const int *num_tchunks,
		SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ret;
	hid_t chunk_space_id;
	void *chunk_data_buf, *compressed_chunk_data_buf = NULL;
	H5Viewport tchunk_vp, middle_vp, dest_vp;
	IntAE *tchunk_midx_buf, *inner_midx_buf;
	long long int tchunk_rank;

	ndim = h5dset->ndim;

	/* Prepare buffers. */

	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}

	/* Allocate 'tchunk_vp', 'middle_vp', and 'dest_vp'.
	   We set 'dest_vp_mode' to ALLOC_OFF_AND_DIM because, in the
	   context of read_data_4_5(), we won't use 'dest_vp.h5off' or
	   'dest_vp.h5dim', only 'dest_vp.off' and 'dest_vp.dim'. */
	if (_alloc_tchunk_vp_middle_vp_dest_vp(ndim,
		&tchunk_vp, &middle_vp, &dest_vp,
		ALLOC_OFF_AND_DIM) < 0)
	{
		H5Sclose(chunk_space_id);
		return -1;
	}

	if (method == 4) {
		chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	} else {
		warning("method 5 is still experimental, use at your own risk");
		chunk_data_buf = malloc(2 * h5dset->chunk_data_buf_size +
					CHUNK_COMPRESSION_OVERHEAD);
	}
	if (chunk_data_buf == NULL) {
		_free_tchunk_vp_middle_vp_dest_vp(&tchunk_vp,
						   &middle_vp,
						   &dest_vp);
		H5Sclose(chunk_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}
	if (method != 4)
		compressed_chunk_data_buf = chunk_data_buf +
					    h5dset->chunk_data_buf_size;

	/* Walk over the chunks touched by the user-supplied array selection. */
	tchunk_rank = 0;
	moved_along = ndim;
	do {
		_update_tchunk_vp_dest_vp(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&tchunk_vp, &dest_vp);
		ret = read_data_from_chunk_4_5(h5dset, method,
			starts,
			ans, ans_dim,
			inner_midx_buf->elts,
			&tchunk_vp, &middle_vp, &dest_vp,
			chunk_data_buf, chunk_space_id,
			compressed_chunk_data_buf);
		if (ret < 0)
			break;
		tchunk_rank++;
		moved_along = _next_midx(ndim, num_tchunks,
					 tchunk_midx_buf->elts);
	} while (moved_along < ndim);
	free(chunk_data_buf);
	_free_tchunk_vp_middle_vp_dest_vp(&tchunk_vp, &middle_vp, &dest_vp);
	H5Sclose(chunk_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to _read_h5selection() (wrapper for H5Dread()) per chunk touched
 * by the user-supplied array selection. No intermediate buffer.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Select the hyperslabs obtained by intersecting the user-supplied
 *     array selection with the current chunk.
 *   - Call _read_h5selection(). This loads the selected data **directly**
 *     to the final destination array (ans).
 *
 * Assumes that 'h5dset->h5chunkdim' and 'h5dset->h5nchunk' are NOT
 * NULL. This is NOT checked!
 */

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP starts,
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
		start = GET_LIST_ELT(starts, along);
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

static void init_inner_vp(int ndim,
		SEXP starts,
		const H5Viewport *tchunk_vp,
		H5Viewport *inner_vp)
{
	int along, h5along;
	SEXP start;
	hsize_t d;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
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
		SEXP starts, const H5Viewport *dest_vp,
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
		if (starts == R_NilValue)
			continue;
		start = VECTOR_ELT(starts, along);
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

/* Return nb of selected elements (or -1 on error). */
static long long int NOT_USED_select_elements_from_chunk(
		const H5DSetDescriptor *h5dset,
		SEXP starts,
		const H5Viewport *dest_vp,
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
			i = dest_vp->off[along] + inner_midx_buf[along];
			start = GET_LIST_ELT(starts, along);
			if (start != R_NilValue) {
				coord = _get_trusted_elt(start, i) - 1;
			} else {
				coord = i;
			}
			*(coord_p++) = (hsize_t) coord;
		}
		inner_moved_along = _next_midx(ndim, dest_vp->dim,
					       inner_midx_buf);
	} while (inner_moved_along < ndim);

	ret = H5Sselect_elements(h5dset->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_intersection_of_chips_with_chunk(
		const H5DSetDescriptor *h5dset,
		SEXP starts,
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

	init_inner_vp(ndim, starts, tchunk_vp, inner_vp);

	/* Walk on the "inner chips" i.e. on the intersections between
	   the "chips" in the user-supplied array selection and the currrent
	   chunk. */
	num_hyperslabs = 0;
	inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_inner_vp(ndim, starts, dest_vp,
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

static int read_data_from_chunk_6(const H5DSetDescriptor *h5dset,
		const int *chunk_midx, int moved_along,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		void *dest, hid_t dest_space_id,
		int *inner_midx_buf,
		const H5Viewport *tchunk_vp,
		H5Viewport *inner_vp,
		const H5Viewport *dest_vp,
		IntAEAE *inner_breakpoint_bufs,
		const IntAE *inner_nchip_buf)
		// hsize_t *coord_buf)
{
	int ret;

	update_inner_breakpoints(h5dset->ndim, moved_along,
			starts, dest_vp,
			inner_breakpoint_bufs, inner_nchip_buf->elts);
	//t0 = clock();
	/* Having 'inner_nchip_buf->elts' identical to 'dest_vp->dim'
	   means that all the inner chips are single elements so we
	   use select_elements_from_chunk() which could be faster than
	   select_intersection_of_chips_with_chunk() in that case.
	   NO IT'S NOT FASTER! */
	//ret = memcmp(inner_nchip_buf->elts, dest_vp->dim,
	//	       ndim * sizeof(int));
	//printf("# chunk %lld: %s\n", tchunk_rank,
	//       ret == 0 ? "select_elements" : "select_hyperslabs");
	//if (ret == 0) {
	//	ret = select_elements_from_chunk(h5dset,
	//		starts, dest_vp,
	//		inner_midx_buf, coord_buf);
	//} else {
		ret = select_intersection_of_chips_with_chunk(
			h5dset, starts, dest_vp, tchunk_vp,
			inner_breakpoint_bufs, inner_nchip_buf->elts,
			inner_midx_buf, inner_vp);
	//}
	if (ret < 0)
		return ret;
	//t_select_elements += clock() - t0;

	//t0 = clock();
	ret = _read_h5selection(h5dset, dest_vp, dest, dest_space_id);
	//t_read_h5selection += clock() - t0;
	return ret;
}

static int read_data_6(const H5DSetDescriptor *h5dset,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		const int *num_tchunks,
		SEXP ans, const int *ans_dim)
{
	void *dest;
	int ndim, moved_along, ret;
	hid_t dest_space_id;
	H5Viewport tchunk_vp, inner_vp, dest_vp;
	//hsize_t *coord_buf;
	IntAE *tchunk_midx_buf, *inner_nchip_buf, *inner_midx_buf;
	IntAEAE *inner_breakpoint_bufs;
	long long int tchunk_rank;

	dest = DATAPTR(ans);
	if (dest == NULL)
		return -1;
	ndim = h5dset->ndim;
	dest_space_id = _create_mem_space(ndim, ans_dim);
	if (dest_space_id < 0)
		return -1;

	/* Allocate 'tchunk_vp', 'inner_vp', and 'dest_vp'. */
	if (_alloc_tchunk_vp_inner_vp_dest_vp(ndim,
		&tchunk_vp, &inner_vp, &dest_vp) < 0)
	{
		H5Sclose(dest_space_id);
		return -1;
	}

	/* Prepare buffers. */
	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_nchip_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	/* Walk over the chunks touched by the user-supplied array selection. */
	tchunk_rank = 0;
	moved_along = ndim;
	//clock_t t_select_elements = 0, t_read_h5selection = 0, t0;
	do {
		_update_tchunk_vp_dest_vp(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&tchunk_vp, &dest_vp);
		ret = read_data_from_chunk_6(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			dest, dest_space_id,
			inner_midx_buf->elts,
			&tchunk_vp, &inner_vp, &dest_vp,
			inner_breakpoint_bufs, inner_nchip_buf);
			//coord_buf);
		if (ret < 0)
			break;
		tchunk_rank++;
		moved_along = _next_midx(ndim, num_tchunks,
					 tchunk_midx_buf->elts);
	} while (moved_along < ndim);
	//free(coord_buf);
	_free_tchunk_vp_inner_vp_dest_vp(&tchunk_vp, &inner_vp, &dest_vp);
	H5Sclose(dest_space_id);
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

static int read_data_7(const H5DSetDescriptor *h5dset,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		const int *num_tchunks,
		SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ok, ret;
	hid_t chunk_space_id, dest_space_id;
	void *dest, *chunk_data_buf;
	H5Viewport tchunk_vp, middle_vp, dest_vp;
	IntAE *tchunk_midx_buf, *inner_midx_buf;
	long long int tchunk_rank;

	ndim = h5dset->ndim;
	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}
	dest = DATAPTR(ans);
	if (dest == NULL) {
		H5Sclose(chunk_space_id);
		return -1;
	}
	dest_space_id = _create_mem_space(ndim, ans_dim);
	if (dest_space_id < 0) {
		H5Sclose(chunk_space_id);
		return -1;
	}

	/* Allocate 'tchunk_vp', 'middle_vp', and 'dest_vp'.
	   Unlike in read_data_4_5(), here we set 'dest_vp_mode' to
	   ALLOC_ALL_FIELDS because, in the context of read_data_7(),
	   we will use 'dest_vp.h5off' and 'dest_vp.h5dim', not
	   just 'dest_vp.off' and 'dest_vp.dim'. */
	if (_alloc_tchunk_vp_middle_vp_dest_vp(ndim,
		&tchunk_vp, &middle_vp, &dest_vp,
		ALLOC_ALL_FIELDS) < 0)
	{
		H5Sclose(dest_space_id);
		H5Sclose(chunk_space_id);
		return -1;
	}

	chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	if (chunk_data_buf == NULL) {
		_free_tchunk_vp_middle_vp_dest_vp(&tchunk_vp,
						   &middle_vp,
						   &dest_vp);
		H5Sclose(dest_space_id);
		H5Sclose(chunk_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}

	/* Prepare buffers. */
	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	/* Walk over the chunks touched by the user-supplied array selection. */
	tchunk_rank = 0;
	moved_along = ndim;
	do {
		_update_tchunk_vp_dest_vp(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&tchunk_vp, &dest_vp);
		ok = _tchunk_is_fully_selected(h5dset->ndim,
					       &tchunk_vp, &dest_vp);
		if (ok) {
			/* Load the chunk **directly** into 'ans' (no
			   intermediate buffer). */
			ret = _read_H5Viewport(h5dset,
				&tchunk_vp, &dest_vp,
				dest, dest_space_id);
		} else {
			/* Load the **entire** chunk to an intermediate
			   buffer then copy the user-selected data from
			   the intermediate buffer to 'ans'. */
			ret = read_data_from_chunk_4_5(h5dset, 4,
				starts,
				ans, ans_dim,
				inner_midx_buf->elts,
				&tchunk_vp, &middle_vp, &dest_vp,
				chunk_data_buf, chunk_space_id, NULL);
		}
		if (ret < 0)
			break;
		tchunk_rank++;
		moved_along = _next_midx(ndim, num_tchunks,
					 tchunk_midx_buf->elts);
	} while (moved_along < ndim);
	free(chunk_data_buf);
	_free_tchunk_vp_middle_vp_dest_vp(&tchunk_vp, &middle_vp, &dest_vp);
	H5Sclose(dest_space_id);
	H5Sclose(chunk_space_id);
	return ret;
}


/****************************************************************************
 * _h5mread_starts()
 *
 * Implements methods 4 to 7.
 * Return an ordinary array or R_NilValue if an error occured.
 */

SEXP _h5mread_starts(const H5DSetDescriptor *h5dset, SEXP starts,
		     int method, int *ans_dim)
{
	int ndim, ret, along;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;  /* touched chunk ids along each dim */
	IntAE *ntchunk_buf;  /* nb of touched chunks along each dim */
	R_xlen_t ans_len;
	SEXP ans;

	ndim = h5dset->ndim;

	/* This call will populate 'ans_dim', 'breakpoint_bufs',
	   and 'tchunkidx_bufs'. */
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	tchunkidx_bufs = new_LLongAEAE(ndim, ndim);
	ret = _map_starts_to_h5chunks(h5dset, starts, ans_dim,
				      breakpoint_bufs, tchunkidx_bufs);
	if (ret < 0)
		return R_NilValue;

	ntchunk_buf = new_IntAE(ndim, ndim, 0);
	_set_num_tchunks(h5dset, starts, tchunkidx_bufs, ntchunk_buf->elts);

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(h5dset->Rtype, ans_len));

	/* ans_len != 0 means that the user-supplied array selection
	   is not empty */
	if (ans_len != 0) {
		if (method <= 5) {
			/* methods 4 and 5 */
			ret = read_data_4_5(h5dset, method, starts,
					breakpoint_bufs, tchunkidx_bufs,
					ntchunk_buf->elts,
					ans, ans_dim);
		} else if (method == 6) {
			/* method 6 */
			ret = read_data_6(h5dset, starts,
					breakpoint_bufs, tchunkidx_bufs,
					ntchunk_buf->elts,
					ans, ans_dim);
		} else {
			/* method 7 */
			ret = read_data_7(h5dset, starts,
					breakpoint_bufs, tchunkidx_bufs,
					ntchunk_buf->elts,
					ans, ans_dim);
		}
		if (ret < 0)
			goto on_error;
	}

	UNPROTECT(1);  /* 'ans' */
	return ans;

    on_error:
	UNPROTECT(1);  /* 'ans' */
	return R_NilValue;
}

