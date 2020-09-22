/****************************************************************************
 *                    Workhorse behind h5mread method 8                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_sparse.h"

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
#include <string.h>  /* for strncpy, memcmp */
#include <limits.h>  /* for INT_MAX */
//#include <time.h>


/****************************************************************************
 * Help manage memory of 'chunkvp', 'middlevp', and 'destvp' buffers
 */

static int alloc_chunkvp_middlevp_destvp_bufs(int ndim,
			H5Viewport *chunkvp_buf,
			H5Viewport *middlevp_buf,
			H5Viewport *destvp_buf, int destvp_mode)
{
	if (_alloc_H5Viewport(chunkvp_buf, ndim, 0) < 0)
		return -1;
	middlevp_buf->h5off =
		_alloc_hsize_t_buf(ndim, 1, "'middlevp_buf->h5off'");
	middlevp_buf->h5dim = chunkvp_buf->h5dim;
	if (middlevp_buf->h5off == NULL) {
		_free_H5Viewport(chunkvp_buf);
		return -1;
	}
	if (_alloc_H5Viewport(destvp_buf, ndim, destvp_mode) < 0) {
		free(middlevp_buf->h5off);
		_free_H5Viewport(chunkvp_buf);
		return -1;
	}
	return 0;
}

static void free_chunkvp_middlevp_destvp_bufs(
			H5Viewport *chunkvp_buf,
			H5Viewport *middlevp_buf,
			H5Viewport *destvp_buf)
{
	_free_H5Viewport(destvp_buf);
	free(middlevp_buf->h5off);
	_free_H5Viewport(chunkvp_buf);
	return;
}


/****************************************************************************
 * Manipulation of 'nzdata' and 'nzindex' buffers
 */

static void *new_nzdata_buf(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return new_IntAE(0, 0, 0);
	    case REALSXP:             return new_DoubleAE(0, 0, 0.0);
	    case STRSXP:              return new_CharAEAE(0, 0);
	    case RAWSXP:              return new_CharAE(0);
	}
	/* Should never happen. The early call to _init_H5DSetDescriptor()
	   in h5mread() already made sure that Rtype is supported. */
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return NULL;
}

static SEXP make_nzdata_from_buf(const void *nzdata_buf, SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP:  return new_LOGICAL_from_IntAE(nzdata_buf);
	    case INTSXP:  return new_INTEGER_from_IntAE(nzdata_buf);
	    case REALSXP: return new_NUMERIC_from_DoubleAE(nzdata_buf);
	    case STRSXP:  return new_CHARACTER_from_CharAEAE(nzdata_buf);
	    case RAWSXP:  return new_RAW_from_CharAE(nzdata_buf);
	}
	/* Should never happen. The early call to _init_H5DSetDescriptor()
	   in h5mread() already made sure that Rtype is supported. */
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return R_NilValue;
}

static SEXP make_nzindex_from_buf(const IntAE *nzindex_buf,
				  R_xlen_t nzdata_len, int ndim)
{
	SEXP nzindex;
	size_t nzindex_len;
	const int *in_p;
	int i, along;

	nzindex_len = IntAE_get_nelt(nzindex_buf);
	if (nzindex_len != (size_t) nzdata_len * ndim) {
		/* Should never happen. */
		PRINT_TO_ERRMSG_BUF("HDF5Array internal error in "
				    "C function make_nzindex_from_buf(): "
				    "length(nzindex) != length(nzdata) * ndim");
		return R_NilValue;
	}
	/* 'nzdata_len' is guaranteed to be <= INT_MAX otherwise the
	   earlier calls to load_nonzero_val_to_nzdata_buf() (see above)
	   would have raised an error. */
	nzindex = PROTECT(allocMatrix(INTSXP, (int) nzdata_len, ndim));
	in_p = nzindex_buf->elts;
	for (i = 0; i < nzdata_len; i++) {
		for (along = 0; along < ndim; along++) {
			INTEGER(nzindex)[i + nzdata_len * along] = *(in_p++);
		}
	}
	UNPROTECT(1);
	return nzindex;
}

/* 'make_nzdata_from_IntAE_bufs()' and 'make_nzindex_from_bufs()' are not used
   at the moment. Would be used if we were using one nzdata and one nzindex
   buffer per touched chunk (e.g. an IntAEAE of length the total nb of touched
   chunks) instead of a single nzdata and nzindex buffer for all the chunks
   like we do at the moment. */
static SEXP make_nzdata_from_IntAE_bufs(const IntAEAE *nzdata_bufs,
					SEXPTYPE Rtype)
{
	size_t total_num_tchunks, nzdata_len, tchunk_rank, buf_nelt;
	const IntAE *buf;
	SEXP nzdata;
	int *out_p;

	total_num_tchunks = IntAEAE_get_nelt(nzdata_bufs);
	nzdata_len = 0;
	for (tchunk_rank = 0; tchunk_rank < total_num_tchunks; tchunk_rank++) {
		buf = nzdata_bufs->elts[tchunk_rank];
		buf_nelt = IntAE_get_nelt(buf);
		nzdata_len += buf_nelt;
	}
	nzdata = PROTECT(allocVector(Rtype, nzdata_len));
	out_p = (int *) DATAPTR(nzdata);
	for (tchunk_rank = 0; tchunk_rank < total_num_tchunks; tchunk_rank++) {
		buf = nzdata_bufs->elts[tchunk_rank];
		buf_nelt = IntAE_get_nelt(buf);
		memcpy(out_p, buf->elts, sizeof(int) * buf_nelt);
		out_p += buf_nelt;
	}
	UNPROTECT(1);
	return nzdata;
}

static SEXP make_nzindex_from_bufs(const IntAEAE *nzindex_bufs,
				   R_xlen_t nzdata_len, int ndim)
{
	size_t total_num_tchunks, tchunk_rank, buf_nelt;
	const IntAE *buf;
	SEXP nzindex;
	int *out_p, i, buf_nrow, along;
	const int *in_p;

	/* 'nzdata_len' is guaranteed to be <= INT_MAX otherwise the
	   earlier calls to load_nonzero_val_to_nzdata_buf() (see above)
	   would have raised an error. */
	nzindex = PROTECT(allocMatrix(INTSXP, (int) nzdata_len, ndim));
	out_p = INTEGER(nzindex);
	total_num_tchunks = IntAEAE_get_nelt(nzindex_bufs);
	for (tchunk_rank = 0; tchunk_rank < total_num_tchunks; tchunk_rank++) {
		buf = nzindex_bufs->elts[tchunk_rank];
		buf_nelt = IntAE_get_nelt(buf);
		buf_nrow = buf_nelt / ndim;
		in_p = buf->elts;
		for (i = 0; i < buf_nrow; i++) {
			for (along = 0; along < ndim; along++)
				out_p[nzdata_len * along] = *(in_p++);
			out_p++;
		}
	}
	UNPROTECT(1);
	return nzindex;
}


/****************************************************************************
 * read_data_8()
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

static void update_chunkvp_buf(const H5DSetDescriptor *h5dset,
			const int *chunk_midx, int moved_along,
			SEXP starts, const LLongAEAE *tchunkidx_bufs,
			H5Viewport *chunkvp_buf)
{
	int ndim, along, h5along, i;
	SEXP start;
	long long int tchunkidx;
	hsize_t chunkd, off, d;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = chunk_midx[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			tchunkidx = tchunkidx_bufs->elts[along]->elts[i];
		} else {
			tchunkidx = i;
		}
		chunkd = h5dset->h5chunkdim[h5along];
		off = tchunkidx * chunkd;
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
			const H5Viewport *chunkvp, H5Viewport *destvp_buf)
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
			const LLongAEAE *tchunkidx_bufs,
			H5Viewport *chunkvp_buf, H5Viewport *destvp_buf)
{
	update_chunkvp_buf(h5dset,
			chunk_midx, moved_along,
			starts, tchunkidx_bufs,
			chunkvp_buf);
	update_destvp_buf(h5dset,
			chunk_midx, moved_along,
			starts, breakpoint_bufs,
			chunkvp_buf, destvp_buf);
	return;
}

static void init_in_offset(int ndim, SEXP starts,
			   const H5Viewport *destvp,
			   const H5Viewport *chunkvp,
			   const hsize_t *h5chunkdim,
			   size_t *in_offset)
{
	size_t in_off;
	int along, h5along, i;
	SEXP start;

	in_off = 0;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunkdim[h5along];
		i = destvp->off[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  chunkvp->h5off[h5along];
	}
	*in_offset = in_off;
	return;
}

static inline void update_in_offset(int ndim,
			const int *inner_midx, int inner_moved_along,
			SEXP starts,
			const H5Viewport *destvp,
			const hsize_t *h5chunkdim,
			size_t *in_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc;

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
	if (inner_moved_along >= 1) {
		along = inner_moved_along - 1;
		h5along = ndim - inner_moved_along;
		do {
			in_off_inc *= h5chunkdim[h5along];
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
			along--;
			h5along++;
		} while (along >= 0);
	}
	*in_offset += in_off_inc;
	return;
}

static inline int load_nonzero_val_to_nzdata_buf(
		const H5DSetDescriptor *h5dset,
		const int *in, size_t in_offset,
		void *nzdata_buf)
{
	size_t nzdata_len;

	/* The maximum 'nzdata' length that we support is INT_MAX. This
	   is because at the end of the read_data_8() call we will need
	   to turn 'nzindex_buf' into an ordinary matrix with 'length(nzdata)'
	   rows. However R does not support ordinary matrices or arrays
	   with dimensions >= INT_MAX yet.
	   The function that will turn 'nzindex_buf' into an ordinary matrix
	   is make_nzindex_from_buf(). See below.
	*/
	switch (h5dset->Rtype) {
	    case LGLSXP:
	    case INTSXP: {
		int val = ((int *) in)[in_offset];
		if (val == 0)
			return 0;
		IntAE *buf = (IntAE *) nzdata_buf;
		nzdata_len = IntAE_get_nelt(buf);
		if (nzdata_len >= INT_MAX) {
			PRINT_TO_ERRMSG_BUF("too many non-zero values to load");
			return -1;
		}
		IntAE_insert_at(buf, nzdata_len, val);
	    } break;
	    case REALSXP: {
		double val = ((double *) in)[in_offset];
		if (val == 0.0)
			return 0;
		DoubleAE *buf = (DoubleAE *) nzdata_buf;
		nzdata_len = DoubleAE_get_nelt(buf);
		if (nzdata_len >= INT_MAX) {
			PRINT_TO_ERRMSG_BUF("too many non-zero values to load");
			return -1;
		}
		DoubleAE_insert_at(buf, nzdata_len, val);
	    } break;
	    case STRSXP: {
		const char *val = ((char *) in) + in_offset * h5dset->H5size;
		int val_len;
		for (val_len = 0; val_len < h5dset->H5size; val_len++)
			if (val[val_len] == 0)
				break;
		if (val_len == 0)
			return 0;
		CharAEAE *buf = (CharAEAE *) nzdata_buf;
		nzdata_len = CharAEAE_get_nelt(buf);
		if (nzdata_len >= INT_MAX) {
			PRINT_TO_ERRMSG_BUF("too many non-zero values to load");
			return -1;
		}
		CharAE *buf_elt = new_CharAE(val_len);
		strncpy(buf_elt->elts, val, val_len);
		CharAEAE_insert_at(buf, nzdata_len, buf_elt);
	    } break;
	    case RAWSXP: {
		char val = ((char *) in)[in_offset];
		if (val == 0)
			return 0;
		CharAE *buf = (CharAE *) nzdata_buf;
		nzdata_len = CharAE_get_nelt(buf);
		if (nzdata_len >= INT_MAX) {
			PRINT_TO_ERRMSG_BUF("too many non-zero values to load");
			return -1;
		}
		CharAE_insert_at(buf, nzdata_len, val);
	    } break;
	    default:
		PRINT_TO_ERRMSG_BUF("unsupported type: %s",
				    CHAR(type2str(h5dset->Rtype)));
		return -1;
	}
	return 1;
}

static inline void load_midx_to_nzindex_buf(int ndim, const H5Viewport *destvp,
		const int *inner_midx_buf,
		IntAE *nzindex_buf)
{
	int along, midx;

	for (along = 0; along < ndim; along++) {
		midx = destvp->off[along] + inner_midx_buf[along] + 1;
		IntAE_insert_at(nzindex_buf, IntAE_get_nelt(nzindex_buf), midx);
	}
	return;
}

static int gather_selected_chunk_data_as_sparse(
		const H5DSetDescriptor *h5dset,
		SEXP starts, const void *in, const H5Viewport *chunkvp,
		IntAE *nzindex_buf, void *nzdata_buf,
		const H5Viewport *destvp, int *inner_midx_buf)
{
	int ndim, inner_moved_along, ret;
	size_t in_offset;

	ndim = h5dset->ndim;

	/* Walk on the selected elements in current chunk and count the
	   non-zero ones. */
	init_in_offset(ndim, starts, destvp, chunkvp,
		       h5dset->h5chunkdim, &in_offset);
	while (1) {
		ret = load_nonzero_val_to_nzdata_buf(h5dset,
					in, in_offset, nzdata_buf);
		if (ret < 0)
			return ret;
		if (ret > 0)
			load_midx_to_nzindex_buf(ndim, destvp, inner_midx_buf,
						 nzindex_buf);
		inner_moved_along = _next_midx(ndim, destvp->dim,
					       inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset(ndim,
				 inner_midx_buf, inner_moved_along,
				 starts,
				 destvp,
				 h5dset->h5chunkdim,
				 &in_offset);
	};
	return 0;
}

static int read_data_from_chunk_8(const H5DSetDescriptor *h5dset,
			SEXP starts,
			IntAE *nzindex_buf, void *nzdata_buf,
			int *inner_midx_buf,
			H5Viewport *chunkvp_buf,
			const H5Viewport *middlevp,
			H5Viewport *destvp_buf,
			void *chunk_data_buf, hid_t chunk_space_id)
{
	int ret;

	ret = _read_H5Viewport(h5dset, chunkvp_buf, middlevp,
			       chunk_data_buf, chunk_space_id);
	if (ret < 0)
		return ret;
	return gather_selected_chunk_data_as_sparse(
			h5dset,
			starts, chunk_data_buf, chunkvp_buf,
			nzindex_buf, nzdata_buf,
			destvp_buf, inner_midx_buf);
}

static int read_data_8(const H5DSetDescriptor *h5dset,
			 SEXP starts,
			 const IntAEAE *breakpoint_bufs,
			 const LLongAEAE *tchunkidx_bufs,
			 const int *ntchunks,
			 SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ret;
	hid_t chunk_space_id;
	void *chunk_data_buf;
	H5Viewport chunkvp_buf, middlevp_buf, destvp_buf;
	IntAE *tchunk_midx_buf, *inner_midx_buf, *nzindex_buf;
	void *nzdata_buf;
	long long int tchunk_rank;

	ndim = h5dset->ndim;

	/* Prepare buffers. */

	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	nzindex_buf = new_IntAE(ndim, 0, 0);
	nzdata_buf = new_nzdata_buf(h5dset->Rtype);
	if (nzdata_buf == NULL)
		return -1;

	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}

	/* Allocate 'chunkvp_buf', 'middlevp_buf', and 'destvp_buf'.
	   We set 'destvp_mode' to 1 because in the context of read_data_8()
	   we won't use 'destvp_buf.h5off' or 'destvp_buf.h5dim'. */
	if (alloc_chunkvp_middlevp_destvp_bufs(ndim,
		&chunkvp_buf, &middlevp_buf, &destvp_buf, 1) < 0)
	{
		H5Sclose(chunk_space_id);
		return -1;
	}

	chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	if (chunk_data_buf == NULL) {
		free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
						  &middlevp_buf,
						  &destvp_buf);
		H5Sclose(chunk_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}

	/* Walk over the chunks touched by the user selection. */
	tchunk_rank = 0;
	moved_along = ndim;
	do {
		update_chunkvp_destvp_bufs(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ret = read_data_from_chunk_8(h5dset,
			starts,
			nzindex_buf, nzdata_buf,
			inner_midx_buf->elts,
			&chunkvp_buf, &middlevp_buf, &destvp_buf,
			chunk_data_buf, chunk_space_id);
		if (ret < 0)
			break;
		tchunk_rank++;
		moved_along = _next_midx(ndim, ntchunks,
					 tchunk_midx_buf->elts);
	} while (moved_along < ndim);
	free(chunk_data_buf);
	free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
					  &middlevp_buf,
					  &destvp_buf);
	H5Sclose(chunk_space_id);

	if (ret < 0)
		return -1;

	/* Move the data in 'nzdata_buf' to an atomic vector. */
	SEXP nzdata = PROTECT(make_nzdata_from_buf(nzdata_buf, h5dset->Rtype));
	SET_VECTOR_ELT(ans, 1, nzdata);
	UNPROTECT(1);
	if (nzdata == R_NilValue)
		return -1;

	/* Move the data in 'nzindex_buf' to an ordinary matrix. */
	SEXP nzindex = PROTECT(make_nzindex_from_buf(nzindex_buf,
						     XLENGTH(nzdata),
						     ndim));
	SET_VECTOR_ELT(ans, 0, nzindex);
	UNPROTECT(1);
	if (nzindex == R_NilValue)
		return -1;

	return 0;
}


/****************************************************************************
 * _h5mread_sparse()
 */

static long long int set_ntchunks(const H5DSetDescriptor *h5dset,
				  const SEXP starts,
				  const LLongAEAE *tchunkidx_bufs,
				  int *ntchunks)
{
	int ndim, along, h5along, ntchunk;
	long long int total_num_tchunks;  /* total nb of touched chunks */
	SEXP start;

	ndim = h5dset->ndim;
	total_num_tchunks = 1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			ntchunk = LLongAE_get_nelt(tchunkidx_bufs->elts[along]);
		} else {
			ntchunk = h5dset->h5nchunk[h5along];
		}
		total_num_tchunks *= ntchunks[along] = ntchunk;
	}
	return total_num_tchunks;
}

/* Implements method 8.
   Return 'list(nzindex, nzdata, NULL)' or R_NilValue if an error occured. */
SEXP _h5mread_sparse(const H5DSetDescriptor *h5dset, SEXP starts, int *ans_dim)
{
	int ndim, ret;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;  /* touched chunk ids along each dim */
	IntAE *ntchunk_buf;  /* nb of touched chunks along each dim */
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
	set_ntchunks(h5dset, starts, tchunkidx_bufs, ntchunk_buf->elts);

	ans = PROTECT(NEW_LIST(3));
	ret = read_data_8(h5dset, starts,
			  breakpoint_bufs, tchunkidx_bufs,
			  ntchunk_buf->elts,
			  ans, ans_dim);
	UNPROTECT(1);  /* 'ans' */
	if (ret < 0)
		return R_NilValue;
	return ans;
}

