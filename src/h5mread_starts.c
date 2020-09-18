/****************************************************************************
 *               Workhorses behind h5mread methods 4, 5, 6, 7               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_starts.h"

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
#include <zlib.h>
//#include <time.h>


/****************************************************************************
 * Low-level helpers
 */

/* Used in read_data_4_5() and read_data_7(). */
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

/* Used in read_data_6(). */
static int alloc_chunkvp_innervp_destvp_bufs(int ndim,
			H5Viewport *chunkvp_buf,
			H5Viewport *innervp_buf,
			H5Viewport *destvp_buf)
{
	if (_alloc_H5Viewport(chunkvp_buf, ndim, 0) < 0)
		return -1;
	if (_alloc_H5Viewport(innervp_buf, ndim, 0) < 0) {
		_free_H5Viewport(chunkvp_buf);
		return -1;
	}
	if (_alloc_H5Viewport(destvp_buf, ndim, 2) < 0) {
		_free_H5Viewport(innervp_buf);
		_free_H5Viewport(chunkvp_buf);
		return -1;
	}
	return 0;
}

static void free_chunkvp_innervp_destvp_bufs(
			H5Viewport *chunkvp_buf,
			H5Viewport *innervp_buf,
			H5Viewport *destvp_buf)
{
	_free_H5Viewport(destvp_buf);
	_free_H5Viewport(innervp_buf);
	_free_H5Viewport(chunkvp_buf);
	return;
}

static int read_H5Viewport(const H5DSetDescriptor *h5dset,
			   const H5Viewport *dsetvp,
			   const H5Viewport *memvp,
			   void *mem, hid_t mem_space_id)
{
	int ret;

	ret = _select_H5Viewport(h5dset->space_id, dsetvp);
	if (ret < 0)
		return -1;
	ret = _select_H5Viewport(mem_space_id, memvp);
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
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
		      const H5Viewport *chunkvp,
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

static void init_in_offset_and_out_offset(int ndim, SEXP starts,
			const int *out_dim, const H5Viewport *destvp,
			const H5Viewport *chunkvp,
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

static inline void update_in_offset_and_out_offset(int ndim,
			const int *inner_midx, int inner_moved_along,
			SEXP starts,
			const int *out_dim, const H5Viewport *destvp,
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

static int load_val_to_nzdata(const H5DSetDescriptor *h5dset,
		const void *in, size_t in_offset,
		SEXP nzdata, R_xlen_t nzdata_offset)
{
	switch (h5dset->Rtype) {
	    case LGLSXP: {
		int val = ((int *) in)[in_offset];
		if (val == 0)
			return 0;
		if (nzdata != R_NilValue)
			LOGICAL(nzdata)[nzdata_offset] = val;
	    } break;
	    case INTSXP: {
		int val = ((int *) in)[in_offset];
		if (val == 0)
			return 0;
		if (nzdata != R_NilValue)
			INTEGER(nzdata)[nzdata_offset] = val;
	    } break;
	    case REALSXP: {
		double val = ((double *) in)[in_offset];
		if (val == 0.0)
			return 0;
		if (nzdata != R_NilValue)
			REAL(nzdata)[nzdata_offset] = val;
	    } break;
	    case STRSXP: {
		const char *val = ((char *) in) + in_offset * h5dset->H5size;
		int val_len;
		for (val_len = 0; val_len < h5dset->H5size; val_len++)
			if (val[val_len] == 0)
				break;
		if (val_len == 0)
			return 0;
		if (nzdata != R_NilValue) {
			int is_na = h5dset->as_na_attr &&
				    val_len == 2 &&
				    val[0] == 'N' && val[1] == 'A';
			if (is_na) {
			    SET_STRING_ELT(nzdata, nzdata_offset, NA_STRING);
			} else {
			    SEXP nzdata_elt = PROTECT(mkCharLen(val, val_len));
			    SET_STRING_ELT(nzdata, nzdata_offset, nzdata_elt);
			    UNPROTECT(1);
			}
		}
	    } break;
	    case RAWSXP: {
		char val = ((char *) in)[in_offset];
		if (val == 0)
			return 0;
		if (nzdata != R_NilValue)
			RAW(nzdata)[nzdata_offset] = val;
	    } break;
	    default:
		PRINT_TO_ERRMSG_BUF("unsupported type: %s",
				    CHAR(type2str(h5dset->Rtype)));
		return -1;
	}
	return 1;
}

static void load_midx_to_nzindex(const H5Viewport *destvp,
		const int *inner_midx_buf,
		SEXP nzindex, R_xlen_t nzindex_nrow, int ndim,
		R_xlen_t nzdata_offset)
{
	int along;

	for (along = 0; along < ndim; along++) {
		INTEGER(nzindex)[nzdata_offset + nzindex_nrow * along] =
			destvp->off[along] + inner_midx_buf[along] + 1;
	}
	return;
}

static int load_nonzero_val_to_nzdata_buf(
		const H5DSetDescriptor *h5dset,
		const int *in, size_t in_offset,
		void *nzdata_buf)
{
	size_t nzdata_len;

	/* The maximum 'nzdata' length that we support is INT_MAX. This
	   is because at the end of the read_data_4_5() call we will need
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

static void load_midx_to_nzindex_buf(int ndim, const H5Viewport *destvp,
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

/* 2-pass, and one call to NEW_LIST(), allocMatrix(), and allocVector() per
   touched chunk --> no very efficient! */
static int gather_selected_chunk_data_as_sparse(
		const H5DSetDescriptor *h5dset,
		SEXP starts, const void *in, const H5Viewport *chunkvp,
		SEXP ans, long long int current_chunk_Lidx,
		const H5Viewport *destvp, int *inner_midx_buf)
{
	int ndim, inner_moved_along, ret;
	size_t in_offset;
	R_xlen_t nzdata_len, nzdata_offset;
	SEXP chunk_data, nzindex, nzdata;

	ndim = h5dset->ndim;

	/* Walk on the selected elements in current chunk and count the
	   non-zero ones. */
	init_in_offset(ndim, starts, destvp, chunkvp,
		       h5dset->h5chunkdim, &in_offset);
	nzdata_len = 0;
	while (1) {
		ret = load_val_to_nzdata(h5dset, in, in_offset, R_NilValue, 0);
		if (ret < 0)
			return ret;
		if (ret > 0)
			nzdata_len++;
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

	/* Alloc 'chunk_data', 'nzindex' and 'nzdata'. */
	chunk_data = PROTECT(NEW_LIST(2));
	nzindex = PROTECT(allocMatrix(INTSXP, (int) nzdata_len, ndim));
	nzdata = PROTECT(allocVector(h5dset->Rtype, nzdata_len));
	SET_VECTOR_ELT(chunk_data, 0, nzindex);
	SET_VECTOR_ELT(chunk_data, 1, nzdata);
	SET_VECTOR_ELT(VECTOR_ELT(ans, 1), current_chunk_Lidx, chunk_data);
	UNPROTECT(3);

	/* Walk again on the selected elements in current chunk and load the
	   non-zero ones in 'chunk_data'. */
	init_in_offset(ndim, starts, destvp, chunkvp,
		       h5dset->h5chunkdim, &in_offset);
	nzdata_offset = 0;
	while (1) {
		ret = load_val_to_nzdata(h5dset, in, in_offset,
					 nzdata, nzdata_offset);
		if (ret < 0)
			return ret;
		if (ret > 0) {
			load_midx_to_nzindex(destvp, inner_midx_buf,
					     nzindex, nzdata_len, ndim,
					     nzdata_offset);
			nzdata_offset++;
		}
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

/* Single pass */
static int gather_selected_chunk_data_as_sparse_in_one_pass(
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

static int gather_selected_chunk_data_as_dense(
		const H5DSetDescriptor *h5dset,
		SEXP starts, const void *in, const H5Viewport *chunkvp,
		SEXP ans, const int *out_dim,
		const H5Viewport *destvp, int *inner_midx_buf)
{
	int ndim, inner_moved_along, ret;
	size_t in_offset, out_offset;
	long long int num_elts;

	ndim = h5dset->ndim;

	/* Walk on the selected elements in current chunk. */
	init_in_offset_and_out_offset(ndim, starts,
			out_dim, destvp,
			chunkvp, h5dset->h5chunkdim,
			&in_offset, &out_offset);
	num_elts = 0;
	while (1) {
		ret = load_val_to_array(h5dset, in, in_offset, ans, out_offset);
		if (ret < 0)
			return ret;
		num_elts++;
		inner_moved_along = _next_midx(ndim, destvp->dim,
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
			const LLongAEAE *tchunkidx_bufs,
			int sparse, long long int current_chunk_Lidx,
			SEXP ans, const int *ans_dim,
			IntAE *nzindex_buf, void *nzdata_buf,
			int *inner_midx_buf,
			H5Viewport *chunkvp_buf,
			const H5Viewport *middlevp,
			H5Viewport *destvp_buf,
			void *chunk_data_buf, hid_t chunk_space_id,
			void *compressed_chunk_data_buf)
{
	int ret;

	if (method == 4) {
		/* It takes about 218s on my laptop to load all the chunks
		   from the EH1040 dataset (big 10x Genomics brain dataset
		   in dense format, chunks of 100x100, wrapped in the
		   TENxBrainData package). That's 60 microseconds per chunk! */
		ret = read_H5Viewport(h5dset,
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
	if (sparse) {
/*
		ret = gather_selected_chunk_data_as_sparse(
				h5dset,
				starts, chunk_data_buf, chunkvp_buf,
				ans, current_chunk_Lidx,
				destvp_buf, inner_midx_buf);
*/
		ret = gather_selected_chunk_data_as_sparse_in_one_pass(
				h5dset,
				starts, chunk_data_buf, chunkvp_buf,
				nzindex_buf, nzdata_buf,
				destvp_buf, inner_midx_buf);
	} else {
		ret = gather_selected_chunk_data_as_dense(
				h5dset,
				starts, chunk_data_buf, chunkvp_buf,
				ans, ans_dim,
				destvp_buf, inner_midx_buf);
	}
	return ret;
}

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
			 const int *ntchunks,
			 int sparse, SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ret;
	hid_t chunk_space_id;
	void *chunk_data_buf, *compressed_chunk_data_buf = NULL;
	H5Viewport chunkvp_buf, middlevp_buf, destvp_buf;
	IntAE *tchunk_midx_buf, *inner_midx_buf, *nzindex_buf;
	void *nzdata_buf;
	long long int current_chunk_Lidx;

	ndim = h5dset->ndim;

	/* Prepare buffers. */

	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	if (sparse) {
		nzindex_buf = new_IntAE(ndim, 0, 0);
		nzdata_buf = new_nzdata_buf(h5dset->Rtype);
		if (nzdata_buf == NULL)
			return -1;
	} else {
		nzindex_buf = nzdata_buf = NULL;
	}

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

	/* Walk over the chunks touched by the user selection. */
	current_chunk_Lidx = 0;
	moved_along = ndim;
	do {
		update_chunkvp_destvp_bufs(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ret = read_data_from_chunk_4_5(h5dset, method,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			sparse, current_chunk_Lidx,
			ans, ans_dim,
			nzindex_buf, nzdata_buf,
			inner_midx_buf->elts,
			&chunkvp_buf, &middlevp_buf, &destvp_buf,
			chunk_data_buf, chunk_space_id,
			compressed_chunk_data_buf);
		if (ret < 0)
			break;
		current_chunk_Lidx++;
		moved_along = _next_midx(ndim, ntchunks,
					 tchunk_midx_buf->elts);
	} while (moved_along < ndim);
	free(chunk_data_buf);
	free_chunkvp_middlevp_destvp_bufs(&chunkvp_buf,
					  &middlevp_buf,
					  &destvp_buf);
	H5Sclose(chunk_space_id);
	if (sparse) {
		/* Move the data in 'nzdata_buf' to an atomic vector. */
		SEXP nzdata = PROTECT(make_nzdata_from_buf(nzdata_buf,
							   h5dset->Rtype));
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
	}
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
		const H5Viewport *destvp,
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
			SEXP starts, const H5Viewport *chunkvp,
			H5Viewport *innervp_buf)
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
			SEXP starts, const H5Viewport *destvp,
			const int *inner_midx, int inner_moved_along,
			const IntAEAE *inner_breakpoint_bufs,
			H5Viewport *innervp_buf)
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
		const H5Viewport *destvp,
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
		inner_moved_along = _next_midx(ndim, destvp->dim,
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
		const H5Viewport *destvp, const H5Viewport *chunkvp,
		const IntAEAE *inner_breakpoint_bufs,
		const int *inner_nblock,
		int *inner_midx_buf,
		H5Viewport *innervp_buf)
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
		ret = _add_H5Viewport_to_selection(h5dset->space_id,
						   innervp_buf);
		if (ret < 0)
			return -1;
		inner_moved_along = _next_midx(ndim, inner_nblock,
					       inner_midx_buf);
	} while (inner_moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(inner_nblock)

	return num_hyperslabs;
}

static int read_selection(const H5DSetDescriptor *h5dset,
			  const H5Viewport *destvp_buf,
			  void *mem, hid_t mem_space_id)
{
	int ret;

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

static int read_data_from_chunk_6(const H5DSetDescriptor *h5dset,
			const int *chunk_midx, int moved_along,
			SEXP starts,
			const IntAEAE *breakpoint_bufs,
			const LLongAEAE *tchunkidx_bufs,
			void *mem, hid_t mem_space_id,
			int *inner_midx_buf,
			H5Viewport *chunkvp_buf,
			H5Viewport *innervp_buf,
			H5Viewport *destvp_buf,
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
	//printf("# chunk %lld: %s\n", current_chunk_Lidx,
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
		       const LLongAEAE *tchunkidx_bufs,
		       const int *ntchunks,
		       SEXP ans, const int *ans_dim)
{
	void *mem;
	int ndim, moved_along, ret;
	hid_t mem_space_id;
	H5Viewport chunkvp_buf, innervp_buf, destvp_buf;
	//hsize_t *coord_buf;
	IntAE *tchunk_midx_buf, *inner_nblock_buf, *inner_midx_buf;
	IntAEAE *inner_breakpoint_bufs;
	long long int current_chunk_Lidx;

	mem = DATAPTR(ans);
	if (mem == NULL)
		return -1;
	ndim = h5dset->ndim;
	mem_space_id = _create_mem_space(ndim, ans_dim);
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
	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_nblock_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	/* Walk over the chunks touched by the user selection. */
	current_chunk_Lidx = 0;
	moved_along = ndim;
	//clock_t t_select_elements = 0, t_read_selection = 0, t0;
	do {
		update_chunkvp_destvp_bufs(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ret = read_data_from_chunk_6(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			mem, mem_space_id,
			inner_midx_buf->elts,
			&chunkvp_buf, &innervp_buf, &destvp_buf,
			inner_breakpoint_bufs, inner_nblock_buf);
			//coord_buf);
		if (ret < 0)
			break;
		current_chunk_Lidx++;
		moved_along = _next_midx(ndim, ntchunks, tchunk_midx_buf->elts);
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
				   const H5Viewport *chunkvp,
				   const H5Viewport *destvp)
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
		       const LLongAEAE *tchunkidx_bufs,
		       const int *ntchunks,
		       SEXP ans, const int *ans_dim)
{
	int ndim, moved_along, ok, ret;
	hid_t chunk_space_id, mem_space_id;
	void *mem, *chunk_data_buf;
	H5Viewport chunkvp_buf, middlevp_buf, destvp_buf;
	IntAE *tchunk_midx_buf, *inner_midx_buf;
	long long int current_chunk_Lidx;

	ndim = h5dset->ndim;
	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}
	mem = DATAPTR(ans);
	if (mem == NULL) {
		H5Sclose(chunk_space_id);
		return -1;
	}
	mem_space_id = _create_mem_space(ndim, ans_dim);
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
	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	/* Walk over the chunks touched by the user selection. */
	current_chunk_Lidx = 0;
	moved_along = ndim;
	do {
		update_chunkvp_destvp_bufs(h5dset,
			tchunk_midx_buf->elts, moved_along,
			starts, breakpoint_bufs, tchunkidx_bufs,
			&chunkvp_buf, &destvp_buf);
		ok = chunk_is_fully_selected(h5dset->ndim,
					     &chunkvp_buf, &destvp_buf);
		if (ok) {
			/* Load the chunk **directly** into 'ans' (no
			   intermediate buffer). */
			ret = read_H5Viewport(h5dset,
				&chunkvp_buf, &destvp_buf,
				mem, mem_space_id);
		} else {
			/* Load the **entire** chunk to an intermediate
			   buffer then copy the user-selected data from
			   the intermediate buffer to 'ans'. */
			ret = read_data_from_chunk_4_5(h5dset, 4,
				tchunk_midx_buf->elts, moved_along,
				starts, breakpoint_bufs, tchunkidx_bufs,
				0, current_chunk_Lidx,
				ans, ans_dim,
				NULL, NULL,
				inner_midx_buf->elts,
				&chunkvp_buf, &middlevp_buf, &destvp_buf,
				chunk_data_buf, chunk_space_id, NULL);
		}
		if (ret < 0)
			break;
		current_chunk_Lidx++;
		moved_along = _next_midx(ndim, ntchunks, tchunk_midx_buf->elts);
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
 * _h5mread_starts()
 */

static int map_starts_to_chunks(const H5DSetDescriptor *h5dset,
				SEXP starts,
				int *nstart_buf,
				IntAEAE *breakpoint_bufs,
				LLongAEAE *tchunkidx_bufs)
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
				     breakpoint_bufs, tchunkidx_bufs);
}

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

/*
 * Return an ordinary array if 'as.sparse' is FALSE, or a list of length 3
 * if it's TRUE. If the latter:
 *   - The first list element in 'ans' will be an atomic vector of length 0
 *     of type 'h5dset.Rtype'.
 *   - The second list element in 'ans' will be a list of length the total nb
 *     of touched chunks where each list element is itself a list of length 2.
 *     This list of length 2 is a sparse representation of the data loaded
 *     from the corresponding touched chunk i.e. its 1st element is the
 *     "nzindex" (matrix) and its 2nd element the "nzdata" (atomic vector).
 *   - The third list element in 'ans' will be NULL.
 * Return R_NilValue on error.
 */
SEXP _h5mread_starts(const H5DSetDescriptor *h5dset, SEXP starts,
		     int sparse, int method, int *ans_dim)
{
	int ndim, ret, along;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;  /* touched chunk ids along each dim */
	IntAE *ntchunk_buf;  /* nb of touched chunks along each dim */
	long long int total_num_tchunks;  /* total nb of touched chunks */
	R_xlen_t ans_len;
	SEXP ans;

	ndim = h5dset->ndim;

	/* This call will populate 'ans_dim', 'breakpoint_bufs',
	   and 'tchunkidx_bufs'. */
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	tchunkidx_bufs = new_LLongAEAE(ndim, ndim);
	ret = map_starts_to_chunks(h5dset, starts,
				   ans_dim, breakpoint_bufs, tchunkidx_bufs);
	if (ret < 0)
		return R_NilValue;

	ntchunk_buf = new_IntAE(ndim, ndim, 0);
	total_num_tchunks = set_ntchunks(h5dset, starts,
					 tchunkidx_bufs, ntchunk_buf->elts);

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	if (sparse) {
/*
		ans = PROTECT(NEW_LIST(3));
		SEXP nzdata0 = PROTECT(allocVector(h5dset->Rtype, 0));
		SET_VECTOR_ELT(ans, 0, nzdata0);
		SEXP sparse_data_list =
			PROTECT(NEW_LIST((R_xlen_t) total_num_tchunks));
		SET_VECTOR_ELT(ans, 1, sparse_data_list);
		UNPROTECT(2);
*/
		/* Will become 'list(nzindex, nzdata, ans_dim)'. */
		ans = PROTECT(NEW_LIST(3));
	} else {
		ans = PROTECT(allocVector(h5dset->Rtype, ans_len));
	}

	/* Note that ans_len == 0 is equivalent to total_num_tchunks == 0.
	   Both mean that the selection is empty i.e. that at least one of
	   its dimensions is 0. */
	if (ans_len != 0) {
		if (method <= 5) {
			/* methods 4 and 5 */
			ret = read_data_4_5(h5dset, method, starts,
					breakpoint_bufs, tchunkidx_bufs,
					ntchunk_buf->elts,
					sparse, ans, ans_dim);
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

