/****************************************************************************
 *               Iterating over the chunks of an HDF5 dataset               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "ChunkIterator.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"

#include <stdlib.h>  /* for malloc, free */
#include <zlib.h>  /* for uncompress(), Z_OK, Z_MEM_ERROR, etc.. */


/****************************************************************************
 * Low-level helpers (non-exported)
 */

static int alloc_tchunk_vp_middle_vp_dest_vp(int ndim,
		H5Viewport *tchunk_vp,
		H5Viewport *middle_vp,
		H5Viewport *dest_vp, int dest_vp_mode)
{
	if (_alloc_H5Viewport(tchunk_vp, ndim, ALLOC_H5OFF_AND_H5DIM) < 0)
		return -1;
	middle_vp->h5off = _alloc_hsize_t_buf(ndim, 1, "'middle_vp->h5off'");
	if (middle_vp->h5off == NULL) {
		_free_H5Viewport(tchunk_vp);
		return -1;
	}
	middle_vp->h5dim = tchunk_vp->h5dim;
	if (_alloc_H5Viewport(dest_vp, ndim, dest_vp_mode) < 0) {
		free(middle_vp->h5off);
		_free_H5Viewport(tchunk_vp);
		return -1;
	}
	return 0;
}

static void free_tchunk_vp_middle_vp_dest_vp(
		H5Viewport *tchunk_vp,
		H5Viewport *middle_vp,
		H5Viewport *dest_vp)
{
	_free_H5Viewport(dest_vp);
	free(middle_vp->h5off);
	_free_H5Viewport(tchunk_vp);
	return;
}

static int map_starts_to_h5chunks(const H5DSetDescriptor *h5dset,
		SEXP starts,
		int *nstart_buf,
		IntAEAE *breakpoint_bufs, LLongAEAE *tchunkidx_bufs)
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

static long long int set_num_tchunks(const H5DSetDescriptor *h5dset,
		const SEXP starts,
		const LLongAEAE *tchunkidx_bufs,
		int *num_tchunks_buf)
{
	int ndim, along, h5along, n;
	long long int total_num_tchunks;  /* total nb of touched chunks */
	SEXP start;

	ndim = h5dset->ndim;
	total_num_tchunks = 1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			n = LLongAE_get_nelt(tchunkidx_bufs->elts[along]);
		} else {
			n = h5dset->h5nchunk[h5along];
		}
		total_num_tchunks *= num_tchunks_buf[along] = n;
	}
	return total_num_tchunks;
}

static void update_tchunk_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const LLongAEAE *tchunkidx_bufs,
		H5Viewport *tchunk_vp)
{
	int ndim, along, h5along, i;
	SEXP start;
	long long int tchunkidx;
	hsize_t chunkd, off, d;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = tchunk_midx[along];
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
		tchunk_vp->h5off[h5along] = off;
		tchunk_vp->h5dim[h5along] = d;
	}
	//printf("# tchunk_vp->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", tchunk_vp->h5off[h5along]);
	//printf("\n");
	//printf("# tchunk_vp->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", tchunk_vp->h5dim[h5along]);
	//printf("\n");
	return;
}

static void update_dest_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const IntAEAE *breakpoint_bufs,
		const H5Viewport *tchunk_vp, H5Viewport *dest_vp)
{
	int ndim, along, h5along, i, off, d;
	SEXP start;
	const int *breakpoint;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = tchunk_midx[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue ) {
			breakpoint = breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = tchunk_vp->h5off[h5along];
			d = tchunk_vp->h5dim[h5along];
		}
		if (dest_vp->h5off != NULL) {
			dest_vp->h5off[h5along] = off;
			dest_vp->h5dim[h5along] = d;
		}
		dest_vp->off[along] = off;
		dest_vp->dim[along] = d;
	}
	//printf("# dest_vp (offsets):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", dest_vp->off[along]);
	//printf("\n");
	//printf("# dest_vp (dims):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", dest_vp->dim[along]);
	//printf("\n");
	return;
}

static void update_tchunk_vp_dest_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		H5Viewport *tchunk_vp, H5Viewport *dest_vp)
{
	update_tchunk_vp(h5dset,
			tchunk_midx, moved_along,
			starts, tchunkidx_bufs,
			tchunk_vp);
	update_dest_vp(h5dset,
			tchunk_midx, moved_along,
			starts, breakpoint_bufs,
			tchunk_vp, dest_vp);
	return;
}


/****************************************************************************
 * read_h5chunk()
 */

static int uncompress_chunk_data(const void *compressed_chunk_data,
				 hsize_t compressed_size,
				 void *uncompressed_chunk_data,
				 size_t uncompressed_size)
{
	int ret;
	uLong destLen;

	destLen = (uLong) uncompressed_size;
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

static void transpose_bytes(const char *in, size_t nrow, size_t ncol, char *out)
{
	size_t i, j, in_offset;

	for (i = 0; i < nrow; i++) {
		in_offset = i;
		for (j = 0; j < ncol; j++) {
			*(out++) = *(in + in_offset);
			in_offset += nrow;
		}
	}
	return;
}

#define	CHUNK_COMPRESSION_OVERHEAD 8  // empirical (increase if necessary)

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
 *             ??
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
static int read_h5chunk(const H5DSetDescriptor *h5dset,
		const H5Viewport *h5chunk_vp,
		void *compressed_chunk_data_buf,
		void *chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(h5dset->dset_id,
					h5chunk_vp->h5off,
					&chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > h5dset->chunk_data_buf_size +
				 CHUNK_COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%lu + %d)",
				    chunk_storage_size,
				    h5dset->chunk_data_buf_size,
				    CHUNK_COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(h5dset->dset_id, H5P_DEFAULT,
			    h5chunk_vp->h5off, &filters,
			    compressed_chunk_data_buf);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dread_chunk() returned an error");
		return -1;
	}

	//printf("filters = %u\n", filters);

	//FIXME: This will error if chunk data is not compressed!
	//TODO: Decompress only if chunk data is compressed. There should be
	//a bit in the returned 'filters' that indicates this.
	ret = uncompress_chunk_data(compressed_chunk_data_buf,
				    chunk_storage_size,
				    chunk_data_buf,
				    h5dset->chunk_data_buf_size);
	if (ret < 0)
		return -1;
	size_t nval = h5dset->chunk_data_buf_size / h5dset->ans_elt_size;
	transpose_bytes(chunk_data_buf, nval, h5dset->ans_elt_size,
			compressed_chunk_data_buf);
	//print_chunk_data(h5dset, compressed_chunk_data_buf);
	return 0;
}


/****************************************************************************
 * Exported functions
 */

void _destroy_ChunkIterator(ChunkIterator *chunk_iter)
{
	if (chunk_iter->tchunk_vp.h5off != NULL)
		free_tchunk_vp_middle_vp_dest_vp(&chunk_iter->tchunk_vp,
						 &chunk_iter->middle_vp,
						 &chunk_iter->dest_vp);
	return;
}

int _init_ChunkIterator(ChunkIterator *chunk_iter,
		const H5DSetDescriptor *h5dset, SEXP index,
		int *selection_dim,
		int alloc_full_dest_vp)
{
	int ndim, ret;

	if (h5dset->h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("'h5dset->h5chunkdim' is NULL");
		return -1;
	}

	chunk_iter->h5dset = h5dset;
	chunk_iter->index = index;
	ndim = h5dset->ndim;

	/* Initialize the members that _destroy_ChunkIterator() will free
	   or close. */
	chunk_iter->tchunk_vp.h5off = NULL;

	/* Set members 'breakpoint_bufs' and 'tchunkidx_bufs'.
	   Also populate 'selection_dim' if not set to NULL. */
	chunk_iter->breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunk_iter->tchunkidx_bufs = new_LLongAEAE(ndim, ndim);
	ret = map_starts_to_h5chunks(h5dset, index, selection_dim,
				     chunk_iter->breakpoint_bufs,
				     chunk_iter->tchunkidx_bufs);
	if (ret < 0)
		goto on_error;

	/* Set members 'num_tchunks' and 'total_num_tchunks'. */
	chunk_iter->num_tchunks = new_IntAE(ndim, ndim, 0)->elts;
	chunk_iter->total_num_tchunks = set_num_tchunks(h5dset, index,
						chunk_iter->tchunkidx_bufs,
						chunk_iter->num_tchunks);

	/* Allocate members 'tchunk_vp', 'middle_vp', and 'dest_vp'. */
	ret = alloc_tchunk_vp_middle_vp_dest_vp(ndim,
				&chunk_iter->tchunk_vp,
				&chunk_iter->middle_vp,
				&chunk_iter->dest_vp,
				alloc_full_dest_vp ? ALLOC_ALL_FIELDS
						   : ALLOC_OFF_AND_DIM);
	if (ret < 0)
		goto on_error;

	/* Set member 'tchunk_midx_buf'. */
	chunk_iter->tchunk_midx_buf = new_IntAE(ndim, ndim, 0)->elts;

	/* Set member 'tchunk_rank'. */
	chunk_iter->tchunk_rank = -1;
	return 0;

    on_error:
	_destroy_ChunkIterator(chunk_iter);
	return -1;
}

/* Return:
 *     1 = if the chunk before the move was not the last chunk and the move to
 *         the next chunk was successful;
 *     0 = if the chunk before the move was the last chunk and so the move to
 *         the next chunk was not possible;
 *   < 0 = if error
 * Typical use:
 *     while (ret = _next_chunk(&chunk_iter)) {
 *         if (ret < 0) {
 *             an error occured
 *         }
 *         handle current chunk
 *     }
 */
int _next_chunk(ChunkIterator *chunk_iter)
{
	const H5DSetDescriptor *h5dset;

	chunk_iter->tchunk_rank++;
	if (chunk_iter->tchunk_rank == chunk_iter->total_num_tchunks)
		return 0;
	h5dset = chunk_iter->h5dset;
	chunk_iter->moved_along = chunk_iter->tchunk_rank == 0 ?
					h5dset->ndim :
					_next_midx(h5dset->ndim,
						chunk_iter->num_tchunks,
						chunk_iter->tchunk_midx_buf);
	update_tchunk_vp_dest_vp(h5dset,
			chunk_iter->tchunk_midx_buf,
			chunk_iter->moved_along,
			chunk_iter->index,
			chunk_iter->breakpoint_bufs, chunk_iter->tchunkidx_bufs,
			&chunk_iter->tchunk_vp, &chunk_iter->dest_vp);
	return 1;
}

void _print_tchunk_info(int ndim,
		const int *num_tchunks_buf, const int *tchunk_midx,
		int tchunk_rank,
		const SEXP index, const LLongAEAE *tchunkidx_bufs,
		const H5Viewport *tchunk_vp)
{
	int along, h5along, i;
	long long int total_num_tchunks, tchunkidx;

	total_num_tchunks = 1;
	for (along = 0; along < ndim; along++)
		total_num_tchunks *= num_tchunks_buf[along];

	printf("processing chunk %d/%lld: [",
	       tchunk_rank + 1, total_num_tchunks);
	for (along = 0; along < ndim; along++) {
		i = tchunk_midx[along] + 1;
		if (along != 0)
			printf(", ");
		printf("%d/%d", i, num_tchunks_buf[along]);
	}
	printf("] -- <<");
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		i = tchunk_midx[along];
		if (GET_LIST_ELT(index, along) != R_NilValue) {
			tchunkidx = tchunkidx_bufs->elts[along]->elts[i];
		} else {
			tchunkidx = i;
		}
		if (along != 0)
			printf(", ");
		printf("#%lld=%llu:%llu", tchunkidx + 1,
		       tchunk_vp->h5off[h5along] + 1,
		       tchunk_vp->h5off[h5along] + tchunk_vp->h5dim[h5along]);
	}
	printf(">>\n");
	return;
}

/* Return 1 if the chunk that 'tchunk_vp' is pointing at is "truncated"
   (a.k.a. "partial edge chunk" in HDF5's terminology), and 0 otherwise
   (i.e. if the new chunk is a full-size chunk). */
int _tchunk_is_truncated(const H5DSetDescriptor *h5dset,
			 const H5Viewport *tchunk_vp)
{
	int ndim, h5along;
	hsize_t chunkd, d;

	ndim = h5dset->ndim;
	for (h5along = 0; h5along < ndim; h5along++) {
		chunkd = h5dset->h5chunkdim[h5along];
		d = tchunk_vp->h5dim[h5along];
		if (d != chunkd)
			return 1;
	}
	return 0;
}

int _tchunk_is_fully_selected(int ndim,
		const H5Viewport *tchunk_vp,
		const H5Viewport *dest_vp)
{
	int along, h5along, not_fully;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		not_fully = tchunk_vp->h5dim[h5along] !=
			    (hsize_t) dest_vp->dim[along];
		if (not_fully)
			return 0;
	}
	return 1;
}

void _destroy_ChunkDataBuffer(ChunkDataBuffer *chunk_data_buf)
{
	if (chunk_data_buf->data_space_id != -1)
		H5Sclose(chunk_data_buf->data_space_id);
	if (chunk_data_buf->data != NULL)
		free(chunk_data_buf->data);
	if (chunk_data_buf->compressed_data != NULL)
		free(chunk_data_buf->compressed_data);
	return;
}

int _init_ChunkDataBuffer(ChunkDataBuffer *chunk_data_buf,
		const H5DSetDescriptor *h5dset)
{
	size_t data_length;
	int h5along;

	if (h5dset->h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("'h5dset->h5chunkdim' is NULL");
		return -1;
	}

	/* Initialize the members that _destroy_ChunkDataBuffer() will free
	   or close. */
	chunk_data_buf->data_space_id = -1;
	chunk_data_buf->data = NULL;
	chunk_data_buf->compressed_data = NULL;

	/* Set members 'data_length' and 'data_size'. */
	data_length = 1;
	for (h5along = 0; h5along < h5dset->ndim; h5along++)
		data_length *= h5dset->h5chunkdim[h5along];
	chunk_data_buf->data_length = data_length;
	chunk_data_buf->data_size = data_length * h5dset->ans_elt_size;

	return 0;
}

int _load_chunk(const ChunkIterator *chunk_iter,
		ChunkDataBuffer *chunk_data_buf,
		int use_H5Dread_chunk)
{
	const H5DSetDescriptor *h5dset;
	hid_t data_space_id;
	int ret;

	if (chunk_data_buf->data == NULL) {
		chunk_data_buf->data = malloc(chunk_data_buf->data_size);
		if (chunk_data_buf->data == NULL) {
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for 'chunk_data_buf->data'");
			return -1;
		}
	}
	h5dset = chunk_iter->h5dset;
	if (!use_H5Dread_chunk) {
		if (chunk_data_buf->data_space_id == -1) {
			data_space_id = H5Screate_simple(h5dset->ndim,
							 h5dset->h5chunkdim,
							 NULL);
			if (data_space_id < 0) {
				PRINT_TO_ERRMSG_BUF("H5Screate_simple() "
						    "returned an error");
				return -1;
			}
			chunk_data_buf->data_space_id = data_space_id;
		}
		ret = _read_H5Viewport(h5dset,
				&chunk_iter->tchunk_vp,
				&chunk_iter->middle_vp,
				chunk_data_buf->data,
				chunk_data_buf->data_space_id);
	} else {
		/* Experimental! */
		if (chunk_data_buf->compressed_data == NULL) {
			chunk_data_buf->compressed_data =
				malloc(chunk_data_buf->data_size +
				       CHUNK_COMPRESSION_OVERHEAD);
			if (chunk_data_buf->compressed_data == NULL) {
				PRINT_TO_ERRMSG_BUF(
					"failed to allocate memory for "
					"'chunk_data_buf->compressed_data'");
				return -1;
			}
		}
		ret = read_h5chunk(h5dset,
				&chunk_iter->tchunk_vp,
				chunk_data_buf->compressed_data,
				chunk_data_buf->data);
	}
	return ret;
}

