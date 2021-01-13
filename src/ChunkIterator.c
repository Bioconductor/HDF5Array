/****************************************************************************
 *               Iterating over the chunks of an HDF5 dataset               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "ChunkIterator.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"

#include <stdlib.h>  /* for malloc, free */

void _destroy_ChunkIterator(ChunkIterator *chunk_iter)
{
	if (chunk_iter->tchunk_vp.h5off != NULL)
		_free_tchunk_vp_middle_vp_dest_vp(&chunk_iter->tchunk_vp,
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
	ret = _map_starts_to_h5chunks(h5dset, index, selection_dim,
				      chunk_iter->breakpoint_bufs,
				      chunk_iter->tchunkidx_bufs);
	if (ret < 0)
		goto on_error;

	/* Set members 'num_tchunks' and 'total_num_tchunks'. */
	chunk_iter->num_tchunks = new_IntAE(ndim, ndim, 0)->elts;
	chunk_iter->total_num_tchunks = _set_num_tchunks(h5dset, index,
						chunk_iter->tchunkidx_bufs,
						chunk_iter->num_tchunks);

	/* Allocate members 'tchunk_vp', 'middle_vp', and 'dest_vp'. */
	ret = _alloc_tchunk_vp_middle_vp_dest_vp(ndim,
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
	int moved_along;

	chunk_iter->tchunk_rank++;
	if (chunk_iter->tchunk_rank == chunk_iter->total_num_tchunks)
		return 0;
	moved_along = chunk_iter->tchunk_rank == 0 ?
				chunk_iter->h5dset->ndim :
				_next_midx(chunk_iter->h5dset->ndim,
					   chunk_iter->num_tchunks,
					   chunk_iter->tchunk_midx_buf);
	_update_tchunk_vp_dest_vp(chunk_iter->h5dset,
			chunk_iter->tchunk_midx_buf,
			moved_along,
			chunk_iter->index,
			chunk_iter->breakpoint_bufs, chunk_iter->tchunkidx_bufs,
			&chunk_iter->tchunk_vp, &chunk_iter->dest_vp);
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
		ret = _read_h5chunk(h5dset,
				&chunk_iter->tchunk_vp,
				chunk_data_buf->compressed_data,
				chunk_data_buf->data);
	}
	return ret;
}

