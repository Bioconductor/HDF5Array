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
	if (chunk_iter->chunk_space_id != -1)
		H5Sclose(chunk_iter->chunk_space_id);
	if (chunk_iter->chunk_data_buf != NULL)
		free(chunk_iter->chunk_data_buf);
	if (chunk_iter->compressed_chunk_data_buf != NULL)
		free(chunk_iter->compressed_chunk_data_buf);
	return;
}

int _init_ChunkIterator(ChunkIterator *chunk_iter,
		const H5DSetDescriptor *h5dset, SEXP index, int *selection_dim,
		int use_H5Dread_chunk)
{
	int ndim, ret;

	chunk_iter->h5dset = h5dset;
	chunk_iter->index = index;
	ndim = h5dset->ndim;

	/* Initialize the members that _destroy_ChunkIterator() will free
	   or close. */
	chunk_iter->tchunk_vp.h5off = NULL;
	chunk_iter->chunk_space_id = -1;
	chunk_iter->chunk_data_buf = NULL;
	chunk_iter->compressed_chunk_data_buf = NULL;

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

	/* Allocate members 'tchunk_vp', 'middle_vp', and 'dest_vp'.
	   We set 'dest_vp_mode' to ALLOC_OFF_AND_DIM because, in the context
	   of h5mread_sparse() or h5summarize(), we won't use 'dest_vp.h5off'
	   or 'dest_vp.h5dim', only 'dest_vp.off' and 'dest_vp.dim'. */
	ret = _alloc_tchunk_vp_middle_vp_dest_vp(ndim,
						 &chunk_iter->tchunk_vp,
						 &chunk_iter->middle_vp,
						 &chunk_iter->dest_vp,
						 ALLOC_OFF_AND_DIM);
	if (ret < 0)
		goto on_error;

	/* Set member 'tchunk_midx_buf'. */
	chunk_iter->tchunk_midx_buf = new_IntAE(ndim, ndim, 0)->elts;

	/* Set member 'chunk_space_id'. */
	chunk_iter->chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim,
						      NULL);
	if (chunk_iter->chunk_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		goto on_error;
	}

	/* Set member 'chunk_data_buf'. */
	chunk_iter->chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	if (chunk_iter->chunk_data_buf == NULL) {
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		goto on_error;
	}

	if (use_H5Dread_chunk) {
		/* Set member 'compressed_chunk_data_buf'. */
		chunk_iter->compressed_chunk_data_buf =
			malloc(h5dset->chunk_data_buf_size +
			       CHUNK_COMPRESSION_OVERHEAD);
		if (chunk_iter->compressed_chunk_data_buf == NULL) {
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for 'compressed_chunk_data_buf'");
			goto on_error;
		}
	}

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
	int moved_along, ret;

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
	if (chunk_iter->compressed_chunk_data_buf == NULL) {
		ret = _read_H5Viewport(chunk_iter->h5dset,
				&chunk_iter->tchunk_vp, &chunk_iter->middle_vp,
				chunk_iter->chunk_data_buf,
				chunk_iter->chunk_space_id);
	} else {
		/* Experimental! */
		ret = _read_h5chunk(chunk_iter->h5dset,
				&chunk_iter->tchunk_vp,
				chunk_iter->compressed_chunk_data_buf,
				chunk_iter->chunk_data_buf);
	}
	return ret < 0 ? ret : 1;
}

