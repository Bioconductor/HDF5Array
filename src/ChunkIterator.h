#ifndef _CHUNKITERATOR_H_
#define _CHUNKITERATOR_H_

#include "H5DSetDescriptor.h"
#include "h5mread_helpers.h"
#include <Rdefines.h>
#include "S4Vectors_interface.h"

/* A data structure for iterating over the chunks of an HDF5 dataset. */
typedef struct chunk_iterator_t {
	const H5DSetDescriptor *h5dset;
	SEXP index;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;  /* touched chunk ids along each dim */
	int *num_tchunks;           /* nb of touched chunks along each dim */
	long long int total_num_tchunks;
	H5Viewport h5dset_vp, mem_vp;
	int *tchunk_midx_buf;
	int moved_along;
	long long int tchunk_rank;
} ChunkIterator;

/* A data structure for storing the data of a full chunk. */
typedef struct chunk_data_buffer_t {
	size_t data_length;
	hid_t data_type_id;
	size_t data_type_size;
	size_t data_size;
	hid_t data_space_id;
	void *data;
	H5Viewport data_vp;
	void *compressed_data;  /* experimental! */
} ChunkDataBuffer;

void _destroy_ChunkIterator(
	ChunkIterator *chunk_iter
);

int _init_ChunkIterator(
	ChunkIterator *chunk_iter,
	const H5DSetDescriptor *h5dset,
	SEXP index,
	int *selection_dim,
	int alloc_full_mem_vp
);

int _next_chunk(
	ChunkIterator *chunk_iter
);

void _print_tchunk_info(
	const ChunkIterator *chunk_iter
);

int _tchunk_is_truncated(
	const H5DSetDescriptor *h5dset,
	const H5Viewport *h5dset_vp
);

int _tchunk_is_fully_selected(
	int ndim,
	const H5Viewport *h5dset_vp,
	const H5Viewport *mem_vp
);

void _destroy_ChunkDataBuffer(
	ChunkDataBuffer *chunk_data_buf
);

int _init_ChunkDataBuffer(
	ChunkDataBuffer *chunk_data_buf,
	const H5DSetDescriptor *h5dset,
	int use_Rtype
);

int _load_chunk(
	const ChunkIterator *chunk_iter,
	ChunkDataBuffer *chunk_data_buf,
	int use_H5Dread_chunk
);

int _reclaim_vlen_bufs(
	ChunkDataBuffer *chunk_data_buf
);

#endif  /* _CHUNKITERATOR_H_ */

