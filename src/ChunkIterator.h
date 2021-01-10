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
	H5Viewport tchunk_vp;
	H5Viewport middle_vp;
	H5Viewport dest_vp;
	int *tchunk_midx_buf;
	int *inner_midx_buf;
	hid_t chunk_space_id;
	void *chunk_data_buf;
	int moved_along;
	long long int tchunk_rank;
} ChunkIterator;

void _destroy_ChunkIterator(
	ChunkIterator *chunk_iter
);

int _init_ChunkIterator(
	ChunkIterator *chunk_iter,
	const H5DSetDescriptor *h5dset,
	SEXP index
);

int _next_chunk(
	ChunkIterator *chunk_iter,
	int verbose
);

#endif  /* _CHUNKITERATOR_H_ */

