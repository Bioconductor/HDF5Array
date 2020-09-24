#ifndef _H5MREAD_HELPERS_H_
#define _H5MREAD_HELPERS_H_

#include "H5DSetDescriptor.h"
#include "hdf5.h"

static inline int _next_midx(int ndim, const int *max_idx_plus_one,
			     int *midx_buf)
{
	int along, i;

	for (along = 0; along < ndim; along++) {
		i = midx_buf[along] + 1;
		if (i < max_idx_plus_one[along]) {
			midx_buf[along] = i;
			break;
		}
		midx_buf[along] = 0;
	}
	return along;
}

/* A data structure for representing a viewport on a HDF5 dataset. */
typedef struct {
	hsize_t *h5off, *h5dim;
	int *off, *dim; // same as h5off and h5dim but stored as int
} H5Viewport;

#define	ALLOC_ALL_FIELDS		0
#define	ALLOC_H5OFF_AND_H5DIM		1
#define	ALLOC_OFF_AND_DIM		2

int _alloc_H5Viewport(
	H5Viewport *vp,
	int ndim,
	int mode
);

void _free_H5Viewport(H5Viewport *vp);

int _alloc_h5chunk_vp_middle_vp_dest_vp(
	int ndim,
	H5Viewport *h5chunk_vp,
	H5Viewport *middle_vp,
	H5Viewport *dest_vp,
	int dest_vp_mode
);

void _free_h5chunk_vp_middle_vp_dest_vp(
	H5Viewport *h5chunk_vp,
	H5Viewport *middle_vp,
	H5Viewport *dest_vp
);

int _alloc_h5chunk_vp_inner_vp_dest_vp(
	int ndim,
	H5Viewport *h5chunk_vp,
	H5Viewport *inner_vp,
	H5Viewport *dest_vp
);

void _free_h5chunk_vp_inner_vp_dest_vp(
	H5Viewport *h5chunk_vp,
	H5Viewport *inner_vp,
	H5Viewport *dest_vp
);

int _map_starts_to_h5chunks(
	const H5DSetDescriptor *h5dset,
	SEXP starts,
	int *nstart_buf,
	IntAEAE *breakpoint_bufs,
	LLongAEAE *tchunkidx_bufs
);

long long int _set_num_tchunks(
	const H5DSetDescriptor *h5dset,
	const SEXP starts,
	const LLongAEAE *tchunkidx_bufs,
	int *num_tchunks_buf
);

hid_t _create_mem_space(
	int ndim,
	const int *dim
);

int _select_H5Viewport(
	hid_t space_id,
	const H5Viewport *vp
);

int _add_H5Viewport_to_selection(
	hid_t space_id,
	const H5Viewport *vp
);

int _read_H5Viewport(
	const H5DSetDescriptor *h5dset,
	const H5Viewport *h5dset_vp,
	const H5Viewport *mem_vp,
	void *mem,
	hid_t mem_space_id
);

void _update_h5chunk_vp_dest_vp(
	const H5DSetDescriptor *h5dset,
	const int *tchunk_midx, int moved_along,
	SEXP starts,
	const IntAEAE *breakpoint_bufs,
	const LLongAEAE *tchunkidx_bufs,
	H5Viewport *h5chunk_vp,
	H5Viewport *dest_vp
);

int _tchunk_is_fully_selected(
	int ndim,
	const H5Viewport *h5chunk_vp,
	const H5Viewport *dest_vp
);

#define CHUNK_COMPRESSION_OVERHEAD 8  // empirical (increase if necessary)

int _read_h5chunk(
	const H5DSetDescriptor *h5dset,
	const H5Viewport *h5chunk_vp,
	void *chunk_data_out,
	void *compressed_chunk_data_buf
);

#endif  /* _H5MREAD_HELPERS_H_ */

