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

int _alloc_H5Viewport(
	H5Viewport *vp,
	int ndim,
	int mode
);

void _free_H5Viewport(H5Viewport *vp);

int _alloc_h5chunkvp_middlevp_destvp_bufs(
	int ndim,
	H5Viewport *h5chunkvp_buf,
	H5Viewport *middlevp_buf,
	H5Viewport *destvp_buf,
	int destvp_mode
);

void _free_h5chunkvp_middlevp_destvp_bufs(
	H5Viewport *h5chunkvp_buf,
	H5Viewport *middlevp_buf,
	H5Viewport *destvp_buf
);

int _alloc_h5chunkvp_innervp_destvp_bufs(
	int ndim,
	H5Viewport *h5chunkvp_buf,
	H5Viewport *innervp_buf,
	H5Viewport *destvp_buf
);

void _free_h5chunkvp_innervp_destvp_bufs(
	H5Viewport *h5chunkvp_buf,
	H5Viewport *innervp_buf,
	H5Viewport *destvp_buf
);

int _map_starts_to_h5chunks(
	const H5DSetDescriptor *h5dset,
	SEXP starts,
	int *nstart_buf,
	IntAEAE *breakpoint_bufs,
	LLongAEAE *tchunkidx_bufs
);

int _select_H5Viewport(
	hid_t space_id,
	const H5Viewport *vp
);

int _add_H5Viewport_to_selection(
	hid_t space_id,
	const H5Viewport *vp
);

hid_t _create_mem_space(
	int ndim,
	const int *dim
);

int _read_H5Viewport(
	const H5DSetDescriptor *h5dset,
	const H5Viewport *dsetvp,
	const H5Viewport *memvp,
	void *mem,
	hid_t mem_space_id
);

void _update_h5chunkvp_destvp_bufs(
	const H5DSetDescriptor *h5dset,
	const int *tchunk_midx, int moved_along,
	SEXP starts,
	const IntAEAE *breakpoint_bufs,
	const LLongAEAE *tchunkidx_bufs,
	H5Viewport *h5chunkvp_buf,
	H5Viewport *destvp_buf
);

int _tchunk_is_fully_selected(
	int ndim,
	const H5Viewport *h5chunkvp,
	const H5Viewport *destvp
);

#define CHUNK_COMPRESSION_OVERHEAD 8  // empirical (increase if necessary)

int _read_h5chunk(
	const H5DSetDescriptor *h5dset,
	const H5Viewport *h5chunkvp,
	void *chunk_data_out,
	void *compressed_chunk_data_buf
);

#endif  /* _H5MREAD_HELPERS_H_ */

