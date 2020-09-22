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

#endif  /* _H5MREAD_HELPERS_H_ */

