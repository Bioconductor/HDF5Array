/****************************************************************************
 *  Some low-level helper functions to support the various h5mread methods  *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_helpers.h"

/*
Some useful links:
- Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
    https://support.hdfgroup.org/HDF5/doc/RM/RM_H5S.html
- Documentation of H5Dread():
    https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
- An H5Dread() example:
    https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
*/

#include "global_errmsg_buf.h"
#include "H5DSetDescriptor.h"

#include <stdlib.h>  /* for malloc, free */


/****************************************************************************
 * Memory management of H5Viewport structs
 */

int _alloc_H5Viewport(H5Viewport *vp, int ndim, int mode)
{
	vp->h5off = NULL;
	vp->off = NULL;
	if (mode != ALLOC_OFF_AND_DIM) {
		/* Allocate memory for members 'h5off' and 'h5dim'. */
		vp->h5off = _alloc_hsize_t_buf(2 * ndim, 0,
					       "H5Viewport members");
		if (vp->h5off == NULL)
			return -1;
		vp->h5dim = vp->h5off + ndim;
	}
	if (mode != ALLOC_H5OFF_AND_H5DIM) {
		/* Allocate memory for members 'off' and 'dim'. */
		vp->off = (int *) malloc(2 * ndim * sizeof(int));
		if (vp->off == NULL) {
			if (mode != ALLOC_OFF_AND_DIM)
				free(vp->h5off);
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for H5Viewport members");
			return -1;
		}
		vp->dim = vp->off + ndim;
	}
	return 0;
}

void _free_H5Viewport(H5Viewport *vp)
{
	if (vp->h5off != NULL)
		free(vp->h5off);
	if (vp->off != NULL)
		free(vp->off);
	return;
}


/****************************************************************************
 * Other helpers
 */

hid_t _create_mem_space(int ndim, const int *dim)
{
	hsize_t *h5dim;
	int along, h5along;
	hid_t mem_space_id;

	/* Allocate and set 'h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		return -1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		h5dim[h5along] = dim[along];
	mem_space_id = H5Screate_simple(ndim, h5dim, NULL);
	if (mem_space_id < 0)
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
	free(h5dim);
	return mem_space_id;
}

int _select_H5Viewport(hid_t space_id, const H5Viewport *vp)
{
	int ret;

	ret = H5Sselect_hyperslab(space_id, H5S_SELECT_SET,
				  vp->h5off, NULL, vp->h5dim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

int _add_H5Viewport_to_h5selection(hid_t space_id, const H5Viewport *vp)
{
	int ret;

	ret = H5Sselect_hyperslab(space_id, H5S_SELECT_OR,
				  vp->h5off, NULL, vp->h5dim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

int _read_h5selection(const H5DSetDescriptor *h5dset,
		hid_t mem_space_id, void *mem,
		const H5Viewport *mem_vp)
{
	int ret;

	if (mem_vp == NULL) {
		ret = H5Sselect_all(mem_space_id);
		if (ret < 0)
		    PRINT_TO_ERRMSG_BUF("H5Sselect_all() returned an error");
	} else {
		ret = _select_H5Viewport(mem_space_id, mem_vp);
	}
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

int _read_H5Viewport(const H5DSetDescriptor *h5dset,
		const H5Viewport *h5dset_vp,
		hid_t mem_space_id, void *mem,
		const H5Viewport *mem_vp)
{
	int ret;

	ret = _select_H5Viewport(h5dset->space_id, h5dset_vp);
	if (ret < 0)
		return -1;
	return _read_h5selection(h5dset, mem_space_id, mem, mem_vp);
}


/****************************************************************************
 * _init_in_offset()
 */

void _init_in_offset(int ndim, SEXP index,
		const hsize_t *h5chunkdim, const H5Viewport *mem_vp,
		const H5Viewport *h5dset_vp,
		size_t *in_offset)
{
	size_t in_off;
	int along, h5along, i;
	SEXP start;

	in_off = 0;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunkdim[h5along];
		i = mem_vp->off[along];
		start = GET_LIST_ELT(index, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  h5dset_vp->h5off[h5along];
	}
	*in_offset = in_off;
	return;
}

