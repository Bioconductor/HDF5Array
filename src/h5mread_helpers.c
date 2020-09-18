/****************************************************************************
 *  Some low-level helper functions to support the various h5mread methods  *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_helpers.h"

#include "global_errmsg_buf.h"
#include "H5DSetDescriptor.h"

#include <stdlib.h>  /* for malloc, free */


/****************************************************************************
 * Basic manipulation of a H5Viewport struct
 */

/* 'mode' controls what fields we should allocate memory for:
 *     mode == 0: 'h5off' and 'h5dim' fields only
 *     mode == 1: 'off' and 'dim' fields only
 *     mode == 2: for all the fields
 */
int _alloc_H5Viewport(H5Viewport *vp, int ndim, int mode)
{
	vp->h5off = NULL;
	vp->off = NULL;
	if (mode != 1) {
		/* Allocate memory for the 'h5off' and 'h5dim' fields. */
		vp->h5off = _alloc_hsize_t_buf(2 * ndim, 0,
					       "H5Viewport fields");
		if (vp->h5off == NULL)
			return -1;
		vp->h5dim = vp->h5off + ndim;
	}
	if (mode != 0) {
		/* Allocate memory for the 'off' and 'dim' fields. */
		vp->off = (int *) malloc(2 * ndim * sizeof(int));
		if (vp->off == NULL) {
			if (mode != 1)
				free(vp->h5off);
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for H5Viewport fields");
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

int _add_H5Viewport_to_selection(hid_t space_id, const H5Viewport *vp)
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


/****************************************************************************
 * Misc helpers
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

