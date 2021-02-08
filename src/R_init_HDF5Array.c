#include <R_ext/Rdynload.h>

#include "lzf/lzf_filter.h"  /* for register_lzf() */

#include "uaselection.h"
#include "H5DSetDescriptor.h"
#include "h5mread.h"
#include "h5dimscales.h"
#include "h5summarize.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* uaselection.c */
	CALLMETHOD_DEF(C_check_uaselection, 3),
	CALLMETHOD_DEF(C_check_ordered_uaselection, 3),
	CALLMETHOD_DEF(C_reduce_uaselection, 3),
	CALLMETHOD_DEF(C_map_starts_to_chunks, 3),

/* H5DSetDescriptor.c */
	CALLMETHOD_DEF(C_destroy_H5DSetDescriptor_xp, 1),
	CALLMETHOD_DEF(C_new_H5DSetDescriptor_xp, 3),
	CALLMETHOD_DEF(C_show_H5DSetDescriptor_xp, 1),
	CALLMETHOD_DEF(C_get_h5mread_returned_type, 3),

/* h5mread.c */
	CALLMETHOD_DEF(C_h5mread, 9),

/* h5dimscales.c */
	CALLMETHOD_DEF(C_h5isdimscale, 2),
	CALLMETHOD_DEF(C_h5getdimscales, 3),
	CALLMETHOD_DEF(C_h5setdimscales, 5),
	CALLMETHOD_DEF(C_h5getdimlabels, 2),
	CALLMETHOD_DEF(C_h5setdimlabels, 3),

/* h5summarize.c */
	CALLMETHOD_DEF(C_h5summarize, 7),

	{NULL, NULL, 0}
};

void R_init_HDF5Array(DllInfo *info)
{
	int ret;

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	/* Register the LZF filter with the hdf5 library. */
	ret = register_lzf();
	if (ret < 0)
		error("failed to register the LZF filter with the hdf5 library");

	return;
}

