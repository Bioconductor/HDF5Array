#include <R_ext/Rdynload.h>

#include "array_selection.h"
#include "H5DSetDescriptor.h"
#include "h5mread.h"
#include "h5dimscales.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* array_selection.c */
	CALLMETHOD_DEF(C_check_selection, 3),
	CALLMETHOD_DEF(C_check_ordered_selection, 3),
	CALLMETHOD_DEF(C_reduce_selection, 3),
	CALLMETHOD_DEF(C_map_starts_to_chunks, 3),

/* H5DSetDescriptor.c */
	CALLMETHOD_DEF(C_destroy_H5DSetDescriptor_xp, 1),
	CALLMETHOD_DEF(C_new_H5DSetDescriptor_xp, 3),
	CALLMETHOD_DEF(C_show_H5DSetDescriptor_xp, 1),
	CALLMETHOD_DEF(C_get_h5mread_returned_type, 3),

/* h5mread.c */
	CALLMETHOD_DEF(C_h5mread, 8),

/* h5dimscales.c */
	CALLMETHOD_DEF(C_h5isdimscale, 2),
	CALLMETHOD_DEF(C_h5getdimscales, 3),
	CALLMETHOD_DEF(C_h5setdimscales, 5),
	CALLMETHOD_DEF(C_h5getdimlabels, 2),
	CALLMETHOD_DEF(C_h5setdimlabels, 3),

	{NULL, NULL, 0}
};

void R_init_HDF5Array(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	return;
}

