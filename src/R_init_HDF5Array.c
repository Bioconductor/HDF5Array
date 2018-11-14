#include "HDF5Array.h"
#include <R_ext/Rdynload.h>

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* simple_h5read.c */
	CALLMETHOD_DEF(C_simple_h5read, 4),

	{NULL, NULL, 0}
};


void R_init_HDF5Array(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);

	return;
}

