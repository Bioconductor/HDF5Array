#ifndef _H5MREAD_STARTSCOUNTS_H_
#define _H5MREAD_STARTSCOUNTS_H_

#include "H5DSetDescriptor.h"
#include <Rdefines.h>

SEXP _h5mread_startscounts(
	const H5DSetDescriptor *h5dset,
	SEXP starts, SEXP counts, int noreduce,
	int method,
	int *ans_dim
);

#endif  /* _H5MREAD_STARTSCOUNTS_H_ */

