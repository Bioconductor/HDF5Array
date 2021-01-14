#ifndef _H5MREAD_SPARSE_H_
#define _H5MREAD_SPARSE_H_

#include "H5DSetDescriptor.h"
#include <Rdefines.h>

SEXP _h5mread_sparse(
	const H5DSetDescriptor *h5dset,
	SEXP index,
	int *ans_dim
);

#endif  /* _H5MREAD_SPARSE_H_ */

