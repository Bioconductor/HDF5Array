#ifndef _H5MREAD_H_
#define _H5MREAD_H_

#include <Rdefines.h>

SEXP C_h5mread(
	SEXP filepath,
	SEXP name,
	SEXP starts,
	SEXP counts,
	SEXP noreduce,
	SEXP as_integer,
	SEXP as_sparse,
	SEXP method
);

#endif  /* _H5MREAD_H_ */

