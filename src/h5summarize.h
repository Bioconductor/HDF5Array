#ifndef _H5SUMMARIZE_H_
#define _H5SUMMARIZE_H_

#include <Rdefines.h>

SEXP C_h5summarize(
	SEXP filepath,
	SEXP name,
	SEXP index,
	SEXP as_integer,
	SEXP op,
	SEXP na_rm
);

#endif  /* _H5SUMMARIZE_H_ */

