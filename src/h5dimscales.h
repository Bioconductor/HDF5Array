#ifndef _H5DIMSCALES_H_
#define _H5DIMSCALES_H_

#include <Rdefines.h>

SEXP C_h5isdimscale(
	SEXP filepath,
	SEXP name
);

SEXP C_h5getdimscales(
	SEXP filepath,
	SEXP name,
	SEXP scalename
);

SEXP C_h5setdimscales(
	SEXP filepath,
	SEXP name,
	SEXP dimscales,
	SEXP scalename,
	SEXP dry_run
);

SEXP C_h5getdimlabels(
	SEXP filepath,
	SEXP name
);

SEXP C_h5setdimlabels(
	SEXP filepath,
	SEXP name,
	SEXP dimlabels
);

#endif  /* _H5DIMSCALES_H_ */

