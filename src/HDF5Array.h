#include <Rdefines.h>


/* h5mread.c */

SEXP C_reduce_starts(SEXP starts);

SEXP C_h5mread(
	SEXP filepath,
	SEXP name,
	SEXP starts,
	SEXP counts
);

