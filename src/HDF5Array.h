#include <Rdefines.h>


/* h5mread.c */

SEXP C_reduce_selection(SEXP starts, SEXP counts, SEXP dim);

SEXP C_h5mread(
	SEXP filepath,
	SEXP name,
	SEXP starts,
	SEXP counts,
	SEXP noreduce,
	SEXP method,
	SEXP as_integer
);

