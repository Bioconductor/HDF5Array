#ifndef _UASELECTION_H_
#define _UASELECTION_H_

#include <Rdefines.h>
#include "S4Vectors_interface.h"

/* Terminology:
   - uaselection: user-supplied array selection.
   - chips: the "chips" in the uaselection are its connected components i.e.
            its contiguous block-like components.
 */

/* Like VECTOR_ELT(x, i) except that 'x' can be R_NilValue. */
#define GET_LIST_ELT(x, i) ((x) != R_NilValue ? VECTOR_ELT(x, i) : R_NilValue)

static inline long long int _get_trusted_elt(SEXP x, int i)
{
	return IS_INTEGER(x) ? (long long int) INTEGER(x)[i] :
			       (long long int) REAL(x)[i];
}

int _shallow_check_uaselection(
	int ndim,
	SEXP starts,
	SEXP counts
);

long long int _check_uaselection(
	int ndim,
	const long long int *dim,
	SEXP starts,
	SEXP counts,
	int *uaselection_dim_buf
);

SEXP C_check_uaselection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

long long int _check_ordered_uaselection(
	int ndim,
	const long long int *dim,
	SEXP starts,
	SEXP counts,
	int *uaselection_dim_buf,
	int *nstart_buf,
	int *nchip_buf,
	long long int *last_chip_start_buf
);

SEXP C_check_ordered_uaselection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _uaselection_can_be_reduced(
	int ndim,
	const int *nstart,
	const int *nchip
);

SEXP _reduce_uaselection(
	int ndim,
	SEXP starts, SEXP counts,
	const int *uaselection_dim,
	const int *nchip,
	const long long int *last_chip_start
);

SEXP C_reduce_uaselection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _map_starts_to_chunks(
	int ndim,
	const long long int *dim,
	const long long int *chunkdim,
	SEXP starts,
	int *nstart_buf,
	IntAEAE *breakpoint_bufs,
	LLongAEAE *tchunkidx_bufs
);

SEXP C_map_starts_to_chunks(
	SEXP starts,
	SEXP dim,
	SEXP chunkdim
);

#endif  /* _UASELECTION_H_ */

