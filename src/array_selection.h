#ifndef _ARRAY_SELECTION_H_
#define _ARRAY_SELECTION_H_

#include <Rdefines.h>
#include "S4Vectors_interface.h"

/* Like VECTOR_ELT(x, i) except that 'x' can be R_NilValue. */
#define GET_LIST_ELT(x, i) ((x) != R_NilValue ? VECTOR_ELT(x, i) : R_NilValue)

static inline long long int _get_trusted_elt(SEXP x, int i)
{
	return IS_INTEGER(x) ? (long long int) INTEGER(x)[i] :
			       (long long int) REAL(x)[i];
}

int _shallow_check_selection(
	int ndim,
	SEXP starts,
	SEXP counts
);

long long int _check_selection(
	int ndim,
	const long long int *dim,
	SEXP starts,
	SEXP counts,
	int *selection_dim_buf
);

SEXP C_check_selection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

long long int _check_ordered_selection(
	int ndim,
	const long long int *dim,
	SEXP starts,
	SEXP counts,
	int *selection_dim_buf,
	int *nstart_buf,
	int *nblock_buf,
	long long int *last_block_start_buf
);

SEXP C_check_ordered_selection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _selection_can_be_reduced(
	int ndim,
	const int *nstart,
	const int *nblock
);

SEXP _reduce_selection(
	int ndim,
	SEXP starts, SEXP counts,
	const int *selection_dim,
	const int *nblock,
	const long long int *last_block_start
);

SEXP C_reduce_selection(
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

#endif  /* _ARRAY_SELECTION_H_ */

