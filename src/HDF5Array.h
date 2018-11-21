#include <Rdefines.h>


#define ERRMSG_BUF_LENGTH 256

#define PRINT_TO_ERRMSG_BUF(...) \
	snprintf(_HDF5Array_errmsg_buf, ERRMSG_BUF_LENGTH, __VA_ARGS__)

inline long long int _get_trusted_elt(SEXP x, int i)
{
	return IS_INTEGER(x) ? (long long int) INTEGER(x)[i] :
			       (long long int) REAL(x)[i];
}


/* array_selection.c */

char _HDF5Array_errmsg_buf[ERRMSG_BUF_LENGTH];

int _shallow_check_selection(
	SEXP starts,
	SEXP counts
);

int _deep_check_selection(
	SEXP starts,
	SEXP counts,
	const long long int *dim,
	int *nstart,
	int *count_sum,
	int *nblock,
	long long int *last_block_start
);

int _selection_can_be_reduced(
	int ndim,
	const int *nstart,
	const int *nblock
);

SEXP _reduce_selection(
	SEXP starts, SEXP counts,
	const int *count_sum,
	const int *nblock,
	const long long int *last_block_start
);

SEXP C_reduce_selection(
	SEXP starts,
	SEXP counts,
	SEXP dim
);


/* h5mread.c */

SEXP C_h5mread(
	SEXP filepath,
	SEXP name,
	SEXP starts,
	SEXP counts,
	SEXP noreduce,
	SEXP method,
	SEXP as_integer
);

