/****************************************************************************
 *                    Manipulation of an array selection                    *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"
#include "S4Vectors_interface.h"

#include <limits.h>  /* for INT_MAX, LLONG_MAX, LLONG_MIN */

//#include <time.h>


char _HDF5Array_errmsg_buf[ERRMSG_BUF_LENGTH];


/****************************************************************************
 * Low-level helpers
 */

#define	NOT_A_FINITE_NUMBER(x) \
	(R_IsNA(x) || R_IsNaN(x) || (x) == R_PosInf || (x) == R_NegInf)

static inline int get_untrusted_elt(SEXP x, int i, long long int *val,
				    const char *what, int along)
{
	int tmp1;
	double tmp2;

	if (IS_INTEGER(x)) {
		tmp1 = INTEGER(x)[i];
		if (tmp1 == NA_INTEGER) {
		    if (along < 0)
			PRINT_TO_ERRMSG_BUF("%s[%d] is NA", what, i + 1);
		    else
			PRINT_TO_ERRMSG_BUF("%s[[%d]][%d] is NA",
					    what, along + 1, i + 1);
		    return -1;
		}
		*val = (long long int) tmp1;
	} else {
		tmp2 = REAL(x)[i];
		if (NOT_A_FINITE_NUMBER(tmp2)) {
		    if (along < 0)
			PRINT_TO_ERRMSG_BUF("%s[%d] is NA or NaN "
					    "or not a finite number",
					    what, i + 1);
		    else
			PRINT_TO_ERRMSG_BUF("%s[[%d]][%d] is NA or NaN "
					    "or not a finite number",
					    what, along + 1, i + 1);
		    return -1;
		}
		if (tmp2 > (double) LLONG_MAX || tmp2 < (double) LLONG_MIN) {
		    if (along < 0)
			PRINT_TO_ERRMSG_BUF("%s[%d] is too large (= %e)",
					    what, i + 1, tmp2);
		    else
			PRINT_TO_ERRMSG_BUF("%s[[%d]][%d] is too large (= %e)",
					    what, along + 1, i + 1, tmp2);
		    return -1;
		}
		*val = (long long int) tmp2;
	}
	return 0;
}

static inline void set_trusted_elt(SEXP x, int i, long long int val)
{
	if (IS_INTEGER(x))
		INTEGER(x)[i] = val;
	else
		REAL(x)[i] = val;
	return;
}


/****************************************************************************
 * Checking the array selection
 */

int _shallow_check_selection(SEXP starts, SEXP counts)
{
	int ndim;

	if (!isVectorList(starts)) {  // IS_LIST() is broken
		PRINT_TO_ERRMSG_BUF("'starts' must be a list");
		return -1;
	}
	ndim = LENGTH(starts);
	if (counts != R_NilValue) {
		if (!isVectorList(counts)) {  // IS_LIST() is broken
			PRINT_TO_ERRMSG_BUF("'counts' must "
					    "be a list (or NULL)");
			return -1;
		}
		if (LENGTH(counts) != ndim) {
			PRINT_TO_ERRMSG_BUF("'starts' and 'counts' must "
					    "have the same length");
			return -1;
		}
	}
	return ndim;
}

static void set_error_for_selection_too_large(int along1)
{
	PRINT_TO_ERRMSG_BUF("too many elements (>= 2^31) selected "
			    "along dimension %d of dataset", along1);
	return;
}

static void set_errmsg_for_non_strictly_ascending_selection(int along1, int i,
							    int starts_only)
{
	const char *msg = "selection must be strictly ascending "
			  "along each dimension, but\n  you have:";
	if (starts_only)
		PRINT_TO_ERRMSG_BUF("%s starts[[%d]][%d] <= starts[[%d]][%d]",
				    msg, along1, i + 1, along1, i);
	else
		PRINT_TO_ERRMSG_BUF("%s starts[[%d]][%d] < starts[[%d]][%d] + "
				    "counts[[%d]][%d]",
				    msg, along1, i + 1, along1, i, along1, i);
	return;
}

static void set_errmsg_for_selection_beyond_dim(int along1, int i,
						int starts_only)
{
	const char *msg = "selection must be within extent of "
			  "dataset, but you\n  have:";
	if (starts_only)
		PRINT_TO_ERRMSG_BUF(
			"%s starts[[%d]][%d] "
			"> dimension %d in dataset",
			msg, along1, i + 1, along1);
	else
		PRINT_TO_ERRMSG_BUF(
			"%s starts[[%d]][%d] + counts[[%d]][%d] - 1 "
			"> dimension %d in dataset",
			msg, along1, i + 1, along1, i + 1, along1);
	return;
}

static inline int get_untrusted_start(SEXP start, int i, long long int *s,
				      int along,
				      long long int e,
				      int *nblock,
				      long long int *last_block_start,
				      int starts_only)
{
	int ret;

	ret = get_untrusted_elt(start, i, s, "starts", along);
	if (ret < 0)
		return -1;
	if (*s <= e) {
		if (e == 0) {
			PRINT_TO_ERRMSG_BUF("starts[[%d]][%d] is <= 0",
					    along + 1, i + 1);
			return -1;
		}
		set_errmsg_for_non_strictly_ascending_selection(
				along + 1, i, starts_only);
		return -1;
	}
	if (i == 0 || *s != e + 1) {
		nblock[along]++;
		last_block_start[along] = *s;
	}
	return 0;
}

static int check_selection_along(int along, SEXP starts, SEXP counts,
				 long long int d,
				 int *nstart, int *count_sum,
				 int *nblock, long long int *last_block_start)
{
	SEXP start, count;
	int n, i, ret;
	long long int e, s, cs, c;

	start = VECTOR_ELT(starts, along);
	count = counts != R_NilValue ? VECTOR_ELT(counts, along) : R_NilValue;
	if (start == R_NilValue) {
		if (count != R_NilValue) {
			PRINT_TO_ERRMSG_BUF(
				"if 'starts[[%d]]' is NULL then 'counts' "
				"or 'counts[[%d]]' must also be NULL",
				along + 1, along + 1);
			return -1;
		}
		if (d >= 0) {
			if (d > INT_MAX) {
				set_error_for_selection_too_large(along + 1);
				return -1;
			}
			nstart[along] = count_sum[along] = d;
			nblock[along] = d != 0;
			last_block_start[along] = 1;
		} else {
			/* 'count_sum' is undefined in that case. */
			nstart[along] = nblock[along] = 1;
			last_block_start[along] = 1;
		}
		return 0;
	}
	if (!(IS_INTEGER(start) || IS_NUMERIC(start))) {
		PRINT_TO_ERRMSG_BUF(
			"'starts[[%d]]' must be an "
			"integer vector (or NULL)", along + 1);
		return -1;
	}
	n = LENGTH(start);
	if (count != R_NilValue) {
		if (!(IS_INTEGER(count) || IS_NUMERIC(count))) {
			PRINT_TO_ERRMSG_BUF(
				"'counts[[%d]]' must be an "
				"integer vector (or NULL)", along + 1);
			return -1;
		}
		if (LENGTH(count) != n) {
			PRINT_TO_ERRMSG_BUF(
				"'starts[[%d]]' and 'counts[[%d]]' "
				"must have the same length",
				along + 1, along + 1);
			return -1;
		}
	}
	nstart[along] = n;
	nblock[along] = 0;
	e = 0;
	if (count == R_NilValue) {
		for (i = 0; i < n; i++) {
			ret = get_untrusted_start(start, i, &s, along,
						  e, nblock, last_block_start,
						  1);
			if (ret < 0)
				return -1;
			e = s;
			if (d >= 0 && e > d) {
				set_errmsg_for_selection_beyond_dim(
					along + 1, i, 1);
				return -1;
			}
		}
		count_sum[along] = n;
	} else {
		cs = 0;
		for (i = 0; i < n; i++) {
			ret = get_untrusted_start(start, i, &s, along,
						  e, nblock, last_block_start,
						  0);
			if (ret < 0)
				return -1;
			ret = get_untrusted_elt(count, i, &c, "counts", along);
			if (ret < 0)
				return -1;
			if (c <= 0) {
				PRINT_TO_ERRMSG_BUF("counts[[%d]][%d] is <= 0",
						    along + 1, i + 1);
				return -1;
			}
			e = s + c - 1;	// could overflow! (FIXME)
			if (d >= 0 && e > d) {
				set_errmsg_for_selection_beyond_dim(
					along + 1, i, 0);
				return -1;
			}
			cs += c;	// could overflow! (FIXME)
			if (cs > INT_MAX) {
				set_error_for_selection_too_large(along + 1);
				return -1;
			}
		}
		count_sum[along] = cs;
	}
	return 0;
}

int _deep_check_selection(SEXP starts, SEXP counts,
			  const long long int *dim,
			  int *nstart, int *count_sum,
			  int *nblock, long long int *last_block_start)
{
	int ndim, along, ret;

	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		ret = check_selection_along(along, starts, counts,
					    dim != NULL ? dim[along] : -1,
					    nstart, count_sum,
					    nblock, last_block_start);
		if (ret < 0)
			return -1;
	}
	return 0;
}


/****************************************************************************
 * Reduce the array selection
 */

int _selection_can_be_reduced(int ndim, const int *nstart, const int *nblock)
{
	int along;

	for (along = 0; along < ndim; along++) {
		/* nblock[along] should always be <= nstart[along] */
		if (nblock[along] < nstart[along])
			return 1;
	}
	return 0;
}

static SEXP dup_or_coerce_to_INTSXP(SEXP x, int dup)
{
	int x_len, i;
	SEXP ans;

	if (dup)
		return duplicate(x);
	x_len = LENGTH(x);
	ans = PROTECT(NEW_INTEGER(x_len));
	for (i = 0; i < x_len; i++)
		INTEGER(ans)[i] = (int) REAL(x)[i];
	UNPROTECT(1);
	return ans;
}

/*
 * Note that this does something similar to what coercion from integer (or
 * numeric) to IRanges does (see .Call entry point "IRanges_from_integer"
 * in IRanges). However we cannot re-use this here because we want to be able
 * to handle start values that are >= 2^31 which this coercion doesn't support
 * at the moment.
 */
static void stitch_selection(SEXP start_in, SEXP count_in,
			     SEXP start_out, int *count_out)
{
	int n, i, j;
	long long int e, s, c;

	n = LENGTH(start_in);
	e = 0;
	j = -1;
	if (count_in == R_NilValue) {
		for (i = 0; i < n; i++) {
			s = _get_trusted_elt(start_in, i);
			if (i == 0 || s != e + 1) {
				j++;
				set_trusted_elt(start_out, j, s);
				count_out[j] = 1;
			} else {
				count_out[j]++;
			}
			e = s;
		}
	} else {
		for (i = 0; i < n; i++) {
			s = _get_trusted_elt(start_in, i);
			c = _get_trusted_elt(count_in, i);
			if (i == 0 || s != e + 1) {
				j++;
				set_trusted_elt(start_out, j, s);
				count_out[j] = c;
			} else {
				count_out[j] += c;
			}
			e = s + c - 1;
		}
	}
	return;
}

static void reduce_selection_along(int along, SEXP starts, SEXP counts,
				   const int *count_sum,
				   const int *nblock,
				   const long long int *last_block_start,
				   SEXP reduced_starts, SEXP reduced_counts)
{
	SEXP start, count, reduced_start, reduced_count;
	int n, dup;
	SEXPTYPE type;

	start = VECTOR_ELT(starts, along);
	if (start == R_NilValue)
		return;
	n = LENGTH(start);
	count = counts != R_NilValue ? VECTOR_ELT(counts, along) : R_NilValue;
	if (nblock[along] == n) {
		/* Nothing to stitch. */
		dup = IS_INTEGER(start) || last_block_start[along] > INT_MAX;
		reduced_start = PROTECT(dup_or_coerce_to_INTSXP(start, dup));
		SET_VECTOR_ELT(reduced_starts, along, reduced_start);
		UNPROTECT(1);
		if (count_sum[along] == n)
			return;
		dup = IS_INTEGER(count);
		reduced_count = PROTECT(dup_or_coerce_to_INTSXP(count, dup));
		SET_VECTOR_ELT(reduced_counts, along, reduced_count);
		UNPROTECT(1);
		return;
	}
	/* Stitch. */
	type = last_block_start[along] <= INT_MAX ? INTSXP : REALSXP;
	reduced_start = PROTECT(allocVector(type, nblock[along]));
	SET_VECTOR_ELT(reduced_starts, along, reduced_start);
	UNPROTECT(1);
	reduced_count = PROTECT(NEW_INTEGER(nblock[along]));
	SET_VECTOR_ELT(reduced_counts, along, reduced_count);
	UNPROTECT(1);
	stitch_selection(start, count, reduced_start, INTEGER(reduced_count));
	return;
}

SEXP _reduce_selection(SEXP starts, SEXP counts,
		       const int *count_sum,
		       const int *nblock,
		       const long long int *last_block_start)
{
	int ndim, along;
	SEXP ans, reduced_starts, reduced_counts;

	//clock_t t0 = clock();
	ndim = LENGTH(starts);
	ans = PROTECT(NEW_LIST(2));
	reduced_starts = PROTECT(NEW_LIST(ndim));
	SET_VECTOR_ELT(ans, 0, reduced_starts);
	UNPROTECT(1);
	reduced_counts = PROTECT(NEW_LIST(ndim));
	SET_VECTOR_ELT(ans, 1, reduced_counts);
	UNPROTECT(1);
	for (along = 0; along < ndim; along++) {
		reduce_selection_along(along, starts, counts,
				       count_sum,
				       nblock, last_block_start,
				       reduced_starts, reduced_counts);
	}
	UNPROTECT(1);
	//printf("time 2nd pass: %e\n", (1.0 * clock() - t0) / CLOCKS_PER_SEC);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Negative values in 'dim' are treated as infinite dimensions.
 * Return a list of length 2 or NULL if the selection could not be reduced.
 * The 1st list element is the list of reduced starts and the 2nd list element
 * the list of reduced 'counts'. The 2 lists have the same shape i.e. same
 * length() and same lengths().
 */
SEXP C_reduce_selection(SEXP starts, SEXP counts, SEXP dim)
{
	int ndim, along, ret;
	LLongAE *dim_buf;
	IntAE *nstart_buf, *count_sum_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	int *nstart, *count_sum, *nblock;
	long long int *dim_p, d, *last_block_start;

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0)
		error(_HDF5Array_errmsg_buf);
	if (dim == R_NilValue) {
		dim_p = NULL;
	} else {
		if (!(IS_INTEGER(dim) || IS_NUMERIC(dim)))
			error("'dim' must be an integer vector (or NULL)");
		if (LENGTH(dim) != ndim)
			error("'starts' and 'dim' must have the same length");
		dim_buf = new_LLongAE(ndim, ndim, 0);
		dim_p = dim_buf->elts;
		for (along = 0; along < ndim; along++) {
			ret = get_untrusted_elt(dim, along, &d, "dim", -1);
			if (ret < 0)
				error(_HDF5Array_errmsg_buf);
			dim_p[along] = d;
		}
	}

	nstart_buf = new_IntAE(ndim, ndim, 0);
	count_sum_buf = new_IntAE(ndim, ndim, 0);
	nblock_buf = new_IntAE(ndim, ndim, 0);
	last_block_start_buf = new_LLongAE(ndim, ndim, 0);

	nstart = nstart_buf->elts;
	count_sum = count_sum_buf->elts;
	nblock = nblock_buf->elts;
	last_block_start = last_block_start_buf->elts;

	/* 1st pass */
	//clock_t t0 = clock();
	ret = _deep_check_selection(starts, counts, dim_p,
				    nstart, count_sum,
				    nblock, last_block_start);
	//printf("time 1st pass: %e\n", (1.0 * clock() - t0) / CLOCKS_PER_SEC);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf);

	if (!_selection_can_be_reduced(ndim, nstart, nblock))
		return R_NilValue;

	/* 2nd pass */
	return _reduce_selection(starts, counts,
				 count_sum, nblock, last_block_start);
}

