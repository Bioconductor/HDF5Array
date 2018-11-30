/****************************************************************************
 *                    Manipulation of an array selection                    *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include <limits.h>  /* for INT_MAX, LLONG_MAX, LLONG_MIN */

//#include <time.h>


char _HDF5Array_errmsg_buf[ERRMSG_BUF_LENGTH];


/****************************************************************************
 * Low-level helpers
 */

static int check_INTEGER_or_NUMERIC(SEXP x, const char *what, int along)
{
	if (!(IS_INTEGER(x) || IS_NUMERIC(x))) {
		PRINT_TO_ERRMSG_BUF("'%s[[%d]]' must be an "
				    "integer vector (or NULL)",
				    what, along + 1);
		return -1;
	}
	return 0;
}

static int shallow_check_count(SEXP count, int n, int along)
{
	if (count == R_NilValue)
		return 0;
	if (check_INTEGER_or_NUMERIC(count, "counts", along) < 0)
		return -1;
	if (LENGTH(count) != n) {
		PRINT_TO_ERRMSG_BUF("'starts[[%d]]' and 'counts[[%d]]' "
				    "must have the same length",
				    along + 1, along + 1);
		return -1;
	}
	return 0;
}

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

/* Called at the very beginning of the various .Call entry points where
   it's used and before any resource is allocated so it's ok to error()
   immediately in case of error. */
static const long long int *check_dim(SEXP dim, int ndim)
{
	long long int *dim_p, d;
	int along, ret;

	if (dim == R_NilValue)
		return NULL;
	if (!(IS_INTEGER(dim) || IS_NUMERIC(dim)))
		error("'dim' must be an integer vector (or NULL)");
	if (LENGTH(dim) != ndim)
		error("'starts' and 'dim' must have the same length");
	dim_p = new_LLongAE(ndim, ndim, 0)->elts;
	for (along = 0; along < ndim; along++) {
		ret = get_untrusted_elt(dim, along, &d, "dim", -1);
		if (ret < 0)
			error(_HDF5Array_errmsg_buf);
		dim_p[along] = d;
	}
	return dim_p;
}


/****************************************************************************
 * Shallow check of an array selection
 */

/* Only check that 'starts' is a list and that 'counts' is NULL or a list
   of the same length as 'starts'.
   Return the nb of list elements in 'starts'. */
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


/****************************************************************************
 * Deep check of an array selection
 */

static void set_error_for_selection_too_large(int along1)
{
	PRINT_TO_ERRMSG_BUF("too many elements (>= 2^31) selected "
			    "along dimension %d of array", along1);
	return;
}

static void set_errmsg_for_selection_beyond_dim(int along1, int i,
						int no_counts)
{
	const char *msg = "selection must be within extent of "
			  "array, but you\n  have:";
	if (no_counts)
		PRINT_TO_ERRMSG_BUF(
			"%s starts[[%d]][%d] "
			"> dimension %d in array",
			msg, along1, i + 1, along1);
	else
		PRINT_TO_ERRMSG_BUF(
			"%s starts[[%d]][%d] + counts[[%d]][%d] - 1 "
			"> dimension %d in array",
			msg, along1, i + 1, along1, i + 1, along1);
	return;
}

static void set_errmsg_for_non_strictly_ascending_selection(int along1, int i,
							    int no_counts)
{
	const char *msg = "selection must be strictly ascending "
			  "along each dimension, but\n  you have:";
	if (no_counts)
		PRINT_TO_ERRMSG_BUF("%s starts[[%d]][%d] <= starts[[%d]][%d]",
				    msg, along1, i + 1, along1, i);
	else
		PRINT_TO_ERRMSG_BUF("%s starts[[%d]][%d] < starts[[%d]][%d] + "
				    "counts[[%d]][%d]",
				    msg, along1, i + 1, along1, i, along1, i);
	return;
}

static inline int get_untrusted_start(SEXP start, int i, long long int *s,
				      long long int min_start,
				      int along, int no_counts)
{
	if (get_untrusted_elt(start, i, s, "starts", along) < 0)
		return -1;
	if (*s <= 0) {
		PRINT_TO_ERRMSG_BUF("starts[[%d]][%d] is <= 0",
				    along + 1, i + 1);
		return -1;
	}
	if (*s < min_start) {
		set_errmsg_for_non_strictly_ascending_selection(
			along + 1, i, no_counts);
		return -1;
	}
	return 0;
}

static int check_selection_along(SEXP start, SEXP count, int along,
				 long long int d)
{
	long long int selection_dim, s, c, e;
	int n, i, ret;

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
			selection_dim = d;
		} else {
			/* The dimension of the selection along the current
			   dimension is undefined in that case.
			   We **arbitrary** set it to INT_MAX. */
			selection_dim = INT_MAX;
		}
		return (int) selection_dim;
	}
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	n = LENGTH(start);
	if (shallow_check_count(count, n, along) < 0)
		return -1;
	/* Walk on the 'start' elements. */
	for (i = 0; i < n; i++) {
		/* Because we've set the 4th argument ('min_start') to 0,
		   the last argument ('no_counts') will be ignored. */
		ret = get_untrusted_start(start, i, &s, 0, along, 0);
		if (ret < 0)
			return -1;
		if (d >= 0 && s > d) {
			set_errmsg_for_selection_beyond_dim(
				along + 1, i, 1);
			return -1;
		}
	}
	if (count == R_NilValue)
		return n;
	/* Walk on the 'count' (and 'start') elements. */
	selection_dim = 0;
	for (i = 0; i < n; i++) {
		ret = get_untrusted_elt(count, i, &c, "counts", along);
		if (ret < 0)
			return -1;
		if (c == 0)
			continue;
		if (c < 0) {
			PRINT_TO_ERRMSG_BUF("counts[[%d]][%d] is < 0",
					    along + 1, i + 1);
			return -1;
		}
		s = _get_trusted_elt(start, i);
		e = s + c - 1;      // could overflow! (FIXME)
		if (d >= 0 && e > d) {
			set_errmsg_for_selection_beyond_dim(
				along + 1, i, 0);
			return -1;
		}
		selection_dim += c; // could overflow! (FIXME)
		if (selection_dim > INT_MAX) {
			set_error_for_selection_too_large(along + 1);
			return -1;
		}
	}
	return (int) selection_dim;
}

/* Assume that 'starts' is a list and that 'counts' is NULL or a list
   of the same length as 'starts'. This should have been checked
   by _shallow_check_selection() already so is not checked again.

   'dim' is assumed to be NULL or to have the same length as 'starts'.

   'selection_dim_buf' is assumed to have the same length as 'starts'.
*/
long long int _check_selection(SEXP starts, SEXP counts,
		     const long long int *dim, int *selection_dim_buf)
{
	long long int selection_len;
	int ndim, along, selection_dim;
	SEXP start, count;

	selection_len = 1;
	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		count = counts != R_NilValue ? VECTOR_ELT(counts, along)
					     : R_NilValue;
		selection_dim = check_selection_along(
					start, count, along,
					dim != NULL ? dim[along] : -1);
		if (selection_dim < 0)
			return -1;
		selection_dim_buf[along] = selection_dim;
		selection_len *= selection_dim;
	}
	return selection_len;
}

/* --- .Call ENTRY POINT ---
 * Return the dimensions of the selection.
 */
SEXP C_check_selection(SEXP starts, SEXP counts, SEXP dim)
{
	int ndim;
	const long long int *dim_p;
	IntAE *selection_dim_buf;
	long long int selection_len;

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0)
		error(_HDF5Array_errmsg_buf);
	dim_p = check_dim(dim, ndim);
	selection_dim_buf = new_IntAE(ndim, ndim, 0);
	selection_len = _check_selection(starts, counts, dim_p,
					 selection_dim_buf->elts);
	if (selection_len < 0)
		error(_HDF5Array_errmsg_buf);
	return new_INTEGER_from_IntAE(selection_dim_buf);
}


/****************************************************************************
 * Deep check of an ordered array selection (in preparation for reduction)
 */

static int check_ordered_selection_along_NULL_start(SEXP count, int along,
			long long int d,
			int *nstart_buf, int *nblock_buf,
			long long int *last_block_start_buf)
{
	int selection_dim;

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
		selection_dim = d;
		nstart_buf[along] = d;
		nblock_buf[along] = d != 0;
		last_block_start_buf[along] = 1;
	} else {
		/* The dimension of the selection along the current
		   dimension is undefined in that case.
		   We **arbitrary** set it to INT_MAX. */
		selection_dim = INT_MAX;
		nstart_buf[along] = nblock_buf[along] = 1;
		last_block_start_buf[along] = 1;
	}
	return selection_dim;
}

static int check_ordered_selection_along(SEXP start, SEXP count, int along,
			long long int d,
			int *nstart_buf, int *nblock_buf,
			long long int *last_block_start_buf)
{
	long long int selection_dim, min_start, s, c;
	int n, i, ret;

	if (start == R_NilValue)
		return check_ordered_selection_along_NULL_start(count, along,
				d,
				nstart_buf, nblock_buf, last_block_start_buf);
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	n = LENGTH(start);
	if (shallow_check_count(count, n, along) < 0)
		return -1;
	nstart_buf[along] = n;
	nblock_buf[along] = 0;
	min_start = 1;
	if (count == R_NilValue) {
		/* Walk on the 'start' elements. */
		for (i = 0; i < n; i++) {
			ret = get_untrusted_start(start, i, &s, min_start,
						  along, 1);
			if (ret < 0)
				return -1;
			if (s != min_start) {
				nblock_buf[along]++;
				last_block_start_buf[along] = s;
			}
			min_start = s + 1;
			if (d >= 0 && s > d) {
				set_errmsg_for_selection_beyond_dim(
					along + 1, i, 1);
				return -1;
			}
		}
		selection_dim = n;
	} else {
		/* Walk on the 'start' and 'count' elements. */
		selection_dim = 0;
		for (i = 0; i < n; i++) {
			ret = get_untrusted_elt(count, i, &c, "counts", along);
			if (ret < 0)
				return -1;
			if (c < 0) {
				PRINT_TO_ERRMSG_BUF("counts[[%d]][%d] is < 0",
						    along + 1, i + 1);
				return -1;
			}
			if (c == 0)
				continue;
			ret = get_untrusted_start(start, i, &s, min_start,
						  along, 0);
			if (ret < 0)
				return -1;
			if (s != min_start) {
				nblock_buf[along]++;
				last_block_start_buf[along] = s;
			}
			min_start = s + c;  // could overflow! (FIXME)
			if (d >= 0 && min_start - 1 > d) {
				set_errmsg_for_selection_beyond_dim(
					along + 1, i, 0);
				return -1;
			}
			selection_dim += c; // could overflow! (FIXME)
			if (selection_dim > INT_MAX) {
				set_error_for_selection_too_large(along + 1);
				return -1;
			}
		}
	}
	return (int) selection_dim;
}

/* Assume that 'starts' is a list and that 'counts' is NULL or a list
   of the same length as 'starts'. This should have been checked
   by _shallow_check_selection() already so is not checked again.

   'dim' is assumed to be NULL or to have the same length as 'starts'.

   'selection_dim_buf', 'nstart_buf', 'nblock_buf', and 'last_block_start_buf'
   are assumed to have the same length as 'starts'.
*/
long long int _check_ordered_selection(SEXP starts, SEXP counts,
			const long long int *dim, int *selection_dim_buf,
			int *nstart_buf, int *nblock_buf,
			long long int *last_block_start_buf)
{
	long long int selection_len;
	int ndim, along, selection_dim;
	SEXP start, count;

	selection_len = 1;
	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		count = counts != R_NilValue ? VECTOR_ELT(counts, along)
					     : R_NilValue;
		selection_dim = check_ordered_selection_along(
					start, count, along,
					dim != NULL ? dim[along] : -1,
					nstart_buf, nblock_buf,
					last_block_start_buf);
		if (selection_dim < 0)
			return -1;
		selection_dim_buf[along] = selection_dim;
		selection_len *= selection_dim;
	}
	return selection_len;
}

/* --- .Call ENTRY POINT ---
 * Return the dimensions of the selection.
 */
SEXP C_check_ordered_selection(SEXP starts, SEXP counts, SEXP dim)
{
	int ndim;
	const long long int *dim_p;
	IntAE *selection_dim_buf, *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	long long int selection_len;

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0)
		error(_HDF5Array_errmsg_buf);
	dim_p = check_dim(dim, ndim);
	selection_dim_buf = new_IntAE(ndim, ndim, 0);
	nstart_buf = new_IntAE(ndim, ndim, 0);
	nblock_buf = new_IntAE(ndim, ndim, 0);
	last_block_start_buf = new_LLongAE(ndim, ndim, 0);
	selection_len = _check_ordered_selection(starts, counts, dim_p,
				selection_dim_buf->elts,
				nstart_buf->elts, nblock_buf->elts,
				last_block_start_buf->elts);
	if (selection_len < 0)
		error(_HDF5Array_errmsg_buf);
	return new_INTEGER_from_IntAE(selection_dim_buf);
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
	long long int min_start, s, c;

	n = LENGTH(start_in);
	min_start = 0;
	j = -1;
	if (count_in == R_NilValue) {
		for (i = 0; i < n; i++) {
			s = _get_trusted_elt(start_in, i);
			if (s != min_start) {
				j++;
				set_trusted_elt(start_out, j, s);
				count_out[j] = 1;
			} else {
				count_out[j]++;
			}
			min_start = s;
		}
	} else {
		for (i = 0; i < n; i++) {
			c = _get_trusted_elt(count_in, i);
			if (c == 0)
				continue;
			s = _get_trusted_elt(start_in, i);
			if (s != min_start) {
				j++;
				set_trusted_elt(start_out, j, s);
				count_out[j] = c;
			} else {
				count_out[j] += c;
			}
			min_start = s + c;
		}
	}
	return;
}

static void reduce_selection_along(SEXP start, SEXP count, int along,
				   const int *selection_dim,
				   const int *nblock,
				   const long long int *last_block_start,
				   SEXP reduced_starts, SEXP reduced_counts)
{
	int n, dup;
	SEXP reduced_start, reduced_count;
	SEXPTYPE type;

	n = LENGTH(start);
	if (nblock[along] == n) {
		/* Nothing to stitch. */
		dup = IS_INTEGER(start) || last_block_start[along] > INT_MAX;
		reduced_start = PROTECT(dup_or_coerce_to_INTSXP(start, dup));
		SET_VECTOR_ELT(reduced_starts, along, reduced_start);
		UNPROTECT(1);
		if (selection_dim[along] == n)
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
		       const int *selection_dim,
		       const int *nblock,
		       const long long int *last_block_start)
{
	int ndim, along;
	SEXP ans, reduced_starts, reduced_counts, start, count;

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
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		count = counts != R_NilValue ? VECTOR_ELT(counts, along) :
					       R_NilValue;
		reduce_selection_along(start, count, along,
				       selection_dim,
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
 * When returning a list of length 2:
 *   - The 1st list element is the list of reduced starts.
 *   - The 2nd list element is the list of reduced counts.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
 */
SEXP C_reduce_selection(SEXP starts, SEXP counts, SEXP dim)
{
	int ndim;
	const long long int *dim_p;
	IntAE *selection_dim_buf, *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	long long int selection_len;

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0)
		error(_HDF5Array_errmsg_buf);
	dim_p = check_dim(dim, ndim);
	selection_dim_buf = new_IntAE(ndim, ndim, 0);
	nstart_buf = new_IntAE(ndim, ndim, 0);
	nblock_buf = new_IntAE(ndim, ndim, 0);
	last_block_start_buf = new_LLongAE(ndim, ndim, 0);

	/* 1st pass */
	selection_len = _check_ordered_selection(starts, counts, dim_p,
				selection_dim_buf->elts,
				nstart_buf->elts, nblock_buf->elts,
				last_block_start_buf->elts);
	if (selection_len < 0)
		error(_HDF5Array_errmsg_buf);
	if (!_selection_can_be_reduced(ndim,
				       nstart_buf->elts,
				       nblock_buf->elts))
		return R_NilValue;

	/* 2nd pass */
	return _reduce_selection(starts, counts,
				 selection_dim_buf->elts, nblock_buf->elts,
				 last_block_start_buf->elts);
}


/****************************************************************************
 * Map array selection to chunks
 */

static int map_start_to_chunks(SEXP start, int along,
		long long int d, long long int chunkd,
		int *nstart_buf,
		IntAE *breakpoint_buf, LLongAE *chunkidx_buf)
{
	int n, i, ret;
	size_t nchunk;
	long long int min_start, s, chunkidx, prev_chunkidx;

	if (start == R_NilValue) {
		if (d > INT_MAX) {
			set_error_for_selection_too_large(along + 1);
			return -1;
		}
		nstart_buf[along] = d;
		return 0;
	}

	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;

	if (IntAE_get_nelt(breakpoint_buf) != 0 ||
	    LLongAE_get_nelt(chunkidx_buf) != 0) {
		/* Should never happen! */
		PRINT_TO_ERRMSG_BUF("internal error: map_start_to_chunks() "
				    "was called with non-empty breakpoint "
				    "or chunkidx buffers");
		return -1;
	}

	n = LENGTH(start);
	nstart_buf[along] = n;

	if (n == 0)
		return 0;

	/* Get 's' and 'chunkidx' for 1st 'start' element. */
	ret = get_untrusted_start(start, 0, &s, 1, along, 1);
	if (ret < 0)
		return -1;
	if (s > d) {
		set_errmsg_for_selection_beyond_dim(along + 1, 0, 1);
		return -1;
	}
	chunkidx = (s - 1) / chunkd;

	/* Walk on the remaining 'start' elements. */
	nchunk = 0;
	for (i = 1; i < n; i++) {
		min_start = s;
		ret = get_untrusted_start(start, i, &s, min_start, along, 1);
		if (ret < 0)
			return -1;
		if (s > d) {
			set_errmsg_for_selection_beyond_dim(along + 1, i, 1);
			return -1;
		}
		prev_chunkidx = chunkidx;
		chunkidx = (s - 1) / chunkd;
		if (chunkidx > prev_chunkidx) {
			IntAE_insert_at(breakpoint_buf, nchunk, i);
			LLongAE_insert_at(chunkidx_buf, nchunk, prev_chunkidx);
			nchunk++;
		}
	}
	IntAE_insert_at(breakpoint_buf, nchunk, n);
	LLongAE_insert_at(chunkidx_buf, nchunk, chunkidx);
	return 0;
}

/* Assume that 'starts' is a list. This should have been checked by
   _shallow_check_selection() already so is not checked again.

   'dim', 'chunk_spacings', 'nstart_buf', 'breakpoint_bufs', and
   'chunkidx_bufs' are assumed to have the same length as 'starts'.
*/
int _map_starts_to_chunks(SEXP starts,
		const long long int *dim,
		const long long int *chunk_spacings,
		int *nstart_buf,
		IntAEAE *breakpoint_bufs, LLongAEAE *chunkidx_bufs)
{
	int ndim, along, ret;
	SEXP start;

	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		ret = map_start_to_chunks(start, along,
					  dim[along], chunk_spacings[along],
					  nstart_buf,
					  breakpoint_bufs->elts[along],
					  chunkidx_bufs->elts[along]);
		if (ret < 0)
			return -1;
	}
	return 0;
}

static SEXP to_integer_LIST(const IntAEAE *aeae, SEXP starts)
{
	int ndim, along;
	SEXP ans, ans_elt, start;

	ndim = LENGTH(starts);
	ans = PROTECT(NEW_LIST(ndim));
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		ans_elt = PROTECT(new_INTEGER_from_IntAE(aeae->elts[along]));
		SET_VECTOR_ELT(ans, along, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP to_numeric_LIST(const LLongAEAE *aeae, SEXP starts)
{
	int ndim, along;
	SEXP ans, ans_elt, start;
	const LLongAE *ae;
	R_xlen_t ans_elt_len, i;

	ndim = LENGTH(starts);
	ans = PROTECT(NEW_LIST(ndim));
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		ae = aeae->elts[along];
		ans_elt_len = LLongAE_get_nelt(ae);
		ans_elt = PROTECT(NEW_NUMERIC(ans_elt_len));
		for (i = 0; i < ans_elt_len; i++)
			REAL(ans_elt)[i] = (double) ae->elts[i];
		SET_VECTOR_ELT(ans, along, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Return a list of length 2:
 *   - The 1st list element is the list of break points.
 *   - The 2nd list element is the list of chunk indices.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
 */
SEXP C_map_starts_to_chunks(SEXP starts, SEXP dim, SEXP chunk_spacings)
{
	int ndim, along, ret;
	LLongAE *dim_buf, *chunk_spacings_buf;
	long long int d, chunkd;
	IntAE *nstart_buf;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *chunkidx_bufs;
	SEXP ans, ans_elt;

	ndim = _shallow_check_selection(starts, R_NilValue);
	if (ndim < 0)
		error(_HDF5Array_errmsg_buf);

	if (!(IS_INTEGER(dim) || IS_NUMERIC(dim)))
		error("'dim' must be an integer vector (or NULL)");
	if (LENGTH(dim) != ndim)
		error("'starts' and 'dim' must have the same length");
	if (!(IS_INTEGER(chunk_spacings) || IS_NUMERIC(chunk_spacings)))
		error("'chunk_spacings' must be an integer vector (or NULL)");
	if (LENGTH(chunk_spacings) != ndim)
		error("'starts' and 'chunk_spacings' must "
		      "have the same length");

	dim_buf = new_LLongAE(ndim, ndim, 0);
	chunk_spacings_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0; along < ndim; along++) {
		ret = get_untrusted_elt(dim, along, &d,
					"dim", -1);
		if (ret < 0)
			error(_HDF5Array_errmsg_buf);
		ret = get_untrusted_elt(chunk_spacings, along, &chunkd,
					"chunk_spacings", -1);
		if (ret < 0)
			error(_HDF5Array_errmsg_buf);
		if (chunkd < 0)
			error("'chunk_spacings' cannot "
			      "contain negative values");
		if (chunkd == 0 && d != 0)
			error("values in 'chunk_spacings' cannot be 0 unless "
			      "their corresponding value\n  in 'dim' is "
			      "also 0");
		dim_buf->elts[along] = d;
		chunk_spacings_buf->elts[along] = chunkd;
	}

	nstart_buf = new_IntAE(ndim, ndim, 0);
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunkidx_bufs = new_LLongAEAE(ndim, ndim);
	ret = _map_starts_to_chunks(starts,
			dim_buf->elts, chunk_spacings_buf->elts,
			nstart_buf->elts,
			breakpoint_bufs, chunkidx_bufs);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf);

	ans = PROTECT(NEW_LIST(2));
	ans_elt = PROTECT(to_integer_LIST(breakpoint_bufs, starts));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);
	ans_elt = PROTECT(to_numeric_LIST(chunkidx_bufs, starts));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(2);
	return ans;
}



