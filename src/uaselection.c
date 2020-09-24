/****************************************************************************
 *             Manipulation of a user-supplied array selection              *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "uaselection.h"

#include "global_errmsg_buf.h"

#include <limits.h>  /* for INT_MAX, LLONG_MAX, LLONG_MIN */
//#include <time.h>


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
		INTEGER(x)[i] = (int) val;
	else
		REAL(x)[i] = (double) val;
	return;
}

/* Called at the very beginning of the various .Call entry points where
   it's used (and before any resource is allocated) so it's ok to error()
   immediately in case of error. */
static const long long int *check_dim(SEXP dim)
{
	int ndim, i, ret;
	long long int *dim_p, d;

	if (!(IS_INTEGER(dim) || IS_NUMERIC(dim)))
		error("'dim' must be an integer vector");
	ndim = LENGTH(dim);
	dim_p = new_LLongAE(ndim, ndim, 0)->elts;
	for (i = 0; i < ndim; i++) {
		ret = get_untrusted_elt(dim, i, &d, "dim", -1);
		if (ret < 0)
			error(_HDF5Array_global_errmsg_buf());
		dim_p[i] = d;
	}
	return dim_p;
}


/****************************************************************************
 * Shallow check of a user-supplied array selection
 */

/* Only check that each of 'starts' and 'counts' is either NULL or a list
   of length as 'ndim'.
   Return 0 if the uaselection is valid and -1 if it's not. */
int _shallow_check_uaselection(int ndim, SEXP starts, SEXP counts)
{
	if (starts == R_NilValue) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF(
				"'counts' must be NULL when 'starts' is NULL");
			return -1;
		}
		return 0;
	}
	if (!isVectorList(starts)) {  // IS_LIST() is broken
		PRINT_TO_ERRMSG_BUF("'starts' must be a list (or NULL)");
		return -1;
	}
	if (LENGTH(starts) != ndim) {
		PRINT_TO_ERRMSG_BUF(
			"Array has %d dimension%s but 'starts' has %d "
			"list element%s.\n  'starts' must have one "
			"list element per dimension in the dataset.",
			ndim, ndim > 1 ? "s" : "",
			LENGTH(starts), LENGTH(starts) > 1 ? "s" : "");
		return -1;
	}
	if (counts == R_NilValue)
		return 0;
	if (!isVectorList(counts)) {  // IS_LIST() is broken
		PRINT_TO_ERRMSG_BUF("'counts' must be a list (or NULL)");
		return -1;
	}
	if (LENGTH(counts) != ndim) {
		PRINT_TO_ERRMSG_BUF("'counts' must have one list element "
				    "per list element in 'starts'");
		return -1;
	}
	return 0;
}


/****************************************************************************
 * Deep check of a user-supplied array selection
 */

static void set_error_for_uaselection_too_large(int along1)
{
	PRINT_TO_ERRMSG_BUF("too many elements (>= 2^31) selected "
			    "along dimension %d of array", along1);
	return;
}

static void set_errmsg_for_uaselection_beyond_dim(
		int along1, int i,
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

static void set_errmsg_for_non_strictly_ascending_uaselection(
		int along1, int i,
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
	if (*s < 1) {
		PRINT_TO_ERRMSG_BUF("starts[[%d]][%d] is < 1",
				    along + 1, i + 1);
		return -1;
	}
	if (*s < min_start) {
		set_errmsg_for_non_strictly_ascending_uaselection(
			along + 1, i, no_counts);
		return -1;
	}
	return 0;
}

static int check_uaselection_along(int along,
				  SEXP start, SEXP count, long long int d)
{
	long long int uaselection_dim, s, c, e;
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
				set_error_for_uaselection_too_large(along + 1);
				return -1;
			}
			uaselection_dim = d;
		} else {
			/* The dimension of the uaselection along the current
			   dimension is undefined in that case.
			   We **arbitrary** set it to INT_MAX. */
			uaselection_dim = INT_MAX;
		}
		return (int) uaselection_dim;
	}
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	n = LENGTH(start);
	if (shallow_check_count(count, n, along) < 0)
		return -1;
	/* Walk on the 'start' elements. */
	for (i = 0; i < n; i++) {
		/* Last arg ('no_counts') is ignored when 4th arg ('min_start')
		   is set to 0. */
		ret = get_untrusted_start(start, i, &s, 0, along, 0);
		if (ret < 0)
			return -1;
		if (d >= 0 && s > d) {
			set_errmsg_for_uaselection_beyond_dim(
				along + 1, i, 1);
			return -1;
		}
	}
	if (count == R_NilValue)
		return n;
	/* Walk on the 'count' (and 'start') elements. */
	uaselection_dim = 0;
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
			set_errmsg_for_uaselection_beyond_dim(
				along + 1, i, 0);
			return -1;
		}
		uaselection_dim += c; // could overflow! (FIXME)
		if (uaselection_dim > INT_MAX) {
			set_error_for_uaselection_too_large(along + 1);
			return -1;
		}
	}
	return (int) uaselection_dim;
}

/* 'dim' must be NULL or point to an array of 'ndim' elements.

   'starts' and 'counts' are **assumed** to be NULL or lists of length 'ndim'.
   This should have been already checked by _shallow_check_uaselection() so is
   not checked again.

   'uaselection_dim_buf' must point to an array of 'ndim' elements.
*/
long long int _check_uaselection(int ndim, const long long int *dim,
			SEXP starts, SEXP counts, int *uaselection_dim_buf)
{
	long long int uaselection_len;
	int along, uaselection_dim;
	SEXP start, count;

	uaselection_len = 1;
	for (along = 0; along < ndim; along++) {
		start = GET_LIST_ELT(starts, along);
		count = GET_LIST_ELT(counts, along);
		uaselection_dim = check_uaselection_along(along, start, count,
							dim[along]);
		if (uaselection_dim < 0)
			return -1;
		uaselection_dim_buf[along] = uaselection_dim;
		uaselection_len *= uaselection_dim;
	}
	return uaselection_len;
}

/* --- .Call ENTRY POINT ---
 * Return the dimensions of the user-supplied array selection.
 */
SEXP C_check_uaselection(SEXP dim, SEXP starts, SEXP counts)
{
	const long long int *dim_p;
	int ndim, ret;
	IntAE *uaselection_dim_buf;
	long long int uaselection_len;

	dim_p = check_dim(dim);
	ndim = LENGTH(dim);
	ret = _shallow_check_uaselection(ndim, starts, counts);
	if (ret < 0)
		error(_HDF5Array_global_errmsg_buf());

	uaselection_dim_buf = new_IntAE(ndim, ndim, 0);
	uaselection_len = _check_uaselection(ndim, dim_p, starts, counts,
					     uaselection_dim_buf->elts);
	if (uaselection_len < 0)
		error(_HDF5Array_global_errmsg_buf());
	return new_INTEGER_from_IntAE(uaselection_dim_buf);
}


/****************************************************************************
 * Deep check of an ordered user-supplied array selection (in preparation for
 * its reduction)
 *
 * The "chips" in the user-supplied array selection are its connected
 * components i.e. its contiguous block-like components.
 */

static int check_ordered_uaselection_along_NULL_start(int along,
			SEXP count, long long int d,
			int *nstart_buf, int *nchip_buf,
			long long int *last_chip_start_buf)
{
	int uaselection_dim;

	if (count != R_NilValue) {
		PRINT_TO_ERRMSG_BUF(
			"if 'starts[[%d]]' is NULL then 'counts' "
			"or 'counts[[%d]]' must also be NULL",
			along + 1, along + 1);
		return -1;
	}
	if (d >= 0) {
		if (d > INT_MAX) {
			set_error_for_uaselection_too_large(along + 1);
			return -1;
		}
		uaselection_dim = d;
		nstart_buf[along] = d;
		nchip_buf[along] = d != 0;
		last_chip_start_buf[along] = 1;
	} else {
		/* The dimension of the uaselection along the current
		   dimension is undefined in that case.
		   We **arbitrary** set it to INT_MAX. */
		uaselection_dim = INT_MAX;
		nstart_buf[along] = nchip_buf[along] = 1;
		last_chip_start_buf[along] = 1;
	}
	return uaselection_dim;
}

static int check_ordered_uaselection_along(int along,
			SEXP start, SEXP count, long long int d,
			int *nstart_buf, int *nchip_buf,
			long long int *last_chip_start_buf)
{
	long long int uaselection_dim, min_start, s, c;
	int n, i, ret;

	if (start == R_NilValue)
		return check_ordered_uaselection_along_NULL_start(along,
				count, d,
				nstart_buf, nchip_buf, last_chip_start_buf);
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	n = LENGTH(start);
	if (shallow_check_count(count, n, along) < 0)
		return -1;
	nstart_buf[along] = n;
	nchip_buf[along] = 0;
	min_start = 0;
	if (count == R_NilValue) {
		/* Walk on the 'start' elements. */
		for (i = 0; i < n; i++) {
			ret = get_untrusted_start(start, i, &s, min_start,
						  along, 1);
			if (ret < 0)
				return -1;
			if (s != min_start) {
				nchip_buf[along]++;
				last_chip_start_buf[along] = s;
			}
			min_start = s + 1;
			if (d >= 0 && s > d) {
				set_errmsg_for_uaselection_beyond_dim(
					along + 1, i, 1);
				return -1;
			}
		}
		uaselection_dim = n;
	} else {
		/* Walk on the 'start' and 'count' elements. */
		uaselection_dim = 0;
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
				nchip_buf[along]++;
				last_chip_start_buf[along] = s;
			}
			min_start = s + c;  // could overflow! (FIXME)
			if (d >= 0 && min_start - 1 > d) {
				set_errmsg_for_uaselection_beyond_dim(
					along + 1, i, 0);
				return -1;
			}
			uaselection_dim += c; // could overflow! (FIXME)
			if (uaselection_dim > INT_MAX) {
				set_error_for_uaselection_too_large(along + 1);
				return -1;
			}
		}
	}
	return (int) uaselection_dim;
}

/* 'dim' must be NULL or point to an array of 'ndim' elements.

   'starts' and 'counts' are **assumed** to be NULL or lists of length 'ndim'.
   This should have been already checked by _shallow_check_uaselection() so is
   not checked again.

   Each of 'uaselection_dim_buf', 'nstart_buf', 'nchip_buf', and
   'last_chip_start_buf' must point to an array of 'ndim' elements.
*/
long long int _check_ordered_uaselection(int ndim, const long long int *dim,
			SEXP starts, SEXP counts, int *uaselection_dim_buf,
			int *nstart_buf, int *nchip_buf,
			long long int *last_chip_start_buf)
{
	long long int uaselection_len;
	int along, uaselection_dim;
	SEXP start, count;

	uaselection_len = 1;
	for (along = 0; along < ndim; along++) {
		start = GET_LIST_ELT(starts, along);
		count = GET_LIST_ELT(counts, along);
		uaselection_dim = check_ordered_uaselection_along(along,
					start, count,
					dim != NULL ? dim[along] : -1,
					nstart_buf, nchip_buf,
					last_chip_start_buf);
		if (uaselection_dim < 0)
			return -1;
		uaselection_dim_buf[along] = uaselection_dim;
		uaselection_len *= uaselection_dim;
	}
	return uaselection_len;
}

/* --- .Call ENTRY POINT ---
 * Return the dimensions of the user-supplied array selection.
 */
SEXP C_check_ordered_uaselection(SEXP dim, SEXP starts, SEXP counts)
{
	const long long int *dim_p;
	int ndim, ret;
	IntAE *uaselection_dim_buf, *nstart_buf, *nchip_buf;
	LLongAE *last_chip_start_buf;
	long long int uaselection_len;

	dim_p = check_dim(dim);
	ndim = LENGTH(dim);
	ret = _shallow_check_uaselection(ndim, starts, counts);
	if (ret < 0)
		error(_HDF5Array_global_errmsg_buf());

	uaselection_dim_buf = new_IntAE(ndim, ndim, 0);
	nstart_buf = new_IntAE(ndim, ndim, 0);
	nchip_buf = new_IntAE(ndim, ndim, 0);
	last_chip_start_buf = new_LLongAE(ndim, ndim, 0);
	uaselection_len = _check_ordered_uaselection(ndim, dim_p,
				starts, counts,
				uaselection_dim_buf->elts,
				nstart_buf->elts, nchip_buf->elts,
				last_chip_start_buf->elts);
	if (uaselection_len < 0)
		error(_HDF5Array_global_errmsg_buf());
	return new_INTEGER_from_IntAE(uaselection_dim_buf);
}


/****************************************************************************
 * Reduce the user-supplied array selection
 */

int _uaselection_can_be_reduced(int ndim, const int *nstart, const int *nchip)
{
	int along;

	for (along = 0; along < ndim; along++) {
		/* nchip[along] should always be <= nstart[along] */
		if (nchip[along] < nstart[along])
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
static void stitch_uaselection(SEXP start_in, SEXP count_in,
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
			min_start = s + 1;
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

static void reduce_uaselection_along(int along,
				    SEXP start, SEXP count,
				    const int *uaselection_dim,
				    const int *nchip,
				    const long long int *last_chip_start,
				    SEXP reduced_starts, SEXP reduced_counts)
{
	int n, dup;
	SEXP reduced_start, reduced_count;
	SEXPTYPE type;

	n = LENGTH(start);
	if (nchip[along] == n) {
		/* Nothing to stitch. */
		dup = IS_INTEGER(start) || last_chip_start[along] > INT_MAX;
		reduced_start = PROTECT(dup_or_coerce_to_INTSXP(start, dup));
		SET_VECTOR_ELT(reduced_starts, along, reduced_start);
		UNPROTECT(1);
		if (uaselection_dim[along] == n)
			return;
		dup = IS_INTEGER(count);
		reduced_count = PROTECT(dup_or_coerce_to_INTSXP(count, dup));
		SET_VECTOR_ELT(reduced_counts, along, reduced_count);
		UNPROTECT(1);
		return;
	}
	/* Stitch. */
	type = last_chip_start[along] <= INT_MAX ? INTSXP : REALSXP;
	reduced_start = PROTECT(allocVector(type, nchip[along]));
	SET_VECTOR_ELT(reduced_starts, along, reduced_start);
	UNPROTECT(1);
	reduced_count = PROTECT(NEW_INTEGER(nchip[along]));
	SET_VECTOR_ELT(reduced_counts, along, reduced_count);
	UNPROTECT(1);
	stitch_uaselection(start, count, reduced_start, INTEGER(reduced_count));
	return;
}

SEXP _reduce_uaselection(int ndim, SEXP starts, SEXP counts,
			 const int *uaselection_dim,
			 const int *nchip,
			 const long long int *last_chip_start)
{
	SEXP reduced_starts, reduced_counts, start, count, ans;
	int along;

	//clock_t t0 = clock();
	reduced_starts = PROTECT(NEW_LIST(ndim));
	reduced_counts = PROTECT(NEW_LIST(ndim));
	if (starts != R_NilValue) {
		for (along = 0; along < ndim; along++) {
			start = VECTOR_ELT(starts, along);
			if (start == R_NilValue)
				continue;
			count = GET_LIST_ELT(counts, along);
			reduce_uaselection_along(along,
					start, count,
					uaselection_dim,
					nchip, last_chip_start,
					reduced_starts, reduced_counts);
		}
	}
	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, reduced_starts);
	SET_VECTOR_ELT(ans, 1, reduced_counts);
	UNPROTECT(3);
	//printf("time 2nd pass: %e\n", (1.0 * clock() - t0) / CLOCKS_PER_SEC);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Negative values in 'dim' are treated as infinite dimensions.
 * Return a list of length 2 or NULL if the user-supplied array selection
 * could not be reduced.
 * When returning a list of length 2:
 *   - The 1st list element is the list of reduced starts.
 *   - The 2nd list element is the list of reduced counts.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
 */
SEXP C_reduce_uaselection(SEXP dim, SEXP starts, SEXP counts)
{
	const long long int *dim_p;
	int ndim, ret;
	IntAE *uaselection_dim_buf, *nstart_buf, *nchip_buf;
	LLongAE *last_chip_start_buf;
	long long int uaselection_len;

	dim_p = check_dim(dim);
	ndim = LENGTH(dim);
	ret = _shallow_check_uaselection(ndim, starts, counts);
	if (ret < 0)
		error(_HDF5Array_global_errmsg_buf());

	uaselection_dim_buf = new_IntAE(ndim, ndim, 0);
	nstart_buf = new_IntAE(ndim, ndim, 0);
	nchip_buf = new_IntAE(ndim, ndim, 0);
	last_chip_start_buf = new_LLongAE(ndim, ndim, 0);

	/* 1st pass */
	uaselection_len = _check_ordered_uaselection(ndim, dim_p,
				starts, counts,
				uaselection_dim_buf->elts,
				nstart_buf->elts, nchip_buf->elts,
				last_chip_start_buf->elts);
	if (uaselection_len < 0)
		error(_HDF5Array_global_errmsg_buf());
	if (!_uaselection_can_be_reduced(ndim,
				       nstart_buf->elts,
				       nchip_buf->elts))
		return R_NilValue;

	/* 2nd pass */
	return _reduce_uaselection(ndim, starts, counts,
				   uaselection_dim_buf->elts, nchip_buf->elts,
				   last_chip_start_buf->elts);
}


/****************************************************************************
 * Map the user-supplied array selection to the physical chunks
 */

static int map_start_to_chunks(int along,
		long long int d, long long int chunkd, SEXP start,
		int *nstart_buf,
		IntAE *breakpoint_buf, LLongAE *tchunkidx_buf)
{
	int n, i, ret;
	size_t ntchunk;
	long long int min_start, s, tchunkidx, prev_tchunkidx;

	if (start == R_NilValue) {
		if (d > INT_MAX) {
			set_error_for_uaselection_too_large(along + 1);
			return -1;
		}
		nstart_buf[along] = d;
		return 0;
	}

	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;

	if (IntAE_get_nelt(breakpoint_buf) != 0 ||
	    LLongAE_get_nelt(tchunkidx_buf) != 0) {
		/* Should never happen! */
		PRINT_TO_ERRMSG_BUF("internal error: map_start_to_chunks() "
				    "was called with non-empty breakpoint "
				    "or tchunkidx buffers");
		return -1;
	}

	n = LENGTH(start);
	nstart_buf[along] = n;

	if (n == 0)
		return 0;

	/* Get 's' and 'tchunkidx' for 1st 'start' element. */
	ret = get_untrusted_start(start, 0, &s, 1, along, 1);
	if (ret < 0)
		return -1;
	if (s > d) {
		set_errmsg_for_uaselection_beyond_dim(along + 1, 0, 1);
		return -1;
	}
	tchunkidx = (s - 1) / chunkd;

	/* Walk on the remaining 'start' elements. */
	ntchunk = 0;
	for (i = 1; i < n; i++) {
		min_start = s + 1;
		ret = get_untrusted_start(start, i, &s, min_start, along, 1);
		if (ret < 0)
			return -1;
		if (s > d) {
			set_errmsg_for_uaselection_beyond_dim(along + 1, i, 1);
			return -1;
		}
		prev_tchunkidx = tchunkidx;
		tchunkidx = (s - 1) / chunkd;
		if (tchunkidx > prev_tchunkidx) {
			IntAE_insert_at(breakpoint_buf, ntchunk, i);
			LLongAE_insert_at(tchunkidx_buf, ntchunk,
					  prev_tchunkidx);
			ntchunk++;
		}
	}
	IntAE_insert_at(breakpoint_buf, ntchunk, n);
	LLongAE_insert_at(tchunkidx_buf, ntchunk, tchunkidx);
	return 0;
}

/* 'dim', 'chunkdim', and 'nstart_buf' must point to arrays of 'ndim' elements.

   'starts' is **assumed** to be NULL or a list of length 'ndim'.
   This should have been already checked by _shallow_check_uaselection() so is
   not checked again.

   'breakpoint_bufs' and 'tchunkidx_bufs' must be of length 'ndim'.
*/
int _map_starts_to_chunks(int ndim, const long long int *dim,
		const long long int *chunkdim,
		SEXP starts,
		int *nstart_buf,
		IntAEAE *breakpoint_bufs, LLongAEAE *tchunkidx_bufs)
{
	int along, ret;
	SEXP start;

	for (along = 0; along < ndim; along++) {
		start = GET_LIST_ELT(starts, along);
		ret = map_start_to_chunks(along,
					  dim[along], chunkdim[along], start,
					  nstart_buf,
					  breakpoint_bufs->elts[along],
					  tchunkidx_bufs->elts[along]);
		if (ret < 0)
			return -1;
	}
	return 0;
}

static SEXP to_integer_LIST(int ndim, const IntAEAE *aeae, SEXP starts)
{
	SEXP ans, ans_elt, start;
	int along;
	const IntAE *ae;

	ans = PROTECT(NEW_LIST(ndim));
	if (starts != R_NilValue) {
		for (along = 0; along < ndim; along++) {
			start = VECTOR_ELT(starts, along);
			if (start == R_NilValue)
				continue;
			ae = aeae->elts[along];
			ans_elt = PROTECT(new_INTEGER_from_IntAE(ae));
			SET_VECTOR_ELT(ans, along, ans_elt);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}

static SEXP to_numeric_LIST(int ndim, const LLongAEAE *aeae, SEXP starts)
{
	SEXP ans, ans_elt, start;
	int along;
	const LLongAE *ae;
	R_xlen_t ans_elt_len, i;

	ans = PROTECT(NEW_LIST(ndim));
	if (starts != R_NilValue) {
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
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Return a list of length 2:
 *   - The 1st list element is the list of break points along each dim.
 *   - The 2nd list element is the list of touched chunk ids along each dim.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
 */
SEXP C_map_starts_to_chunks(SEXP starts, SEXP dim, SEXP chunkdim)
{
	const long long int *dim_p;
	int ndim, ret, i;
	LLongAE *chunkdim_buf;
	long long int chunkd;
	IntAE *nstart_buf;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;
	SEXP ans, ans_elt;

	dim_p = check_dim(dim);
	ndim = LENGTH(dim);
	ret = _shallow_check_uaselection(ndim, starts, R_NilValue);
	if (ret < 0)
		error(_HDF5Array_global_errmsg_buf());

	if (!(IS_INTEGER(chunkdim) || IS_NUMERIC(chunkdim)))
		error("'chunkdim' must be an integer vector");
	if (LENGTH(chunkdim) != ndim)
		error("'chunkdim' must have the same length as 'dim'");
	chunkdim_buf = new_LLongAE(ndim, ndim, 0);
	for (i = 0; i < ndim; i++) {
		ret = get_untrusted_elt(chunkdim, i, &chunkd, "chunkdim", -1);
		if (ret < 0)
			error(_HDF5Array_global_errmsg_buf());
		if (chunkd < 0)
			error("'chunkdim' cannot contain negative values");
		if (chunkd == 0 && dim_p[i] != 0)
			error("values in 'chunkdim' cannot be 0 unless "
			      "their corresponding value\n  in 'dim' is "
			      "also 0");
		chunkdim_buf->elts[i] = chunkd;
	}

	nstart_buf = new_IntAE(ndim, ndim, 0);
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	tchunkidx_bufs = new_LLongAEAE(ndim, ndim);
	ret = _map_starts_to_chunks(ndim, dim_p, chunkdim_buf->elts,
			starts,
			nstart_buf->elts,
			breakpoint_bufs, tchunkidx_bufs);
	if (ret < 0)
		error(_HDF5Array_global_errmsg_buf());

	ans = PROTECT(NEW_LIST(2));
	ans_elt = PROTECT(to_integer_LIST(ndim, breakpoint_bufs, starts));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);
	ans_elt = PROTECT(to_numeric_LIST(ndim, tchunkidx_bufs, starts));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(2);
	return ans;
}

