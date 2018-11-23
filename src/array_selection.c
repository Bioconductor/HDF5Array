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
 * Shallow check the array selection
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
 * Deep check the array selection
 */

static void set_error_for_selection_too_large(int along1)
{
	PRINT_TO_ERRMSG_BUF("too many elements (>= 2^31) selected "
			    "along dimension %d of array", along1);
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
			  "array, but you\n  have:";
	if (starts_only)
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

static inline int get_untrusted_start(SEXP start, int i, long long int *s,
				      int along,
				      long long int e,
				      int starts_only)
{
	if (get_untrusted_elt(start, i, s, "starts", along) < 0)
		return -1;
	if (*s > e)
		return 0;
	if (e == 0) {
		PRINT_TO_ERRMSG_BUF("starts[[%d]][%d] is <= 0",
				    along + 1, i + 1);
	} else {
		set_errmsg_for_non_strictly_ascending_selection(
			along + 1, i, starts_only);
	}
	return -1;
}

static int check_selection_along_NULL_start(SEXP count, int along,
			long long int d,
			int *nstart, int *count_sum,
			int *nblock, long long int *last_block_start)
{
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

static int check_selection_along(SEXP start, SEXP count, int along,
			long long int d,
			int *nstart, int *count_sum,
			int *nblock, long long int *last_block_start)
{
	int n, i, ret;
	long long int e, s, cs, c;

	if (start == R_NilValue)
		return check_selection_along_NULL_start(count, along,
				d,
				nstart, count_sum,
				nblock, last_block_start);
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	n = LENGTH(start);
	if (count != R_NilValue) {
		if (check_INTEGER_or_NUMERIC(count, "counts", along) < 0)
			return -1;
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
			ret = get_untrusted_start(start, i, &s, along, e, 1);
			if (ret < 0)
				return -1;
			if (i == 0 || s != e + 1) {
				nblock[along]++;
				last_block_start[along] = s;
			}
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
			ret = get_untrusted_start(start, i, &s, along, e, 0);
			if (ret < 0)
				return -1;
			if (i == 0 || s != e + 1) {
				nblock[along]++;
				last_block_start[along] = s;
			}
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

/* Assume that 'starts' is a list and that 'counts' is NULL or a list
   of the same length as 'starts'. This should have been checked
   by _shallow_check_selection() already so is not checked again.

   'dim' is assumed to be NULL or to have the same length as 'starts'.

   'nstart', 'count_sum', 'nblock', and 'last_block_start' are assumed
   to have the same length as 'starts'.
*/
int _deep_check_selection(SEXP starts, SEXP counts,
			  const long long int *dim,
			  int *nstart, int *count_sum,
			  int *nblock, long long int *last_block_start)
{
	int ndim, along, ret;
	SEXP start, count;

	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		count = counts != R_NilValue ? VECTOR_ELT(counts, along)
					     : R_NilValue;
		ret = check_selection_along(start, count, along,
					    dim != NULL ? dim[along] : -1,
					    nstart, count_sum,
					    nblock, last_block_start);
		if (ret < 0)
			return -1;
	}
	return 0;
}


/****************************************************************************
 * Map array selection to chunks
 */

static int map_start_to_chunks(SEXP start, int along,
		long long int d, long long int chunkd,
		int *nstart,
		IntAE *breakpoint_buf, IntAE *chunkidx_buf)
{
	int n, i, ret;
	size_t buf_nelt;
	long long int e, s, chunk_idx, prev_chunk_idx;

	if (start == R_NilValue) {
		if (d > INT_MAX) {
			set_error_for_selection_too_large(along + 1);
			return -1;
		}
		nstart[along] = d;
		return 0;
	}
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	if (IntAE_get_nelt(breakpoint_buf) != 0 ||
	    IntAE_get_nelt(chunkidx_buf) != 0) {
		/* Should never happen! */
		PRINT_TO_ERRMSG_BUF("internal error: map_start_to_chunks() "
				    "was called with non-empty breakpoint_buf "
				    "or chunkidx_buf");
		return -1;
	}
	n = LENGTH(start);
	nstart[along] = n;
	buf_nelt = 0;
	e = 0;
	prev_chunk_idx = -1;
	for (i = 0; i < n; i++) {
		ret = get_untrusted_start(start, i, &s, along, e, 1);
		if (ret < 0)
			return -1;
		e = s;
		if (e > d) {
			set_errmsg_for_selection_beyond_dim(along + 1, i, 1);
			return -1;
		}
		chunk_idx = (s - 1) / chunkd;
		if (chunk_idx > INT_MAX) {
			/* If we wanted to support this then we'd need to use
			   something like a LLongAE for 'chunkidx_buf'. */
			PRINT_TO_ERRMSG_BUF("selecting array elements from a "
					    "chunk whose index is > INT_MAX "
					    "is not\n  supported "
					    "at the moment");
			return -1;
		}
		if (chunk_idx > prev_chunk_idx) {
			IntAE_insert_at(breakpoint_buf, buf_nelt, i + 1);
			IntAE_insert_at(chunkidx_buf, buf_nelt, chunk_idx);
			buf_nelt++;
		} else {
			breakpoint_buf->elts[buf_nelt - 1]++;
		}
		prev_chunk_idx = chunk_idx;
	}
	return 0;
}

/* Assume that 'starts' is a list. This should have been checked by
   _shallow_check_selection() already so is not checked again.

   'dim', 'chunkdim', 'nstart', 'breakpoint_bufs', and 'chunkidx_bufs'
   are assumed to have the same length as 'starts'.

   We're storing the chunk indices in a IntAEAE struct which means that we
   cannot support arrays with more than INT_MAX chunks along any dimension.
   If we wanted to support such arrays then we'd need to use something like
   a LLongAEAE instead. (IntAE, IntAEAE, and LLongAE structs are implemented
   in S4Vectors/src/AEbufs.c but LLongAEAE is not implemented yet.)
*/
int _map_starts_to_chunks(SEXP starts,
		const long long int *dim,
		const long long int *chunkdim,
		int *nstart,
		IntAEAE *breakpoint_bufs, IntAEAE *chunkidx_bufs)
{
	int ndim, along, ret;
	SEXP start;

	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		ret = map_start_to_chunks(start, along,
					  dim[along], chunkdim[along],
					  nstart,
					  breakpoint_bufs->elts[along],
					  chunkidx_bufs->elts[along]);
		if (ret < 0)
			return -1;
	}
	return 0;
}

static SEXP as_LIST(const IntAEAE *buf, SEXP starts)
{
	int ndim, along;
	SEXP ans, ans_elt, start;

	ndim = LENGTH(starts);
	ans = PROTECT(NEW_LIST(ndim));
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		ans_elt = PROTECT(new_INTEGER_from_IntAE(buf->elts[along]));
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
SEXP C_map_starts_to_chunks(SEXP starts, SEXP dim, SEXP chunkdim)
{
	int ndim, along, ret;
	LLongAE *dim_buf, *chunkdim_buf;
	long long int d, chunkd;
	IntAE *nstart_buf;
	IntAEAE *breakpoint_bufs, *chunkidx_bufs;
	SEXP ans, ans_elt;

	ndim = _shallow_check_selection(starts, R_NilValue);
	if (ndim < 0)
		error(_HDF5Array_errmsg_buf);

	if (!(IS_INTEGER(dim) || IS_NUMERIC(dim)))
		error("'dim' must be an integer vector (or NULL)");
	if (LENGTH(dim) != ndim)
		error("'starts' and 'dim' must have the same length");
	if (!(IS_INTEGER(chunkdim) || IS_NUMERIC(chunkdim)))
		error("'chunkdim' must be an integer vector (or NULL)");
	if (LENGTH(chunkdim) != ndim)
		error("'starts' and 'chunkdim' must have the same length");

	dim_buf = new_LLongAE(ndim, ndim, 0);
	chunkdim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0; along < ndim; along++) {
		ret = get_untrusted_elt(dim, along, &d,
					"dim", -1);
		if (ret < 0)
			error(_HDF5Array_errmsg_buf);
		ret = get_untrusted_elt(chunkdim, along, &chunkd,
					"chunkdim", -1);
		if (ret < 0)
			error(_HDF5Array_errmsg_buf);
		if (chunkd < 0)
			error("'chunkdim' cannot contain negative values");
		if (chunkd == 0 && d != 0)
			error("values in 'chunkdim' cannot be 0 unless their "
			      "corresponding value\n  in 'dim' is also 0");
		dim_buf->elts[along] = d;
		chunkdim_buf->elts[along] = chunkd;
	}

	nstart_buf = new_IntAE(ndim, ndim, 0);
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunkidx_bufs = new_IntAEAE(ndim, ndim);
	ret = _map_starts_to_chunks(starts,
			dim_buf->elts, chunkdim_buf->elts,
			nstart_buf->elts,
			breakpoint_bufs, chunkidx_bufs);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf);

	ans = PROTECT(NEW_LIST(2));
	ans_elt = PROTECT(as_LIST(breakpoint_bufs, starts));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);
	ans_elt = PROTECT(as_LIST(chunkidx_bufs, starts));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(2);
	return ans;
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

static void reduce_selection_along(SEXP start, SEXP count, int along,
				   const int *count_sum,
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
 * When returning a list of length 2:
 *   - The 1st list element is the list of reduced starts.
 *   - The 2nd list element is the list of reduced counts.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
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

