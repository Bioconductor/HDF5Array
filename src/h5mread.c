/****************************************************************************
 *       Experimenting with alternate rhdf5::h5read() implementations       *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"
#include "S4Vectors_interface.h"
#include "hdf5.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcpy */
#include <limits.h>  /* for INT_MAX, LLONG_MAX, LLONG_MIN */

//#include <time.h>

static char errmsg_buf[256];

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
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "%s[%d] is NA",
					 what, i + 1);
			else
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "%s[[%d]][%d] is NA",
					 what, along + 1, i + 1);
			return -1;
		}
		*val = (long long int) tmp1;
	} else {
		tmp2 = REAL(x)[i];
		if (NOT_A_FINITE_NUMBER(tmp2)) {
			if (along < 0)
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "%s[%d] is NA or NaN "
					 "or not a finite number",
					 what, i + 1);
			else
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "%s[[%d]][%d] is NA or NaN "
					 "or not a finite number",
					 what, along + 1, i + 1);
			return -1;
		}
		if (tmp2 > (double) LLONG_MAX ||
		    tmp2 < (double) LLONG_MIN) {
			if (along < 0)
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "%s[%d] is too large (= %e)",
					 what, i + 1, tmp2);
			else
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "%s[[%d]][%d] is too large (= %e)",
					 what, along + 1, i + 1, tmp2);
			return -1;
		}
		*val = (long long int) tmp2;
	}
	return 0;
}

/* Unlike get_untrusted_elt() above, doesn't check anything. */
static inline long long int get_trusted_elt(SEXP x, int i)
{
	return IS_INTEGER(x) ? (long long int) INTEGER(x)[i] :
			       (long long int) REAL(x)[i];
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
 * Low-level helpers for checking the selection
 */

static int shallow_check_starts_counts(SEXP starts, SEXP counts)
{
	int ndim;

	if (!isVectorList(starts)) {  // IS_LIST() is broken
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'starts' must be a list");
		return -1;
	}
	ndim = LENGTH(starts);
	if (counts != R_NilValue) {
		if (!isVectorList(counts)) {  // IS_LIST() is broken
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'counts' must be a list (or NULL)");
			return -1;
		}
		if (LENGTH(counts) != ndim) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'starts' and 'counts' "
				 "must have the same length");
			return -1;
		}
	}
	return ndim;
}

static void set_error_for_selection_too_large(int along1)
{
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "too many elements (>= 2^31) selected "
		 "along dimension %d of dataset",
		 along1);
	return;
}

static void set_errmsg_for_non_strictly_ascending_selection(int along1, int i,
							    int starts_only)
{
	const char *msg = "selection must be strictly ascending "
			  "along each dimension, but\n  you have:";
	if (starts_only)
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] <= starts[[%d]][%d]",
			 msg, along1, i + 1, along1, i);
	else
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] < starts[[%d]][%d] + "
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
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] "
			 "> dimension %d in dataset",
			 msg, along1, i + 1, along1);
	else
		snprintf(errmsg_buf, sizeof(errmsg_buf),
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
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "starts[[%d]][%d] is <= 0",
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

static int check_starts_counts_along(int along, SEXP starts, SEXP counts,
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
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "if 'starts[[%d]]' is NULL then 'counts' "
				 "or 'counts[[%d]]' must also be NULL",
				 along + 1, along + 1);
			return -1;
		}
		if (d > INT_MAX) {
			set_error_for_selection_too_large(along + 1);
			return -1;
		}
		nstart[along] = 1;
		count_sum[along] = d;
		nblock[along] = d != 0;
		last_block_start[along] = 1;
		return 0;
	}
	if (!(IS_INTEGER(start) || IS_NUMERIC(start))) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'starts[[%d]]' must be "
			 "an integer vector (or NULL)", along + 1);
		return -1;
	}
	n = LENGTH(start);
	if (count != R_NilValue) {
		if (!(IS_INTEGER(count) || IS_NUMERIC(count))) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'counts[[%d]]' must be an integer "
				 "vector (or NULL)", along + 1);
			return -1;
		}
		if (LENGTH(count) != n) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
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
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "counts[[%d]][%d] is <= 0",
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

static int check_starts_counts(SEXP starts, SEXP counts,
			       const long long int *dim,
			       int *nstart, int *count_sum,
			       int *nblock, long long int *last_block_start)
{
	int ndim, along, ret;

	ndim = LENGTH(starts);
	for (along = 0; along < ndim; along++) {
		ret = check_starts_counts_along(along, starts, counts,
						dim != NULL ? dim[along] : -1,
						nstart, count_sum,
						nblock, last_block_start);
		if (ret < 0)
			return -1;
	}
	return 0;
}


/****************************************************************************
 * C_reduce_selection()
 */

static int selection_can_be_reduced(int ndim,
				    const int *nstart, const int *nblock)
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
			s = get_trusted_elt(start_in, i);
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
			s = get_trusted_elt(start_in, i);
			c = get_trusted_elt(count_in, i);
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

static SEXP reduce_selection(SEXP starts, SEXP counts,
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

	ndim = shallow_check_starts_counts(starts, counts);
	if (ndim < 0)
		error(errmsg_buf);
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
				error(errmsg_buf);
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
	ret = check_starts_counts(starts, counts, dim_p,
				  nstart, count_sum,
				  nblock, last_block_start);
	//printf("time 1st pass: %e\n", (1.0 * clock() - t0) / CLOCKS_PER_SEC);
	if (ret < 0)
		error(errmsg_buf);

	if (!selection_can_be_reduced(ndim, nstart, nblock))
		return R_NilValue;

	/* 2nd pass */
	return reduce_selection(starts, counts,
				count_sum, nblock, last_block_start);
}


/****************************************************************************
 * C_h5mread()
 *
 * Some useful links:
 * - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
 *     https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
 * - Documentation of H5Dread():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
 * - A useful example:
 *     https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
*/

static int get_ans_type(hid_t dset_id, SEXPTYPE *ans_type)
{
	hid_t type_id;
	H5T_class_t class;
	const char *classname;

	type_id = H5Dget_type(dset_id);
	if (type_id < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Dget_type() returned an error");
		return -1;
	}
	class = H5Tget_class(type_id);
	H5Tclose(type_id);
	if (class == H5T_NO_CLASS) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Tget_class() returned an error");
		return -1;
	}
	if (class == H5T_INTEGER) {
		*ans_type = INTSXP;
		return 0;
	}
	if (class == H5T_FLOAT) {
		*ans_type = REALSXP;
		return 0;
	}
	switch (class) {
		case H5T_TIME: classname = "H5T_TIME"; break;
		case H5T_STRING: classname = "H5T_STRING"; break;
		case H5T_BITFIELD: classname = "H5T_BITFIELD"; break;
		case H5T_OPAQUE: classname = "H5T_OPAQUE"; break;
		case H5T_COMPOUND: classname = "H5T_COMPOUND"; break;
		case H5T_REFERENCE: classname = "H5T_REFERENCE"; break;
		case H5T_ENUM: classname = "H5T_ENUM"; break;
		case H5T_VLEN: classname = "H5T_VLEN"; break;
		case H5T_ARRAY: classname = "H5T_ARRAY"; break;
		default:
		    snprintf(errmsg_buf, sizeof(errmsg_buf),
			     "unknown dataset class identifier: %d", class);
		return -1;
	}
	snprintf(errmsg_buf, sizeof(errmsg_buf),
		 "unsupported dataset class: %s", classname);
	return -1;
}

static int deep_check_starts_counts(int ndim, const hsize_t *dset_dims,
			SEXP starts, SEXP counts,
			int *nstart, int *ans_dim,
			int *nblock, long long int *last_block_start)
{
	LLongAE *dim_buf;
	int along, h5along, ret;

	/* We already know that 'starts' is a list. */
	if (LENGTH(starts) != ndim) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'starts' must have one list element "
			 "per dimension in the dataset");
		return -1;
	}

	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] = (long long int) dset_dims[h5along];

	ret = check_starts_counts(starts, counts, dim_buf->elts,
				  nstart, ans_dim,
				  nblock, last_block_start);
	return ret;
}

static int next_midx(int ndim, const int *nstart, int *midx)
{
	int along, i;

	for (along = 0; along < ndim; along++) {
		i = midx[along] + 1;
		if (i < nstart[along]) {
			midx[along] = i;
			return 1;
		}
		midx[along] = 0;
	}
	return 0;
}

static int select_hyperslab(hid_t file_space_id, int ndim,
			    SEXP starts, SEXP counts,
			    const int *midx,
			    hsize_t *offset_buf, hsize_t *count_buf)
{
	int along, h5along, i, ret;
	SEXP start, count;

	/* Set 'offset_buf' and 'count_buf'. */
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		offset_buf[h5along] = get_trusted_elt(start, i) - 1;
		if (counts != R_NilValue) {
			count = VECTOR_ELT(counts, along);
			if (count != R_NilValue)
				count_buf[h5along] = get_trusted_elt(count, i);
		}
	}

	/* Add to current selection. */
	ret = H5Sselect_hyperslab(file_space_id, H5S_SELECT_OR,
				  offset_buf, NULL, count_buf, NULL);
	if (ret < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

/* Return nb of hyperslabs (or -1 if error). */
static long long int select_hyperslabs(hid_t file_space_id,
			int ndim, const hsize_t *dset_dims,
			SEXP starts, SEXP counts, const int *nstart)
{
	long long int num_hyperslabs;
	int along, h5along, ret;
	hsize_t *offset_buf, *count_buf;  // hyperslab offsets and dims
	IntAE *midx_buf;

	for (along = 0; along < ndim; along++)
		if (nstart[along] == 0)
			return 0;  // empty region (no hyperslab)

	/* Allocate 'offset_buf' and 'count_buf'. */
	offset_buf = (hsize_t *) malloc(2 * ndim * sizeof(hsize_t));
	if (offset_buf == NULL) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'offset_buf' "
			 "and 'count_buf'");
		return -1;
	}
	count_buf = offset_buf + ndim;

	/* Initialize 'offset_buf' and 'count_buf'. */
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			offset_buf[h5along] = 0;
			count_buf[h5along] = dset_dims[h5along];
		} else {
			count_buf[h5along] = 1;
		}
	}

	/* Walk on the blocks. */
	midx_buf = new_IntAE(ndim, ndim, 0);
	num_hyperslabs = 0;
	do {
		num_hyperslabs++;
		ret = select_hyperslab(file_space_id, ndim,
				       starts, counts,
				       midx_buf->elts,
				       offset_buf, count_buf);
		if (ret < 0)
			break;
	} while (next_midx(ndim, nstart, midx_buf->elts));

	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(nstart)
	free(offset_buf);
	return num_hyperslabs;
}

/* Return nb of elements (or -1 if error). */
static long long int select_elements(hid_t file_space_id,
			int ndim, const hsize_t *dset_dims,
			SEXP starts, const int *ans_dim)
{
	size_t num_elements;
	int along, i, ret;
	hsize_t *coord_buf, *coord_p;
	IntAE *midx_buf;
	SEXP start;
	long long int coord;

	num_elements = 1;
	for (along = 0; along < ndim; along++)
		num_elements *= ans_dim[along];
	//printf("nb of elements = %lu\n", num_elements);  // = length(ans)
	if (num_elements == 0)
		return 0;  // no elements

	/* Allocate 'coord_buf'. */
	coord_buf = (hsize_t *) malloc(num_elements * ndim * sizeof(hsize_t));
	if (coord_buf == NULL) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'coord_buf'");
		return -1;
	}

	/* Walk on the elements. */
	midx_buf = new_IntAE(ndim, ndim, 0);
	coord_p = coord_buf;
	do {
		for (along = ndim - 1; along >= 0; along--) {
			i = midx_buf->elts[along];
			start = VECTOR_ELT(starts, along);
			if (start == R_NilValue) {
				coord = i;
			} else {
				coord = get_trusted_elt(start, i) - 1;
			}
			*(coord_p++) = (hsize_t) coord;
		}
	} while (next_midx(ndim, ans_dim, midx_buf->elts));

	ret = H5Sselect_elements(file_space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	free(coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

static hid_t prepare_file_space(hid_t dset_id, SEXP starts, SEXP counts,
				int *ans_dim, int noreduce)
{
	hid_t file_space_id;
	int ndim, ret;
	hsize_t *dset_dims;
	IntAE *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	int *nstart, *nblock;
	long long int *last_block_start, nselection;
	SEXP reduced;

	file_space_id = H5Dget_space(dset_id);
	if (file_space_id < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Dget_space() returned an error");
		return -1;
	}

	ndim = H5Sget_simple_extent_ndims(file_space_id);

	/* Allocate and set 'dset_dims'. */
	dset_dims = (hsize_t *) malloc(ndim * sizeof(hsize_t));
	if (dset_dims == NULL) {
		H5Sclose(file_space_id);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'dset_dims'");
		return -1;
	}
	if (H5Sget_simple_extent_dims(file_space_id, dset_dims, NULL) != ndim) {
		free(dset_dims);
		H5Sclose(file_space_id);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "H5Sget_simple_extent_dims() returned "
			 "an unexpected value");
		return -1;
	}

	nstart_buf = new_IntAE(ndim, ndim, 0);
	nblock_buf = new_IntAE(ndim, ndim, 0);
	last_block_start_buf = new_LLongAE(ndim, ndim, 0);

	nstart = nstart_buf->elts;
	nblock = nblock_buf->elts;
	last_block_start = last_block_start_buf->elts;

	/* This call will populate 'nstart', 'ans_dim', 'nblock',
	   and 'last_block_start'. */
	ret = deep_check_starts_counts(ndim, dset_dims,
				starts, counts,
				nstart, ans_dim,
				nblock, last_block_start);
	if (ret < 0) {
		free(dset_dims);
		H5Sclose(file_space_id);
		return -1;
	}

	if (!noreduce && selection_can_be_reduced(ndim, nstart, nblock)) {
		reduced = PROTECT(reduce_selection(starts, counts,
						   ans_dim,
						   nblock, last_block_start));
		starts = VECTOR_ELT(reduced, 0);
		counts = VECTOR_ELT(reduced, 1);
		nstart = nblock;
	}

	ret = H5Sselect_none(file_space_id);
	if (ret < 0)
		return -1;

	//clock_t t0 = clock();
	nselection = counts != R_NilValue ?
		select_hyperslabs(file_space_id,
				ndim, dset_dims,
				starts, counts, nstart) :
		select_elements(file_space_id,
				ndim, dset_dims,
				starts, ans_dim);
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("time for setting selection: %e\n", dt);
	//printf("nselection: %lld, time per selection unit: %e\n",
	//	nselection, dt / nselection);
	if (ret < 0)
		return -1;

	if (nstart == nblock)
		UNPROTECT(1);  // unprotect 'reduced'

	free(dset_dims);
	if (ret < 0) {
		H5Sclose(file_space_id);
		return -1;
	}
	return file_space_id;
}

static hid_t prepare_mem_space(int ndim, const int *ans_dim)
{
	hsize_t *dims;
	int along, h5along, ret;
	hid_t mem_space_id;

	/* Allocate and set 'dims'. */
	dims = (hsize_t *) malloc(ndim * sizeof(hsize_t));
	if (dims == NULL) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'dims'");
		return -1;
	}
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dims[h5along] = ans_dim[along];

	mem_space_id = H5Screate_simple(ndim, dims, NULL);
	if (mem_space_id < 0) {
		free(dims);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Screate_simple() returned an error");
		return -1;
	}
	ret = H5Sselect_all(mem_space_id);
	if (ret < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Sselect_hyperslab() returned an error");
		return -1;
	}
	free(dims);
	return mem_space_id;
}

/* Return R_NilValue if error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce)
{
	int ret, ndim, along;
	SEXPTYPE ans_type;
	SEXP ans_dim, ans;
	hid_t file_space_id, mem_space_id, mem_type_id;
	R_xlen_t ans_len;
	void *buf;

	ret = get_ans_type(dset_id, &ans_type);
	if (ret < 0)
		return R_NilValue;

	ndim = shallow_check_starts_counts(starts, counts);
	if (ndim < 0)
		return R_NilValue;

	ans_dim = PROTECT(NEW_INTEGER(ndim));

	/* Prepare 'file_space_id'. */
	file_space_id = prepare_file_space(dset_id, starts, counts,
					   INTEGER(ans_dim), noreduce);
	if (file_space_id < 0) {
		UNPROTECT(1);
		return R_NilValue;
	}

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= INTEGER(ans_dim)[along];

	if (ans_len != 0) {
		/* Prepare 'mem_space_id'. */
		mem_space_id = prepare_mem_space(ndim, INTEGER(ans_dim));
		if (mem_space_id < 0) {
			H5Sclose(file_space_id);
			UNPROTECT(1);
			return R_NilValue;
		}
	}

	ans = PROTECT(allocVector(ans_type, ans_len));
	SET_DIM(ans, ans_dim);

	if (ans_len != 0) {
		if (ans_type == INTSXP) {
			buf = INTEGER(ans);
			mem_type_id = H5T_NATIVE_INT;
		} else {
			buf = REAL(ans);
			mem_type_id = H5T_NATIVE_DOUBLE;
		}
		//clock_t t0 = clock();
		ret = H5Dread(dset_id, mem_type_id,
			      mem_space_id, file_space_id,
			      H5P_DEFAULT, buf);
		//printf("time for reading data from selection: %e\n",
		//	(1.0 * clock() - t0) / CLOCKS_PER_SEC);
		H5Sclose(mem_space_id);
	}

	H5Sclose(file_space_id);
	UNPROTECT(2);
	if (ans_len != 0 && ret < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Dread() returned an error");
		return R_NilValue;
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce)
{
	SEXP filepath0, name0, ans;
	int noreduce0;
	hid_t file_id, dset_id;

	/* Check 'filepath'. */
	if (!(IS_CHARACTER(filepath) && LENGTH(filepath) == 1))
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' cannot be NA");

	/* Check 'name'. */
	if (!(IS_CHARACTER(name) && LENGTH(name) == 1))
		error("'name' must be a single string");
	name0 = STRING_ELT(name, 0);
	if (name0 == NA_STRING)
		error("'name' cannot be NA");

	/* Check 'noreduce'. */
	if (!(IS_LOGICAL(noreduce) && LENGTH(noreduce) == 1))
		error("'noreduce' must be TRUE or FALSE");
	noreduce0 = LOGICAL(noreduce)[0];

	file_id = H5Fopen(CHAR(filepath0), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		error("failed to open file %s", CHAR(filepath0));
	dset_id = H5Dopen(file_id, CHAR(name0), H5P_DEFAULT);
	if (dset_id < 0) {
		H5Fclose(file_id);
		error("failed to open dataset %s from file %s",
		      CHAR(name0), CHAR(filepath0));
	}
	ans = PROTECT(h5mread(dset_id, starts, counts, noreduce0));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ans == R_NilValue) {
		UNPROTECT(1);
		error(errmsg_buf);
	}
	UNPROTECT(1);
	return ans;
}

