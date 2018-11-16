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


static char errmsg_buf[256];

#define	NOT_A_FINITE_NUMBER(x) \
	(R_IsNA(x) || R_IsNaN(x) || (x) == R_PosInf || (x) == R_NegInf)

static int get_untrusted_elt(SEXP x, int i, long long int *val,
			     const char *what, int along)
{
	int tmp1;
	double tmp2;

	if (IS_INTEGER(x)) {
		tmp1 = INTEGER(x)[i];
		if (tmp1 == NA_INTEGER) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "%s[[%d]][%d] is NA",
				 what, along + 1, i + 1);
			return -1;
		}
		*val = (long long int) tmp1;
	} else {
		tmp2 = REAL(x)[i];
		if (NOT_A_FINITE_NUMBER(tmp2)) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "%s[[%d]][%d] is NA or NaN "
				 "or not a finite number",
				 what, along + 1, i + 1);
			return -1;
		}
		if (tmp2 > (double) LLONG_MAX ||
		    tmp2 < (double) LLONG_MIN) {
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
	int ndim, counts_len;

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
		counts_len = LENGTH(counts);
		if (counts_len != ndim) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'starts' and 'counts' must have "
				 "the same length");
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
							    SEXP count)
{
	const char *msg = "selection must be strictly ascending "
			  "along each dimension, but\n  you have:";
	if (count != R_NilValue)
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] < starts[[%d]][%d] + "
			 "counts[[%d]][%d]",
			 msg, along1, i + 1, along1, i, along1, i);
	else
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] <= starts[[%d]][%d]",
			 msg, along1, i + 1, along1, i);
	return;
}

static void set_errmsg_for_selection_beyond_dim(int along1, int i, SEXP count)
{
	const char *msg = "selection must be within extent of "
			  "dataset, but you\n  have:";
	if (count != R_NilValue)
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] + counts[[%d]][%d] - 1 "
			 "> dimension %d in dataset",
			 msg, along1, i + 1, along1, i + 1, along1);
	else
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "%s starts[[%d]][%d] "
			 "> dimension %d in dataset",
			 msg, along1, i + 1, along1);
	return;
}

static int check_starts_counts_along(int along, long long int d,
			SEXP starts, SEXP counts,
			int *nstart, int *count_sum,
			int *nblock, long long int *last_block_start)
{
	SEXP start, count;
	int n, i, ret;
	long long int cs, c, e, s;

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
		if (nblock != NULL)
			nblock[along] = d != 0;
		if (last_block_start != NULL)
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
				 "'counts[[%d]]' must have the "
				 "same length as 'starts[[%d]]'",
				 along + 1, along + 1);
			return -1;
		}
	}
	nstart[along] = n;
	cs = 0;
	if (nblock != NULL)
		nblock[along] = 0;
	c = 1;
	e = 0;
	for (i = 0; i < n; i++) {
		ret = get_untrusted_elt(start, i, &s, "starts", along);
		if (ret < 0)
			return -1;
		if (s <= e) {
			if (e == 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "starts[[%d]][%d] is <= 0",
					 along + 1, i + 1);
				return -1;
			}
			set_errmsg_for_non_strictly_ascending_selection(
				along + 1, i, count);
			return -1;
		}
		if (i == 0 || s != e + 1) {
			if (nblock != NULL)
				nblock[along]++;
			if (last_block_start != NULL)
				last_block_start[along] = s;
		}
		if (count != R_NilValue) {
			ret = get_untrusted_elt(count, i, &c, "counts", along);
			if (ret < 0)
				return -1;
			if (c <= 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "counts[[%d]][%d] is <= 0",
					 along + 1, i + 1);
				return -1;
			}
		}
		e = s + c - 1;  // could overflow! (FIXME)
		if (d >= 0 && e > d) {
			set_errmsg_for_selection_beyond_dim(
				along + 1, i, count);
			return -1;
		}
		cs += c;  // could overflow! (FIXME)
		if (cs > INT_MAX) {
			set_error_for_selection_too_large(along + 1);
			return -1;
		}
	}
	count_sum[along] = cs;
	return 0;
}


/****************************************************************************
 * C_reduce_selection()
 */

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
	for (i = 0; i < n; i++) {
		s = get_trusted_elt(start_in, i);
		c = count_in != R_NilValue ? get_trusted_elt(count_in, i) : 1;
		if (i == 0 || s != e + 1) {
			j++;
			set_trusted_elt(start_out, j, s);
			count_out[j] = c;
		} else {
			count_out[j] += c;
		}
		e = s + c - 1;
	}
	return;
}

static void reduce_selection_along(int along,
				   const int *count_sum,
				   const int *nblock,
				   const long long int *last_block_start,
				   SEXP starts, SEXP counts,
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

/*
 * Return a list of length 2.
 * The 1st list element is the list of reduced starts and the 2nd list element
 * the list of reduced 'counts'. The 2 lists have the same shape i.e. same
 * length() and same lengths().
 */
static SEXP reduce_selection(SEXP starts, SEXP counts, int *status)
{
	int ndim, along, ret;
	IntAE *nstart_buf, *count_sum_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	long long unsigned int total_reduction;
	SEXP ans, reduced_starts, reduced_counts;

	ndim = shallow_check_starts_counts(starts, counts);
	if (ndim < 0) {
		*status = -1;
		return R_NilValue;
	}
	nstart_buf = new_IntAE(ndim, ndim, 0);
	count_sum_buf = new_IntAE(ndim, ndim, 0);
	nblock_buf = new_IntAE(ndim, ndim, 0);
	last_block_start_buf = new_LLongAE(ndim, ndim, 0);

	/* 1st pass */
	total_reduction = 0;
	for (along = 0; along < ndim; along++) {
		ret = check_starts_counts_along(along, -1,
				starts, counts,
				nstart_buf->elts, count_sum_buf->elts,
				nblock_buf->elts, last_block_start_buf->elts);
		if (ret < 0) {
			*status = -1;
			return R_NilValue;
		}
		total_reduction += nstart_buf->elts[along] -
				   nblock_buf->elts[along];
	}

	if (total_reduction == 0) {
		*status = 0;
		return R_NilValue;
	}

	/* 2nd pass */
	ans = PROTECT(NEW_LIST(2));
	reduced_starts = PROTECT(NEW_LIST(ndim));
	SET_VECTOR_ELT(ans, 0, reduced_starts);
	UNPROTECT(1);
	reduced_counts = PROTECT(NEW_LIST(ndim));
	SET_VECTOR_ELT(ans, 1, reduced_counts);
	UNPROTECT(1);
	for (along = 0; along < ndim; along++) {
		reduce_selection_along(along,
			count_sum_buf->elts,
			nblock_buf->elts, last_block_start_buf->elts,
			starts, counts,
			reduced_starts, reduced_counts);
	}
	UNPROTECT(1);
	*status = 1;
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Return NULL if the selection could not be reduced.
 */
SEXP C_reduce_selection(SEXP starts, SEXP counts, SEXP dim)
{
	SEXP ans;
	int status;

	if (dim != R_NilValue)
		error("'dim' not supported yet, sorry!");
	/* No need to PROTECT()! */
	ans = reduce_selection(starts, counts, &status);
	if (status < 0)
		error(errmsg_buf);
	return ans;
}


/****************************************************************************
 * C_h5mread()
 */

static int deep_check_starts_counts(int dset_rank, const hsize_t *dset_dims,
				    SEXP starts, SEXP counts,
				    int *nstart, int *count_sum)
{
	int along, ret;
	long long int d;

	/* We already know that 'starts' is a list. */
	if (LENGTH(starts) != dset_rank) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'starts' must have one list element "
			 "per dimension in the dataset");
		return -1;
	}
	/* We already know that 'counts' is a list (or NULL). */
	if (counts != R_NilValue) {
		if (LENGTH(counts) != dset_rank) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'counts' must have one list element "
				 "per dimension in the dataset");
			return -1;
		}
	}
	for (along = 0; along < dset_rank; along++) {
		d = dset_dims[dset_rank - 1 - along];
		ret = check_starts_counts_along(
				along, d,
				starts, counts,
				nstart, count_sum,
				NULL, NULL);
		if (ret < 0)
			return -1;
	}
	return 0;
}

/* Should we use H5Sselect_hyperslab() or H5Sselect_elements() for this?
   Useful links:
   - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
       https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
   - Documentation of H5Dread():
       https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
   - A useful example:
       https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
*/

static int add_block_to_read(hid_t file_space_id, int dset_rank,
			     SEXP starts, SEXP counts, const int *blockidx,
			     hsize_t *offset_buf, hsize_t *count_buf)
{
	int along, h5along, i, ret;
	SEXP start, count;

	/* Set 'offset_buf' and 'count_buf'. */
	for (along = 0; along < dset_rank; along++) {
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = dset_rank - 1 - along;
		i = blockidx[along];
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

static int next_block(int dset_rank, const int *nstart, int *blockidx)
{
	int along, i;

	for (along = 0; along < dset_rank; along++) {
		i = blockidx[along] + 1;
		if (i < nstart[along]) {
			blockidx[along] = i;
			return 1;
		}
		blockidx[along] = 0;
	}
	return 0;
}

static int set_region_to_read(hid_t file_space_id,
			      int dset_rank, const hsize_t *dset_dims,
			      SEXP starts, SEXP counts, const int *nstart)
{
	int ret, along, h5along;
	hsize_t *offset_buf, *count_buf;  // hyperslab offsets and dims
	IntAE *blockidx_buf;
	long long unsigned n;

	ret = H5Sselect_none(file_space_id);
	if (ret < 0)
		return -1;
	for (along = 0; along < dset_rank; along++)
		if (nstart[along] == 0)
			return 0;  // empty region (no block)

	/* Allocate 'offset_buf' and 'count_buf'. */
	offset_buf = (hsize_t *) malloc(2 * dset_rank * sizeof(hsize_t));
	if (offset_buf == NULL) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'offset_buf' "
			 "and 'count_buf'");
		return -1;
	}
	count_buf = offset_buf + dset_rank;

	/* Initialize 'offset_buf' and 'count_buf'. */
	for (along = 0, h5along = dset_rank - 1;
	     along < dset_rank;
	     along++, h5along--)
	{
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			offset_buf[h5along] = 0;
			count_buf[h5along] = dset_dims[h5along];
		} else {
			count_buf[h5along] = 1;
		}
	}

	/* Add one block at a time. */
	blockidx_buf = new_IntAE(dset_rank, dset_rank, 0);
	n = 0;
	do {
		n++;
		//printf("- indices for block %llu:", n);
		//for (along = 0; along < dset_rank; along++)
		//	printf(" %d", blockidx_buf->elts[along]);
		//printf("\n");
		ret = add_block_to_read(file_space_id, dset_rank,
					starts, counts, blockidx_buf->elts,
					offset_buf, count_buf);
		if (ret < 0)
			break;
	} while (next_block(dset_rank, nstart, blockidx_buf->elts));

	printf("nb of hyperslabs = %llu\n", n);
	free(offset_buf);
	return ret;
}

static hid_t prepare_file_space(hid_t dset_id, SEXP starts, SEXP counts,
				int *count_sum)
{
	hid_t file_space_id;
	int dset_rank, ret;
	hsize_t *dset_dims;
	IntAE *nstart_buf;

	file_space_id = H5Dget_space(dset_id);
	if (file_space_id < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Dget_space() returned an error");
		return -1;
	}

	dset_rank = H5Sget_simple_extent_ndims(file_space_id);

	/* Allocate and set 'dset_dims'. */
	dset_dims = (hsize_t *) malloc(dset_rank * sizeof(hsize_t));
	if (dset_dims == NULL) {
		H5Sclose(file_space_id);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'dset_dims'");
		return -1;
	}
	if (H5Sget_simple_extent_dims(file_space_id, dset_dims, NULL) !=
	    dset_rank)
	{
		free(dset_dims);
		H5Sclose(file_space_id);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "H5Sget_simple_extent_dims() returned "
			 "an unexpected value");
		return -1;
	}

	nstart_buf = new_IntAE(dset_rank, dset_rank, 0);

	/* This call will populate 'nstart' and 'count_sum'. */
	ret = deep_check_starts_counts(dset_rank, dset_dims,
				       starts, counts,
				       nstart_buf->elts, count_sum);
	if (ret < 0) {
		free(dset_dims);
		H5Sclose(file_space_id);
		return -1;
	}

	ret = set_region_to_read(file_space_id,
				 dset_rank, dset_dims,
				 starts, counts, nstart_buf->elts);
	free(dset_dims);
	if (ret < 0) {
		H5Sclose(file_space_id);
		return -1;
	}
	return file_space_id;
}

static hid_t prepare_mem_space(int dset_rank, const int *count_sum)
{
	hsize_t *dims;
	int along, h5along, ret;
	hid_t mem_space_id;

	/* Allocate and set 'dims'. */
	dims = (hsize_t *) malloc(dset_rank * sizeof(hsize_t));
	if (dims == NULL) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'dims'");
		return -1;
	}
	for (along = 0, h5along = dset_rank - 1;
	     along < dset_rank;
	     along++, h5along--)
	{
		dims[h5along] = count_sum[along];
	}

	mem_space_id = H5Screate_simple(dset_rank, dims, NULL);
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
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts)
{
	int ndim, along, ret;
	IntAE *count_sum_buf;
	hid_t file_space_id, mem_space_id;
	R_xlen_t ans_len;
	SEXP ans, ans_dim;

	ndim = shallow_check_starts_counts(starts, counts);
	if (ndim < 0)
		return R_NilValue;

	count_sum_buf = new_IntAE(ndim, ndim, 0);

	/* Prepare 'file_space_id'. */
	file_space_id = prepare_file_space(dset_id, starts, counts,
					   count_sum_buf->elts);
	if (file_space_id < 0)
		return R_NilValue;

	ans_dim = PROTECT(NEW_INTEGER(ndim));
	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= INTEGER(ans_dim)[along] = count_sum_buf->elts[along];

	if (ans_len != 0) {
		/* Prepare 'mem_space_id'. */
		mem_space_id = prepare_mem_space(ndim, count_sum_buf->elts);
		if (mem_space_id < 0) {
			H5Sclose(file_space_id);
			return R_NilValue;
		}
	}

	ans = PROTECT(NEW_INTEGER(ans_len));
	SET_DIM(ans, ans_dim);

	if (ans_len != 0) {
		ret = H5Dread(dset_id, H5T_NATIVE_INT,
			      mem_space_id, file_space_id,
			      H5P_DEFAULT, INTEGER(ans));
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
SEXP C_h5mread(SEXP filepath, SEXP name, SEXP starts, SEXP counts)
{
	SEXP filepath0, name0, ans;
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

	file_id = H5Fopen(CHAR(filepath0), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		error("failed to open file %s", CHAR(filepath0));
	dset_id = H5Dopen(file_id, CHAR(name0), H5P_DEFAULT);
	if (dset_id < 0) {
		H5Fclose(file_id);
		error("failed to open dataset %s from file %s",
		      CHAR(name0), CHAR(filepath0));
	}
	ans = PROTECT(h5mread(dset_id, starts, counts));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ans == R_NilValue) {
		UNPROTECT(1);
		error(errmsg_buf);
	}
	UNPROTECT(1);
	return ans;
}

