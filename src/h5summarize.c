/****************************************************************************
 *                     Summarization of an HDF5 dataset                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5summarize.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "h5mread_helpers.h"
#include "ChunkIterator.h"

#include <string.h>  /* for strcmp */
#include <limits.h>  /* for INT_MAX */
//#include <time.h>

#define	MIN_OPCODE	1
#define	MAX_OPCODE	2
#define	RANGE_OPCODE	3
#define	SUM_OPCODE	4
#define	PROD_OPCODE	5
#define	ANY_OPCODE	6
#define	ALL_OPCODE	7

static int get_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "min") == 0)
		return MIN_OPCODE;
	if (strcmp(s, "max") == 0)
		return MAX_OPCODE;
	if (strcmp(s, "range") == 0)
		return RANGE_OPCODE;
	if (strcmp(s, "sum") == 0)
		return SUM_OPCODE;
	if (strcmp(s, "prod") == 0)
		return PROD_OPCODE;
	if (strcmp(s, "any") == 0)
		return ANY_OPCODE;
	if (strcmp(s, "all") == 0)
		return ALL_OPCODE;
	error("'op' must be one of: \"min\", \"max\", \"range\", "
	      "\"sum\", \"prod\", \"any\", \"all\"");
	return 0;
}

static int check_Rtype(SEXPTYPE Rtype, int opcode)
{
	switch (Rtype) {
	    case LGLSXP: return 0;
	    case INTSXP: case REALSXP:
		if (opcode != ANY_OPCODE && opcode != ALL_OPCODE)
			return 0;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return -1;
}


/****************************************************************************
 * Callback functions used for summarization of int data.
 *
 * Return new status:
 *   0 = 'init' has not been set yet
 *   1 = 'init' has been set
 *   2 = 'init' has been set and we don't need to continue (break condition)
 */

typedef int (*IntOP)(void *init, int x, int na_rm, int status);
typedef int (*DoubleOP)(void *init, double x, int na_rm, int status);

static inline int int_min(void *init, int x, int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		*int_init = x;
		return 2;
	}
	if (status == 0 || x < *int_init)
		*int_init = x;
	return 1;
}
static inline int double_min(void *init, double x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (R_IsNA(x) || R_IsNaN(x)) {
		if (na_rm)
			return status;
		*double_init = x;
		return R_IsNA(x) ? 2 : 1;
	}
	if (!R_IsNaN(*double_init) && x < *double_init)
		*double_init = x;
	return 1;
}

static inline int int_max(void *init, int x, int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		*int_init = x;
		return 2;
	}
	if (status == 0 || x > *int_init)
		*int_init = x;
	return 1;
}
static inline int double_max(void *init, double x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (R_IsNA(x) || R_IsNaN(x)) {
		if (na_rm)
			return status;
		*double_init = x;
		return R_IsNA(x) ? 2 : 1;
	}
	if (!R_IsNaN(*double_init) && x > *double_init)
		*double_init = x;
	return 1;
}

static inline int int_range(void *init, int x, int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		int_init[0] = int_init[1] = x;
		return 2;
	}
	if (status == 0) {
		int_init[0] = int_init[1] = x;
		return 1;
	}
	if (x < int_init[0])
		int_init[0] = x;
	if (x > int_init[1])
		int_init[1] = x;
	return 1;
}
static inline int double_range(void *init, double x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (R_IsNA(x) || R_IsNaN(x)) {
		if (na_rm)
			return status;
		double_init[0] = double_init[1] = x;
		return R_IsNA(x) ? 2 : 1;
	}
	if (!R_IsNaN(double_init[0])) {
		if (x < double_init[0])
			double_init[0] = x;
		if (x > double_init[1])
			double_init[1] = x;
	}
	return 1;
}

static inline int int_sum(void *init, int x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		*double_init = NA_REAL;
		return 2;
	}
	*double_init += x;
	return 1;
}
static inline int double_sum(void *init, double x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (R_IsNA(x) || R_IsNaN(x)) {
		if (na_rm)
			return status;
		*double_init = x;
		return R_IsNA(x) ? 2 : 1;
	}
	if (!R_IsNaN(*double_init))
		*double_init += x;
	return 1;
}

static inline int int_prod(void *init, int x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		*double_init = NA_REAL;
		return 2;
	}
	*double_init *= x;
	return 1;
}
static inline int double_prod(void *init, double x, int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	if (R_IsNA(x) || R_IsNaN(x)) {
		if (na_rm)
			return status;
		*double_init = x;
		return R_IsNA(x) ? 2 : 1;
	}
	if (!R_IsNaN(*double_init))
		*double_init *= x;
	return 1;
}

static inline int int_any(void *init, int x, int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		*int_init = x;
		return 1;
	}
	if (x == 0)
		return 1;
	*int_init = x;
	return 2;
}

static inline int int_all(void *init, int x, int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	if (x == NA_INTEGER) {
		if (na_rm)
			return status;
		*int_init = x;
		return 1;
	}
	if (x != 0)
		return 1;
	*int_init = x;
	return 2;
}

/* One of '*int_OP' or '*double_OP' will be set to NULL and the other one
   to a non-NULL value. */
static void select_OP(int opcode, SEXPTYPE Rtype,
		IntOP *int_OP, DoubleOP *double_OP, void *init)
{
	int *int_init = (int *) init;
	double *double_init = (double *) init;

	*int_OP = NULL;
	*double_OP = NULL;
	if (opcode == ANY_OPCODE) {
		*int_OP = int_any;
		int_init[0] = 0;
		return;
	}
	if (opcode == ALL_OPCODE) {
		*int_OP = int_all;
		int_init[0] = 1;
		return;
	}
	if (opcode == SUM_OPCODE) {
		if (Rtype == REALSXP) {
			*double_OP = double_sum;
		} else {
			*int_OP = int_sum;
		}
		double_init[0] = 0.0;
		return;
	}
	if (opcode == PROD_OPCODE) {
		if (Rtype == REALSXP) {
			*double_OP = double_prod;
		} else {
			*int_OP = int_prod;
		}
		double_init[0] = 1.0;
		return;
	}
	if (Rtype == REALSXP) {
		switch (opcode) {
		    case MIN_OPCODE:
			*double_OP = double_min;
			double_init[0] = R_PosInf;
			break;
		    case MAX_OPCODE:
			*double_OP = double_max;
			double_init[0] = R_NegInf;
			break;
		    case RANGE_OPCODE:
			*double_OP = double_range;
			double_init[0] = R_PosInf;
			double_init[1] = R_NegInf;
			break;
		}
		return;
	}
	/* NO initial value! */
	switch (opcode) {
	    case MIN_OPCODE:   *int_OP = int_min;   break;
	    case MAX_OPCODE:   *int_OP = int_max;   break;
	    case RANGE_OPCODE: *int_OP = int_range; break;
	}
	return;
}

static SEXP init2SEXP(int opcode, SEXPTYPE Rtype, void *init, int status)
{
	int *int_init = (int *) init;
	double *double_init = (double *) init;
	SEXP ans;

	if (opcode == ANY_OPCODE || opcode == ALL_OPCODE)
		return ScalarLogical(int_init[0]);
	if (opcode == MIN_OPCODE || opcode == MAX_OPCODE) {
		if (Rtype == REALSXP)
			return ScalarReal(double_init[0]);
		if (status == 0) {
			return ScalarReal(opcode == MIN_OPCODE ? R_PosInf
							       : R_NegInf);
		} else {
			return ScalarInteger(int_init[0]);
		}
	}
	if (opcode == RANGE_OPCODE) {
		if (Rtype == REALSXP) {
			ans = PROTECT(NEW_NUMERIC(2));
			REAL(ans)[0] = double_init[0];
			REAL(ans)[1] = double_init[1];
		} else {
			if (status == 0) {
				ans = PROTECT(NEW_NUMERIC(2));
				REAL(ans)[0] = R_PosInf;
				REAL(ans)[1] = R_NegInf;
			} else {
				ans = PROTECT(NEW_INTEGER(2));
				INTEGER(ans)[0] = int_init[0];
				INTEGER(ans)[1] = int_init[1];
			}
		}
		UNPROTECT(1);
		return ans;
	}
	/* 'opcode' is either SUM_OPCODE or PROD_OPCODE. */
	if (Rtype == REALSXP)
		return ScalarReal(double_init[0]);
	/* Direct comparison with NA_REAL is safe. No need to use R_IsNA(). */
	if (double_init[0] == NA_REAL)
		return ScalarInteger(NA_INTEGER);
	if (double_init[0] <= INT_MAX && double_init[0] >= -INT_MAX)
		return ScalarInteger((int) double_init[0]);
	return ScalarReal(double_init[0]);
}


/****************************************************************************
 * Summarize chunk data
 */

/* Does NOT work properly on a truncated chunk! Works properly only if the
   chunk data fills the full 'chunk_data_buf.data', that is, if the current
   chunk is a full-size chunk and not a "truncated" chunk (a.k.a. "partial
   edge chunk" in HDF5's terminology). */
static void summarize_full_chunk_int_data(
		const int *data, size_t data_length,
		IntOP int_OP, void *init, int na_rm, int *status)
{
	size_t offset;

	/* Walk on the **all** elements in current chunk. */
	for (offset = 0; offset < data_length; offset++) {
		*status = int_OP(init, data[offset], na_rm, *status);
		if (*status == 2)
			return;
	}
	return;
}

/* Does NOT work properly on a truncated chunk! Works properly only if the
   chunk data fills the full 'chunk_data_buf.data', that is, if the current
   chunk is a full-size chunk and not a "truncated" chunk (a.k.a. "partial
   edge chunk" in HDF5's terminology). */
static void summarize_full_chunk_double_data(
		const double *data, size_t data_length,
		DoubleOP double_OP, void *init, int na_rm, int *status)
{
	size_t offset;

	/* Walk on the **all** elements in current chunk. */
	for (offset = 0; offset < data_length; offset++) {
		*status = double_OP(init, data[offset], na_rm, *status);
		if (*status == 2)
			return;
	};
	return;
}

static void summarize_selected_chunk_int_data(
		const H5DSetDescriptor *h5dset, SEXP index,
		const int *data, const H5Viewport *h5dset_vp,
		const H5Viewport *mem_vp, int *inner_midx_buf,
		IntOP int_OP, void *init, int na_rm, int *status)
{
	int ndim, inner_moved_along;
	size_t offset;

	ndim = h5dset->ndim;
	_init_in_offset(ndim, index, h5dset->h5chunkdim,
			mem_vp, h5dset_vp,
			&offset);
	/* Walk on the **selected** elements in current chunk. */
	while (1) {
		*status = int_OP(init, data[offset], na_rm, *status);
		if (*status == 2)
			break;
		inner_moved_along = _next_midx(ndim, mem_vp->dim,
					       inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		_update_in_offset(ndim, index, h5dset->h5chunkdim, mem_vp,
				  inner_midx_buf, inner_moved_along,
				  &offset);
	};
	return;
}

static void summarize_selected_chunk_double_data(
		const H5DSetDescriptor *h5dset, SEXP index,
		const double *data, const H5Viewport *h5dset_vp,
		const H5Viewport *mem_vp, int *inner_midx_buf,
		DoubleOP double_OP, void *init, int na_rm, int *status)
{
	int ndim, inner_moved_along;
	size_t offset;

	ndim = h5dset->ndim;
	_init_in_offset(ndim, index, h5dset->h5chunkdim,
			mem_vp, h5dset_vp,
			&offset);
	/* Walk on the **selected** elements in current chunk. */
	while (1) {
		*status = double_OP(init, data[offset], na_rm, *status);
		if (*status == 2)
			break;
		inner_moved_along = _next_midx(ndim, mem_vp->dim,
					       inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		_update_in_offset(ndim, index, h5dset->h5chunkdim, mem_vp,
				  inner_midx_buf, inner_moved_along,
				  &offset);
	};
	return;
}

static SEXP h5summarize(const H5DSetDescriptor *h5dset, SEXP index,
		int opcode, int na_rm, int verbose)
{
	IntOP int_OP;
	DoubleOP double_OP;
	double init[2];  /* 'init' will store 1 or 2 ints or doubles */
	//double init2[2];
	int ndim, status, ret, go_fast;
	IntAE *inner_midx_buf;
	ChunkIterator chunk_iter;
	ChunkDataBuffer chunk_data_buf;
	//ChunkDataBuffer chunk_data_buf2;

	/* Set one of 'int_OP' or 'double_OP' to NULL and the other one to
	   an IntOP or DoubleOP function. Also initializes 'init' with 1 or
	   2 ints or doubles. */
	select_OP(opcode, h5dset->h5type->Rtype, &int_OP, &double_OP, init);
	//select_OP(opcode, h5dset->h5type->Rtype, &int_OP, &double_OP, init2);

	ndim = h5dset->ndim;
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	status = 0;

	/* In the context of h5summarize(), we won't use
	   'chunk_iter.mem_vp.h5off' or 'chunk_iter.mem_vp.h5dim', only
	   'chunk_iter.mem_vp.off' and 'chunk_iter.mem_vp.dim', so we
	   set 'alloc_full_mem_vp' (last arg) to 0. */
	ret = _init_ChunkIterator(&chunk_iter, h5dset, index, NULL, 0);
	if (ret < 0)
		return R_NilValue;

	ret = _init_ChunkDataBuffer(&chunk_data_buf, h5dset, 1);
	if (ret < 0) {
		_destroy_ChunkIterator(&chunk_iter);
		return R_NilValue;
	}

	//ret = _init_ChunkDataBuffer(&chunk_data_buf2, h5dset, 0);
	//if (ret < 0) {
	//	_destroy_ChunkIterator(&chunk_iter);
	//	return R_NilValue;
	//}

	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(&chunk_iter))) {
		if (ret < 0)
			break;
		if (verbose)
			_print_tchunk_info(&chunk_iter);

		//clock_t t0 = clock();
		ret = _load_chunk(&chunk_iter, &chunk_data_buf, 0);
		if (ret < 0)
			break;
		//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
		//printf("_load_chunk() to chunk_data_buf: %2.3f ms\n", dt);

		//t0 = clock();
		//ret = _load_chunk(&chunk_iter, &chunk_data_buf2, 0);
		//if (ret < 0)
		//	break;
		//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
		//printf("_load_chunk() to chunk_data_buf2: %2.3f ms\n", dt);

		go_fast = _tchunk_is_fully_selected(ndim,
				&chunk_iter.h5dset_vp,
				&chunk_iter.mem_vp) &&
			  ! _tchunk_is_truncated(h5dset,
				&chunk_iter.h5dset_vp);
		if (go_fast) {
			if (int_OP != NULL) {
				//t0 = clock();
				summarize_full_chunk_int_data(
					chunk_data_buf.data,
					chunk_data_buf.data_length,
					int_OP, init, na_rm, &status);
				//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
				//printf("summarize_full_chunk_int_data: %2.3f ms\n", dt);
				//t0 = clock();
				//summarize_full_chunk_int16_data(
				//	chunk_data_buf2.data,
				//	chunk_data_buf2.data_length,
				//	int_OP, init2, na_rm, &status);
				//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
				//printf("summarize_full_chunk_int16_data: %2.3f ms\n", dt);
			} else {
				summarize_full_chunk_double_data(
					chunk_data_buf.data,
					chunk_data_buf.data_length,
					double_OP, init, na_rm, &status);
			}
		} else {
			if (int_OP != NULL) {
				summarize_selected_chunk_int_data(
					h5dset, index,
					chunk_data_buf.data,
					&chunk_iter.h5dset_vp,
					&chunk_iter.mem_vp,
					inner_midx_buf->elts,
					int_OP, init, na_rm, &status);
			} else {
				summarize_selected_chunk_double_data(
					h5dset, index,
					chunk_data_buf.data,
					&chunk_iter.h5dset_vp,
					&chunk_iter.mem_vp,
					inner_midx_buf->elts,
					double_OP, init, na_rm, &status);
			}
		}
		if (status == 2)
			break;
	}
	//_destroy_ChunkDataBuffer(&chunk_data_buf2);
	_destroy_ChunkDataBuffer(&chunk_data_buf);
	_destroy_ChunkIterator(&chunk_iter);
	if (ret < 0)
		return R_NilValue;
	return init2SEXP(opcode, h5dset->h5type->Rtype, init, status);
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5summarize(SEXP filepath, SEXP name, SEXP index, SEXP as_integer,
		   SEXP op, SEXP na_rm, SEXP verbose)
{
	int as_int, opcode, narm0, verbose0, ret, is_supported;
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;
	const H5TypeDescriptor *h5type;
	SEXP ans;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	/* Check 'op'. */
	opcode = get_opcode(op);

	/* Check 'na_rm'. */
	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm0 = LOGICAL(na_rm)[0];

	/* Check 'verbose'. */
	if (!(IS_LOGICAL(verbose) && LENGTH(verbose) == 1))
		error("'verbose' must be TRUE or FALSE");
	verbose0 = LOGICAL(verbose)[0];

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	ret = _init_H5DSetDescriptor(&h5dset, dset_id, as_int, 0);
	if (ret < 0) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_global_errmsg_buf());
	}

	h5type = h5dset.h5type;
	is_supported = h5type->Rtype_is_set && !h5type->is_variable_str;
	if (!is_supported) {
		_destroy_H5DSetDescriptor(&h5dset);
		H5Dclose(dset_id);
		H5Fclose(file_id);
		PRINT_TO_ERRMSG_BUF(
			"h5summarize() does not support this type "
			"of dataset yet, sorry. You can\n  "
			"use 'H5DSetDescriptor(filepath, name)' "
			"to see details about the dataset.");
		error(_HDF5Array_global_errmsg_buf());
	}

	ret = check_Rtype(h5type->Rtype, opcode);
	if (ret < 0) {
		_destroy_H5DSetDescriptor(&h5dset);
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_global_errmsg_buf());
	}

	ans = PROTECT(h5summarize(&h5dset, index,
				  opcode, narm0, verbose0));

	_destroy_H5DSetDescriptor(&h5dset);
	H5Dclose(dset_id);
	H5Fclose(file_id);
	UNPROTECT(1);
	if (ans == R_NilValue)
		error(_HDF5Array_global_errmsg_buf());
	return ans;
}

