/****************************************************************************
 *            Exploring alternate rhdf5::h5read() implementations           *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread.h"

#include "global_errmsg_buf.h"
#include "array_selection.h"
#include "H5DSetDescriptor.h"
#include "h5mread_startscounts.h"
#include "h5mread_starts.h"
#include "h5mread_sparse.h"

#include "hdf5.h"

/* Return -1 on error. */
static int select_method(const H5DSetDescriptor *h5dset,
			 SEXP starts, SEXP counts, int sparse, int method)
{
	int along;

	if (sparse) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "'as.sparse' is set to TRUE");
			return -1;
		}
		if (h5dset->h5chunkdim == NULL) {
			PRINT_TO_ERRMSG_BUF("'as.sparse=TRUE' is not supported "
					    "on a contiguous dataset");
			return -1;
		}
		if (method == 0) {
			method = 8;
		} else if (method != 8) {
			PRINT_TO_ERRMSG_BUF("only method 8 is supported "
					    "when 'as.sparse' is set to TRUE");
			return -1;
		}
		return method;
	}
	if (h5dset->Rtype == STRSXP) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "reading string data");
			return -1;
		}
		if (method == 0) {
			method = 4;
		} else if (method != 4 && method != 5) {
			PRINT_TO_ERRMSG_BUF("only methods 4 and 5 are "
					    "supported when reading "
					    "string data");
			return -1;
		}
		return method;
	}
	if (method == 0) {
		method = 1;
		/* March 27, 2019: My early testing (from Nov 2018) seemed
		   to indicate that method 6 was a better choice over method 4
		   when the layout is chunked and 'counts' is NULL. Turns out
		   that doing more testing today seems to indicate the opposite
		   i.e. method 4 now seems to perform better than method 6 on
		   all the datasets I've tested so far, including those used
		   by Pete Hickey here:
		     https://github.com/Bioconductor/DelayedArray/issues/13
		   and those used in the examples in man/h5mread.Rd.
		   Note sure what happened between Nov 2018 and today. Did I
		   do something stupid in my early testing? Did something
		   change in Rhdf5lib?
		   Anyway thanks to Pete for providing such a useful report.

		   Nov 26, 2019: I added method 7. Is like method 4 but
		   bypasses the intermediate buffer if a chunk is fully
		   selected. This is now preferred over methods 4 or 6. */
		if (h5dset->h5chunkdim != NULL &&
		    counts == R_NilValue &&
		    starts != R_NilValue)
		{
			for (along = 0; along < h5dset->ndim; along++) {
				if (VECTOR_ELT(starts, along) != R_NilValue) {
					//method = 6;
					//method = 4;
					method = 7;
					break;
				}
			}
		}
	} else if (method < 0 || method > 7) {
		PRINT_TO_ERRMSG_BUF("'method' must be >= 0 and <= 7");
		return -1;
	} else if (method >= 4) {
		/* Make sure the data is chunked and 'counts' is NULL. */
		if (h5dset->h5chunkdim == NULL) {
			PRINT_TO_ERRMSG_BUF("methods 4, 5, 6, and 7 cannot "
				"be used on a contiguous dataset (unless\n  "
				"it contains string data in which case "
				"methods 4 and 5 can be used)");
			return -1;
		}
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("methods 4, 5, 6, and 7 can "
				"only be used when 'counts' is NULL");
			return -1;
		}
	}
	return method;
}

/* If the H5 datatype that was used to store the logical data is an 8-bit
   or 16-bit signed integer type (e.g. H5T_STD_I8LE or H5T_STD_I16BE) then
   NA values got loaded as negative values that are not equal to NA_LOGICAL
   (e.g. as -128 for H5T_STD_I8LE and -2^16 for H5T_STD_I16BE).
   These values must be replaced with NA_LOGICAL. */
static void fix_logical_NAs(SEXP x)
{
	R_xlen_t x_len, i;
	int *x_p;

	x_len = XLENGTH(x);
	for (i = 0, x_p = LOGICAL(x); i < x_len; i++, x_p++) {
		if (*x_p < 0)
			*x_p = NA_LOGICAL;
	}
	return;
}

/* Replace "NA" strings in 'x' with character NAs (NA_character_).
   'x' is assumed to be NA-free. */
static void set_character_NAs(SEXP x)
{
	R_xlen_t x_len, i;
	int x_elt_len;
	SEXP x_elt;

	x_len = XLENGTH(x);
	for (i = 0; i < x_len; i++) {
		x_elt = STRING_ELT(x, i);
		x_elt_len = LENGTH(x_elt);
		if (x_elt_len == 2) {
			const char *s = CHAR(x_elt);
			if (s[0] == 'N' && s[1] == 'A')
				SET_STRING_ELT(x, i, NA_STRING);
		}
	}
	return;
}

/* Return R_NilValue on error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce,
		    int as_int, int sparse, int method)
{
	SEXP ans, ans_dim;
	H5DSetDescriptor h5dset;
	int ret;

	ans = R_NilValue;

	if (_init_H5DSetDescriptor(&h5dset, dset_id, as_int, 0) < 0)
		return ans;

	ret = _shallow_check_selection(h5dset.ndim, starts, counts);
	if (ret < 0)
		goto on_error;

	method = select_method(&h5dset, starts, counts, sparse, method);
	if (method < 0)
		goto on_error;

	ans_dim = PROTECT(NEW_INTEGER(h5dset.ndim));

	if (method <= 3) {
		/* Implements methods 1 to 3. */
		ans = _h5mread_startscounts(&h5dset, starts, counts, noreduce,
					    method, INTEGER(ans_dim));
	} else if (method <= 7) {
		/* Implements methods 4 to 7. */
		ans = _h5mread_starts(&h5dset, starts,
				      method, INTEGER(ans_dim));
	} else {
		/* Implements method 8.
		   Return 'list(nzindex, nzdata, NULL)' or R_NilValue if
		   an error occured. */
		ans = _h5mread_sparse(&h5dset, starts, INTEGER(ans_dim));
	}

	if (ans != R_NilValue) {
		PROTECT(ans);
		if (sparse) {
			if (h5dset.Rtype == LGLSXP)
				fix_logical_NAs(VECTOR_ELT(ans, 1));
			else if (h5dset.Rtype == STRSXP && h5dset.as_na_attr)
				set_character_NAs(VECTOR_ELT(ans, 1));
			/* Final 'ans' is 'list(nzindex, nzdata, ans_dim)'. */
			SET_VECTOR_ELT(ans, 2, ans_dim);
		} else {
			if (h5dset.Rtype == LGLSXP)
				fix_logical_NAs(ans);
			SET_DIM(ans, ans_dim);
		}
		UNPROTECT(1);  /* 'ans' */
	}

	UNPROTECT(1);  /* 'ans_dim' */

    on_error:
	_destroy_H5DSetDescriptor(&h5dset);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce,
	       SEXP as_integer, SEXP as_sparse, SEXP method)
{
	int noreduce0, as_int, sparse, method0;
	hid_t file_id, dset_id;
	SEXP ans;

	/* Check 'noreduce'. */
	if (!(IS_LOGICAL(noreduce) && LENGTH(noreduce) == 1))
		error("'noreduce' must be TRUE or FALSE");
	noreduce0 = LOGICAL(noreduce)[0];

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	/* Check 'as_sparse'. */
	if (!(IS_LOGICAL(as_sparse) && LENGTH(as_sparse) == 1))
		error("'as_sparse' must be TRUE or FALSE");
	sparse = LOGICAL(as_sparse)[0];

	/* Check 'method'. */
	if (!(IS_INTEGER(method) && LENGTH(method) == 1))
		error("'method' must be a single integer");
	method0 = INTEGER(method)[0];

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	ans = PROTECT(h5mread(dset_id, starts, counts, noreduce0,
			      as_int, sparse, method0));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	UNPROTECT(1);
	if (ans == R_NilValue)
		error(_HDF5Array_global_errmsg_buf());
	return ans;
}

