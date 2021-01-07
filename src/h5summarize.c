/****************************************************************************
 *                     Summarization of an HDF5 dataset                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5summarize.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "H5DSetDescriptor.h"
#include "h5mread_helpers.h"

#include "hdf5.h"
#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for strcmp */

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

static SEXP summarize0(int opcode, SEXPTYPE Rtype)
{
	SEXP ans;

	switch (opcode) {
	    case MIN_OPCODE:
		ans = PROTECT(NEW_NUMERIC(1));
		REAL(ans)[0] = R_PosInf;
		break;
	    case MAX_OPCODE:
		ans = PROTECT(NEW_NUMERIC(1));
		REAL(ans)[0] = R_NegInf;
		break;
	    case RANGE_OPCODE:
		ans = PROTECT(NEW_NUMERIC(2));
		REAL(ans)[0] = R_PosInf;
		REAL(ans)[1] = R_NegInf;
		break;
	    case SUM_OPCODE:
		ans = PROTECT(NEW_NUMERIC(1));
		REAL(ans)[0] = 0.0;
		break;
	    case PROD_OPCODE:
		ans = PROTECT(NEW_NUMERIC(1));
		REAL(ans)[0] = 1.0;
		break;
	    case ANY_OPCODE:
		ans = PROTECT(NEW_LOGICAL(1));
		LOGICAL(ans)[0] = 0;
		break;
	    case ALL_OPCODE:
		ans = PROTECT(NEW_LOGICAL(1));
		LOGICAL(ans)[0] = 1;
		break;
	    default:
		/* Should never happen. */
		PRINT_TO_ERRMSG_BUF("invalid opcode: %d", opcode);
		return R_NilValue;
	}
	UNPROTECT(1);
	return ans;
}

/* Return non-zero value if a break condition was encountered. */
static inline int summarize_val(
		const H5DSetDescriptor *h5dset,
		const int *in, size_t in_offset,
		int opcode, int na_rm, SEXP ans, int *done)
{
	int x;
	double xx;

	switch (h5dset->Rtype) {
	    case LGLSXP:
		x = ((int *) in)[in_offset];
		/* See fix_logical_NAs() in h5mread.c for why we replace
		   negative values with NA_LOGICAL. */
		if (x < 0)
			x = NA_LOGICAL;
		if (opcode == ANY_OPCODE) {
			if (x == NA_LOGICAL) {
				if (!na_rm) {
					LOGICAL(ans)[0] = NA_LOGICAL;
					*done = 1;
				}
			} else if (x != 0) {
				LOGICAL(ans)[0] = 1;
				*done = 1;
			}
			return 0;
		}
		if (opcode == ALL_OPCODE) {
			if (x == NA_LOGICAL) {
				if (!na_rm) {
					LOGICAL(ans)[0] = NA_LOGICAL;
					*done = 1;
				}
			} else if (x == 0) {
				LOGICAL(ans)[0] = 0;
				*done = 1;
			}
			return 0;
		}
		xx = x == NA_LOGICAL ? NA_REAL : (double) x;
		break;
	    case INTSXP:
		x = ((int *) in)[in_offset];
		xx = x == NA_INTEGER ? NA_REAL : (double) x;
		break;
	    case REALSXP:
		xx = ((double *) in)[in_offset];
		break;
	    default:
		/* Should never happen. Early call to check_Rtype() should
		   have caught this. */
		PRINT_TO_ERRMSG_BUF("unsupported type: %s",
				    CHAR(type2str(h5dset->Rtype)));
		return -1;
	}
	if (R_IsNA(xx) || R_IsNaN(xx)) {
		if (!na_rm) {
			if (opcode == RANGE_OPCODE) {
				REAL(ans)[0] = REAL(ans)[1] = NA_REAL;
			} else {
				REAL(ans)[0] = NA_REAL;
			}
			*done = 1;
		}
		return 0;
	}
	switch (opcode) {
	    case MIN_OPCODE:
		if (xx < REAL(ans)[0]) {
			REAL(ans)[0] = xx;
			if (na_rm && xx == R_NegInf)
				*done = 1;
		}
		break;
	    case MAX_OPCODE:
		if (xx > REAL(ans)[0]) {
			REAL(ans)[0] = xx;
			if (na_rm && xx == R_PosInf)
				*done = 1;
		}
		break;
	    case RANGE_OPCODE:
		if (xx < REAL(ans)[0])
			REAL(ans)[0] = xx;
		if (xx > REAL(ans)[1])
			REAL(ans)[1] = xx;
		if (na_rm
		 && REAL(ans)[0] == R_NegInf
		 && REAL(ans)[1] == R_PosInf)
			*done = 1;
		break;
	    case SUM_OPCODE:
		REAL(ans)[0] += xx;
		break;
	    case PROD_OPCODE:
		REAL(ans)[0] *= xx;
		break;
	    default:
		/* Should never happen. */
		PRINT_TO_ERRMSG_BUF("invalid opcode: %d", opcode);
		return -1;
	}
	return 0;
}

static int summarize_selected_chunk_data(
		const H5DSetDescriptor *h5dset,
		SEXP index, const void *in, const H5Viewport *tchunk_vp,
		const H5Viewport *dest_vp, int *inner_midx_buf,
		int opcode, int na_rm, SEXP ans)
{
	int ndim, done, inner_moved_along, ret;
	size_t in_offset;

	ndim = h5dset->ndim;
	_init_in_offset(ndim, index, h5dset->h5chunkdim, dest_vp,
			tchunk_vp,
			&in_offset);
	done = 0;
	/* Walk on the **selected** elements in current chunk. */
	while (1) {
		ret = summarize_val(h5dset, in, in_offset,
				    opcode, na_rm, ans, &done);
		if (ret < 0)
			return -1;
		if (done)
			break;
		inner_moved_along = _next_midx(ndim, dest_vp->dim,
					       inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		_update_in_offset(ndim, index, h5dset->h5chunkdim, dest_vp,
				  inner_midx_buf, inner_moved_along,
				  &in_offset);
	};
	return 0;
}

static SEXP h5summarize(const H5DSetDescriptor *h5dset,
		SEXP index,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		const int *num_tchunks,
		int opcode, int na_rm, int verbose)
{
	int ndim, moved_along, ret;
	IntAE *tchunk_midx_buf, *inner_midx_buf;
	void *chunk_data_buf;
	hid_t chunk_space_id;
	H5Viewport tchunk_vp, middle_vp, dest_vp;
	long long int tchunk_rank;
	SEXP ans;

	ndim = h5dset->ndim;

	/* Prepare buffers. */

	tchunk_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	chunk_data_buf = malloc(h5dset->chunk_data_buf_size);
	if (chunk_data_buf == NULL) {
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return R_NilValue;
	}
	chunk_space_id = H5Screate_simple(ndim, h5dset->h5chunkdim, NULL);
	if (chunk_space_id < 0) {
		free(chunk_data_buf);
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return R_NilValue;
	}

	/* Allocate 'tchunk_vp', 'middle_vp', and 'dest_vp'.
	   We set 'dest_vp_mode' to ALLOC_OFF_AND_DIM because, in the
	   context of h5summarize(), we won't use 'dest_vp.h5off' or
	   'dest_vp.h5dim', only 'dest_vp.off' and 'dest_vp.dim'. */
	if (_alloc_tchunk_vp_middle_vp_dest_vp(ndim,
	        &tchunk_vp, &middle_vp, &dest_vp,
	        ALLOC_OFF_AND_DIM) < 0)
	{
		H5Sclose(chunk_space_id);
		free(chunk_data_buf);
		return R_NilValue;
	}

	ans = PROTECT(summarize0(opcode, h5dset->Rtype));
	if (ans == R_NilValue)
		goto done;

	/* Walk over the chunks touched by the user-supplied array selection. */

	tchunk_rank = 0;
	moved_along = ndim;
	do {
		_update_tchunk_vp_dest_vp(h5dset,
				tchunk_midx_buf->elts, moved_along,
				index, breakpoint_bufs, tchunkidx_bufs,
				&tchunk_vp, &dest_vp);
		if (verbose)
			_print_tchunk_info(ndim,
					num_tchunks, tchunk_midx_buf->elts,
					tchunk_rank,
					index,
					tchunkidx_bufs, &tchunk_vp);
		ret = _read_H5Viewport(h5dset,
				&tchunk_vp, &middle_vp,
				chunk_data_buf, chunk_space_id);
		if (ret < 0) {
			ans = R_NilValue;
			break;
		}
		ret = summarize_selected_chunk_data(h5dset,
				index, chunk_data_buf, &tchunk_vp,
				&dest_vp, inner_midx_buf->elts,
				opcode, na_rm, ans);
		if (ret < 0) {
			ans = R_NilValue;
			break;
		}
		tchunk_rank++;
		moved_along = _next_midx(ndim, num_tchunks,
					 tchunk_midx_buf->elts);
	} while (moved_along < ndim);

    done:
	_free_tchunk_vp_middle_vp_dest_vp(&tchunk_vp, &middle_vp, &dest_vp);
	H5Sclose(chunk_space_id);
	free(chunk_data_buf);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5summarize(SEXP filepath, SEXP name, SEXP index, SEXP as_integer,
		   SEXP op, SEXP na_rm, SEXP verbose)
{
	int as_int, opcode, narm0, verbose0, ret;
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;  /* touched chunk ids along each dim */
	IntAE *ntchunk_buf;  /* nb of touched chunks along each dim */
	long long int total_num_tchunks;
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
	ret = check_Rtype(h5dset.Rtype, opcode);
	if (ret < 0) {
		_destroy_H5DSetDescriptor(&h5dset);
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_global_errmsg_buf());
	}

	/* This call will populate 'breakpoint_bufs' and 'tchunkidx_bufs'. */
	breakpoint_bufs = new_IntAEAE(h5dset.ndim, h5dset.ndim);
	tchunkidx_bufs = new_LLongAEAE(h5dset.ndim, h5dset.ndim);
	ret = _map_starts_to_h5chunks(&h5dset, index, NULL,
				      breakpoint_bufs, tchunkidx_bufs);
	if (ret < 0) {
		_destroy_H5DSetDescriptor(&h5dset);
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_global_errmsg_buf());
	}
	ntchunk_buf = new_IntAE(h5dset.ndim, h5dset.ndim, 0);
	total_num_tchunks = _set_num_tchunks(&h5dset, index,
					     tchunkidx_bufs, ntchunk_buf->elts);

	ans = PROTECT(h5summarize(&h5dset, index,
				  breakpoint_bufs, tchunkidx_bufs,
				  ntchunk_buf->elts,
				  opcode, narm0, verbose0));

	_destroy_H5DSetDescriptor(&h5dset);
	H5Dclose(dset_id);
	H5Fclose(file_id);
	UNPROTECT(1);
	if (ans == R_NilValue)
		error(_HDF5Array_global_errmsg_buf());
	return ans;
}

