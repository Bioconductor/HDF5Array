/****************************************************************************
 *                 Workhorses behind h5mread methods 4, 5, 6                *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_index.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "h5mread_helpers.h"
#include "ChunkIterator.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */
//#include <time.h>


/****************************************************************************
 * copy_selected_chunk_data_to_Rarray()
 */

static void init_in_offset_and_out_offset(int ndim, SEXP index,
			const int *out_dim, const H5Viewport *out_vp,
			const H5Viewport *h5dset_vp,
			const hsize_t *h5chunkdim,
			size_t *in_offset, size_t *out_offset)
{
	size_t in_off, out_off;
	int along, h5along, i;
	SEXP start;

	in_off = out_off = 0;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunkdim[h5along];
		out_off *= out_dim[along];
		i = out_vp->off[along];
		start = GET_LIST_ELT(index, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  h5dset_vp->h5off[h5along];
		out_off += i;
	}
	*in_offset = in_off;
	*out_offset = out_off;
	//printf("# in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static inline void update_in_offset_and_out_offset(int ndim,
		SEXP index,
		const hsize_t *h5chunkdim,
		const H5Viewport *out_vp,
		const int *inner_midx, int inner_moved_along,
		const int *out_dim,
		size_t *in_offset, size_t *out_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc, out_off_inc;

	start = GET_LIST_ELT(index, inner_moved_along);
	if (start != R_NilValue) {
		i1 = out_vp->off[inner_moved_along] +
		     inner_midx[inner_moved_along];
		i0 = i1 - 1;
		in_off_inc = _get_trusted_elt(start, i1) -
			     _get_trusted_elt(start, i0);
	} else {
		in_off_inc = 1;
	}
	out_off_inc = 1;
	if (inner_moved_along >= 1) {
		along = inner_moved_along - 1;
		h5along = ndim - inner_moved_along;
		do {
			in_off_inc *= h5chunkdim[h5along];
			out_off_inc *= out_dim[along];
			di = 1 - out_vp->dim[along];
			start = GET_LIST_ELT(index, along);
			if (start != R_NilValue) {
				i1 = out_vp->off[along];
				i0 = i1 - di;
				in_off_inc += _get_trusted_elt(start, i1) -
					      _get_trusted_elt(start, i0);
			} else {
				in_off_inc += di;
			}
			out_off_inc += di;
			along--;
			h5along++;
		} while (along >= 0);
	}
	*in_offset += in_off_inc;
	*out_offset += out_off_inc;
	//printf("## in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static inline void copy_vlen_string_to_character_Rarray(
		const H5DSetDescriptor *h5dset,
		const char * const *in, size_t in_offset,
		SEXP Rarray, size_t Rarray_offset)
{
	const char *s;
	int is_na;
	SEXP Rarray_elt;

	/* Variable length strings are always null-terminated. */
	s = in[in_offset];
	is_na = h5dset->as_na_attr && s[0] == 'N' && s[1] == 'A' && s[2] == 0;
	if (is_na) {
		SET_STRING_ELT(Rarray, Rarray_offset, NA_STRING);
	} else {
		Rarray_elt = PROTECT(mkChar(s));
		SET_STRING_ELT(Rarray, Rarray_offset, Rarray_elt);
		UNPROTECT(1);
	}
	return;
}

static inline void copy_string_to_character_Rarray(
		const H5DSetDescriptor *h5dset,
		const char *in, size_t in_offset,
		SEXP Rarray, size_t Rarray_offset)
{
	size_t h5type_size;
	const char *s;
	int s_len, is_na;
	SEXP Rarray_elt;

	h5type_size = h5dset->h5type->h5type_size;
	s = in + in_offset * h5type_size;
	for (s_len = 0; s_len < h5type_size; s_len++)
		if (s[s_len] == 0)
			break;
	is_na = h5dset->as_na_attr && s_len == 2 && s[0] == 'N' && s[1] == 'A';
	if (is_na) {
		SET_STRING_ELT(Rarray, Rarray_offset, NA_STRING);
	} else {
		Rarray_elt = PROTECT(mkCharLen(s, s_len));
		SET_STRING_ELT(Rarray, Rarray_offset, Rarray_elt);
		UNPROTECT(1);
	}
	return;
}

static long long int copy_selected_string_chunk_data_to_character_Rarray(
		const ChunkIterator *chunk_iter, int *inner_midx_buf,
		const void *in, size_t in_offset,
		const int *Rarray_dim, SEXP Rarray, size_t Rarray_offset)
{
	const H5DSetDescriptor *h5dset;
	int ndim, inner_moved_along;
	long long int nvals;

	h5dset = chunk_iter->h5dset;
	ndim = h5dset->ndim;
	nvals = 0;
	while (1) {
		if (h5dset->h5type->is_variable_str)
			copy_vlen_string_to_character_Rarray(h5dset,
						in, in_offset,
						Rarray, Rarray_offset);
		else
			copy_string_to_character_Rarray(h5dset,
						in, in_offset,
						Rarray, Rarray_offset);
		nvals++;
		inner_moved_along = _next_midx(ndim,
					       chunk_iter->mem_vp.dim,
					       inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset_and_out_offset(ndim,
				chunk_iter->index,
				h5dset->h5chunkdim,
				&chunk_iter->mem_vp,
				inner_midx_buf,
				inner_moved_along,
				Rarray_dim,
				&in_offset, &Rarray_offset);
	};
	return nvals;
}

#define	ARGS_AND_BODY_OF_COPY_FUNCTION(in_type, out_type)(		\
		const ChunkIterator *chunk_iter, int *inner_midx_buf,	\
		const in_type *in, size_t in_offset,			\
		const int *out_dim, out_type *out, size_t out_offset)	\
{									\
	const H5DSetDescriptor *h5dset;					\
	int ndim, inner_moved_along;					\
	long long int nvals;						\
									\
	h5dset = chunk_iter->h5dset;					\
	ndim = h5dset->ndim;						\
	nvals = 0;							\
	while (1) {							\
		out[out_offset] = in[in_offset];			\
		nvals++;						\
		inner_moved_along = _next_midx(ndim,			\
					chunk_iter->mem_vp.dim,		\
					inner_midx_buf);		\
		if (inner_moved_along == ndim)				\
			break;						\
		update_in_offset_and_out_offset(ndim,			\
				chunk_iter->index,			\
				h5dset->h5chunkdim,			\
				&chunk_iter->mem_vp,			\
				inner_midx_buf,				\
				inner_moved_along,			\
				out_dim,				\
				&in_offset, &out_offset);		\
	};								\
	return nvals;							\
}

/* copy_selected_XXX_chunk_data_to_int_array() functions: copy ints and
   any smaller standard native type to an array of ints. */
static long long int copy_selected_int_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(int, int)
static long long int copy_selected_char_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(char, int)
static long long int copy_selected_schar_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(signed char, int)
static long long int copy_selected_uchar_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, int)
static long long int copy_selected_short_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(short, int)
static long long int copy_selected_ushort_chunk_data_to_int_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned short, int)

/* copy_selected_XXX_chunk_data_to_double_array() functions: copy doubles
   and any smaller standard native type to an array of doubles. */
static long long int copy_selected_double_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(double, double)
static long long int copy_selected_char_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(char, double)
static long long int copy_selected_schar_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(signed char, double)
static long long int copy_selected_uchar_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, double)
static long long int copy_selected_short_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(short, double)
static long long int copy_selected_ushort_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned short, double)
static long long int copy_selected_int_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(int, double)
static long long int copy_selected_uint_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned int, double)
static long long int copy_selected_long_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(long, double)
static long long int copy_selected_ulong_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned long, double)
static long long int copy_selected_llong_chunk_data_to_double_array // be safe
	ARGS_AND_BODY_OF_COPY_FUNCTION(long long, double)
static long long int copy_selected_ullong_chunk_data_to_double_array // be safe
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned long long, double)
static long long int copy_selected_float_chunk_data_to_double_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(float, double)

/* Copy unsigned chars (a.k.a Rbytes) to an array of unsigned chars. */
static long long int copy_selected_uchar_chunk_data_to_uchar_array
	ARGS_AND_BODY_OF_COPY_FUNCTION(unsigned char, unsigned char)

static int copy_selected_chunk_data_to_Rarray(
		const ChunkIterator *chunk_iter,
		ChunkDataBuffer *chunk_data_buf,
		int *inner_midx_buf,
		const int *Rarray_dim, SEXP Rarray)
{
	const H5DSetDescriptor *h5dset;
	size_t in_offset, out_offset;
	SEXPTYPE Rtype;
	long long int nvals;
	int copy_without_type_casting;
	void *out;

	//clock_t t0;
	//double dt;

	//t0 = clock();
	h5dset = chunk_iter->h5dset;
	init_in_offset_and_out_offset(h5dset->ndim, chunk_iter->index,
			Rarray_dim, &chunk_iter->mem_vp,
			&chunk_iter->h5dset_vp, h5dset->h5chunkdim,
			&in_offset, &out_offset);
	Rtype = h5dset->h5type->Rtype;
	if (Rtype == STRSXP) {
		//printf("- copying selected chunk character data ... ");
		nvals = copy_selected_string_chunk_data_to_character_Rarray(
				chunk_iter, inner_midx_buf,
				chunk_data_buf->data, in_offset,
				Rarray_dim, Rarray, out_offset);
		//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
		//printf("ok (%lld value%s copied in %3.3f ms)\n",
		//       nvals, nvals == 1 ? "" : "s", dt);
		return 0;
	}
	copy_without_type_casting = chunk_data_buf->data_type_id ==
				    h5dset->h5type->native_type_id_for_Rtype;
	//printf("- copying selected chunk data %s type casting ... ",
	//       copy_without_type_casting ? "WITHOUT" : "WITH");
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		out = Rtype == INTSXP ? INTEGER(Rarray) : LOGICAL(Rarray);
		if (copy_without_type_casting) {
			nvals = copy_selected_int_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_CHAR) {
			nvals = copy_selected_char_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SCHAR) {
			nvals = copy_selected_schar_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UCHAR) {
			nvals = copy_selected_uchar_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SHORT) {
			nvals = copy_selected_short_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_USHORT) {
			nvals = copy_selected_ushort_chunk_data_to_int_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	    case REALSXP:
		out = REAL(Rarray);
		if (copy_without_type_casting) {
			nvals = copy_selected_double_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_CHAR) {
			nvals = copy_selected_char_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SCHAR) {
			nvals = copy_selected_schar_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UCHAR) {
			nvals = copy_selected_uchar_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_SHORT) {
			nvals = copy_selected_short_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_USHORT) {
			nvals = copy_selected_ushort_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_INT) {
			nvals = copy_selected_int_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_UINT) {
			nvals = copy_selected_uint_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_LONG) {
			nvals = copy_selected_long_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_ULONG) {
			nvals = copy_selected_ulong_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_LLONG) {
			nvals = copy_selected_llong_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_ULLONG) {
			nvals = copy_selected_ullong_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		if (chunk_data_buf->data_type_id == H5T_NATIVE_FLOAT) {
			nvals = copy_selected_float_chunk_data_to_double_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	    case RAWSXP:
		out = RAW(Rarray);
		if (copy_without_type_casting) {
			nvals = copy_selected_uchar_chunk_data_to_uchar_array(
					chunk_iter, inner_midx_buf,
					chunk_data_buf->data, in_offset,
					Rarray_dim, out, out_offset);
			break;
		}
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	    default:
		PRINT_TO_ERRMSG_BUF("unsupported dataset type");
		return -1;
	}
	if (nvals < 0)
		return -1;
	//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
	//printf("ok (%lld value%s copied in %3.3f ms)\n",
	//       nvals, nvals == 1 ? "" : "s", dt);
	return 0;
}


/****************************************************************************
 * read_data_4_5()
 *
 * method 4: One call to _read_H5Viewport() or _read_h5chunk() (wrappers for
 * H5Dread() or H5Dread_chunk(), respectively) per chunk touched by the
 * user-supplied array selection.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Make one call to _read_H5Viewport() or _read_h5chunk() to load the
 *     **entire** chunk data to an intermediate buffer.
 *   - Copy the user-selected data from the intermediate buffer to 'Rarray'.
 *
 * method 5: Like method 4 but bypasses the intermediate buffer if a
 * chunk is fully selected.
 */

static int read_data_4_5(ChunkIterator *chunk_iter,
		const int *Rarray_dim, SEXP Rarray,
		int method, int use_H5Dread_chunk)
{
	const H5DSetDescriptor *h5dset;
	int ndim, ret, direct_load;
	IntAE *inner_midx_buf;
	void *out;
	hid_t out_space_id;
	ChunkDataBuffer chunk_data_buf;

	if (use_H5Dread_chunk)
		warning("using 'use.H5Dread_chunk=TRUE' is still "
			"experimental, use at your own risk");

	h5dset = chunk_iter->h5dset;
	ndim = h5dset->ndim;
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	out = DATAPTR(Rarray);
	if (out == NULL)
		return -1;

	out_space_id = _create_mem_space(ndim, Rarray_dim);
	if (out_space_id < 0)
		return -1;

	ret = _init_ChunkDataBuffer(&chunk_data_buf, h5dset, 0);
	if (ret < 0) {
		H5Sclose(out_space_id);
		return ret;
	}
	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		//_print_tchunk_info(chunk_iter);
		direct_load = method == 5 && _tchunk_is_fully_selected(ndim,
							&chunk_iter->h5dset_vp,
							&chunk_iter->mem_vp);
		if (direct_load) {
			/* Load the chunk **directly** into 'Rarray' (no
			   intermediate buffer). */
			ret = _read_H5Viewport(h5dset,
				&chunk_iter->h5dset_vp,
				h5dset->h5type->native_type_id_for_Rtype,
				out_space_id, out,
				&chunk_iter->mem_vp);
		} else {
			/* Load the **entire** chunk to an intermediate
			   buffer then copy the user-selected chunk data
			   from the intermediate buffer to 'Rarray'. */
			//clock_t t0 = clock();
			ret = _load_chunk(chunk_iter,
					&chunk_data_buf,
					use_H5Dread_chunk);
			if (ret < 0)
				break;
			//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
			//printf("- load chunk: %3.3f ms\n", dt);

			ret = copy_selected_chunk_data_to_Rarray(
					chunk_iter,
					&chunk_data_buf,
					inner_midx_buf->elts,
					Rarray_dim, Rarray);
			if (ret < 0)
				break;
			if (h5dset->h5type->is_variable_str)
				_reclaim_vlen_bufs(&chunk_data_buf);
		}
		if (ret < 0)
			break;
	}
	_destroy_ChunkDataBuffer(&chunk_data_buf);
	H5Sclose(out_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to _read_h5selection() (wrapper for H5Dread()) per chunk touched
 * by the user-supplied array selection. No intermediate buffer.
 *
 * More precisely, walk over the chunks touched by 'index'. For each chunk:
 *   - Select the hyperslabs obtained by intersecting the user-supplied
 *     array selection with the current chunk.
 *   - Call _read_h5selection(). This loads the selected data **directly**
 *     to the final R array.
 */

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP index,
		const H5Viewport *out_vp,
		IntAEAE *inner_breakpoint_bufs, int *inner_nchip_buf)
{
	int along, d, off, nchip, i;
	IntAE *inner_breakpoint_buf;
	SEXP start;
	long long int s0, s1;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		inner_breakpoint_buf = inner_breakpoint_bufs->elts[along];
		IntAE_set_nelt(inner_breakpoint_buf, 0);
		d = out_vp->dim[along];
		start = GET_LIST_ELT(index, along);
		if (start == R_NilValue) {
			IntAE_insert_at(inner_breakpoint_buf, 0, d);
			inner_nchip_buf[along] = 1;
			continue;
		}
		off = out_vp->off[along];
		s1 = _get_trusted_elt(start, off);
		nchip = 0;
		for (i = 1; i < d; i++) {
			s0 = s1;
			s1 = _get_trusted_elt(start, off + i);
			if (s1 != s0 + 1)
				IntAE_insert_at(inner_breakpoint_buf,
						nchip++, i);
		}
		IntAE_insert_at(inner_breakpoint_buf, nchip++, d);
		inner_nchip_buf[along] = nchip;
	}
	return;
}

static void init_inner_vp(int ndim, SEXP index,
		const H5Viewport *h5dset_vp,
		H5Viewport *inner_vp)
{
	int along, h5along;
	SEXP start;
	hsize_t d;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(index, along);
		if (start == R_NilValue) {
			inner_vp->h5off[h5along] =
					h5dset_vp->h5off[h5along];
			d = h5dset_vp->h5dim[h5along];
		} else {
			d = 1;
		}
		inner_vp->h5dim[h5along] = d;
	}
	return;
}

static void update_inner_vp(int ndim,
		SEXP index, const H5Viewport *out_vp,
		const int *inner_midx, int inner_moved_along,
		const IntAEAE *inner_breakpoint_bufs,
		H5Viewport *inner_vp)
{
	int along, h5along, idx, off, d, i;
	SEXP start;
	const int *inner_breakpoint;

	for (along = 0; along < ndim; along++) {
		if (along > inner_moved_along)
			break;
		if (index == R_NilValue)
			continue;
		start = VECTOR_ELT(index, along);
		if (start == R_NilValue)
			continue;
		inner_breakpoint = inner_breakpoint_bufs->elts[along]->elts;
		idx = inner_midx[along];
		off = idx == 0 ? 0 : inner_breakpoint[idx - 1];
		d = inner_breakpoint[idx] - off;
		i = out_vp->off[along] + off;
		h5along = ndim - 1 - along;
		inner_vp->h5off[h5along] = _get_trusted_elt(start, i) - 1;
		inner_vp->h5dim[h5along] = d;
	}
	return;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_intersection_of_chips_with_chunk(
		const H5DSetDescriptor *h5dset, SEXP index,
		const H5Viewport *out_vp, const H5Viewport *h5dset_vp,
		const IntAEAE *inner_breakpoint_bufs,
		const int *inner_nchip,
		int *inner_midx_buf,
		H5Viewport *inner_vp)
{
	int ret, ndim, inner_moved_along;
	long long int num_hyperslabs;

	ret = H5Sselect_none(h5dset->h5space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = h5dset->ndim;

	init_inner_vp(ndim, index, h5dset_vp, inner_vp);

	/* Walk on the "inner chips" i.e. on the intersections between
	   the "chips" in the user-supplied array selection and the currrent
	   chunk. */
	num_hyperslabs = 0;
	inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_inner_vp(ndim, index, out_vp,
				inner_midx_buf, inner_moved_along,
				inner_breakpoint_bufs,
				inner_vp);
		ret = _add_H5Viewport_to_h5selection(h5dset->h5space_id,
						     inner_vp);
		if (ret < 0)
			return -1;
		inner_moved_along = _next_midx(ndim, inner_nchip,
					       inner_midx_buf);
	} while (inner_moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(inner_nchip)

	return num_hyperslabs;
}

static int direct_load_selected_chunk_data(
		const ChunkIterator *chunk_iter,
		int *inner_midx_buf,
		H5Viewport *inner_vp,
		IntAEAE *inner_breakpoint_bufs,
		const IntAE *inner_nchip_buf,
		hid_t out_space_id, void *out)
{
	const H5DSetDescriptor *h5dset;
	int ret;

	h5dset = chunk_iter->h5dset;
	update_inner_breakpoints(h5dset->ndim, chunk_iter->moved_along,
			chunk_iter->index, &chunk_iter->mem_vp,
			inner_breakpoint_bufs, inner_nchip_buf->elts);
	ret = select_intersection_of_chips_with_chunk(
			h5dset, chunk_iter->index,
			&chunk_iter->mem_vp, &chunk_iter->h5dset_vp,
			inner_breakpoint_bufs, inner_nchip_buf->elts,
			inner_midx_buf, inner_vp);
	if (ret < 0)
		return ret;
	ret = _read_h5selection(h5dset,
				h5dset->h5type->native_type_id_for_Rtype,
				out_space_id, out,
				&chunk_iter->mem_vp);
	return ret;
}

static int read_data_6(ChunkIterator *chunk_iter,
		const int *Rarray_dim, SEXP Rarray)
{
	const H5DSetDescriptor *h5dset;
	int ndim, ret;
	IntAE *inner_midx_buf, *inner_nchip_buf;
	IntAEAE *inner_breakpoint_bufs;
	void *out;
	hid_t out_space_id;
	H5Viewport inner_vp;

	h5dset = chunk_iter->h5dset;
	ndim = h5dset->ndim;
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_nchip_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	out = DATAPTR(Rarray);
	if (out == NULL)
		return -1;

	out_space_id = _create_mem_space(ndim, Rarray_dim);
	if (out_space_id < 0)
		return -1;

	ret = _alloc_H5Viewport(&inner_vp, ndim, ALLOC_H5OFF_AND_H5DIM);
	if (ret < 0) {
		H5Sclose(out_space_id);
		return ret;
	}

	/* Walk over the chunks touched by the user-supplied array selection. */
	while ((ret = _next_chunk(chunk_iter))) {
		if (ret < 0)
			break;
		ret = direct_load_selected_chunk_data(
			chunk_iter,
			inner_midx_buf->elts,
			&inner_vp,
			inner_breakpoint_bufs, inner_nchip_buf,
			out_space_id, out);
		if (ret < 0)
			break;
	}
	_free_H5Viewport(&inner_vp);
	H5Sclose(out_space_id);
	return ret;
}


/****************************************************************************
 * _h5mread_index()
 *
 * Implements methods 4 to 6.
 * Return an ordinary array or R_NilValue if an error occured.
 */

SEXP _h5mread_index(const H5DSetDescriptor *h5dset, SEXP index,
		    int method, int use_H5Dread_chunk, int *ans_dim)
{
	ChunkIterator chunk_iter;
	int ret, ndim, along;
	R_xlen_t ans_len;
	SEXP ans;

	/* In the context of methods 5 & 6, 'chunk_iter.mem_vp.h5off'
	   and 'chunk_iter.mem_vp.h5dim' will be used, not just
	   'chunk_iter.mem_vp.off' and 'chunk_iter.mem_vp.dim',
	   so we set 'alloc_full_mem_vp' (last arg) to 1. */
	ret = _init_ChunkIterator(&chunk_iter, h5dset, index,
				  ans_dim, method == 5 || method == 6);
	if (ret < 0)
		return R_NilValue;

	ndim = h5dset->ndim;
	for (along = 0, ans_len = 1; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(h5dset->h5type->Rtype, ans_len));

	if (method <= 5) {
		/* methods 4 and 5 */
		ret = read_data_4_5(&chunk_iter, ans_dim, ans,
				    method, use_H5Dread_chunk);
	} else {
		/* method 6 */
		ret = read_data_6(&chunk_iter, ans_dim, ans);
	}
	_destroy_ChunkIterator(&chunk_iter);
	UNPROTECT(1);
	return ret < 0 ? R_NilValue : ans;
}

