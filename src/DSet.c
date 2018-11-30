/****************************************************************************
 *                    Basic manipulation of a DSet struct                   *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include <stdlib.h>  /* for malloc, free */
#include <limits.h>  /* for INT_MAX */


hsize_t *_alloc_hsize_t_buf(size_t buflength, int zeroes, const char *what)
{
	hsize_t *buf;
	int i;

	buf = (hsize_t *) malloc(buflength * sizeof(hsize_t));
	if (buf == NULL)
		PRINT_TO_ERRMSG_BUF("failed to allocate memory for %s", what);
	if (zeroes) {
		for (i = 0; i < buflength; i++)
			buf[i] = 0;
	}
	return buf;
}

/* See hdf5-1.10.3/src/H5Tpublic.h for the list of datatype classes. We only
   support H5T_INTEGER, H5T_FLOAT, and H5T_STRING for now. */
static int get_ans_type_from_class(H5T_class_t class, int as_int, size_t size,
				   SEXPTYPE *ans_type)
{
	const char *classname;

	switch (class) {
	    case H5T_INTEGER:
		*ans_type = (as_int || size <= sizeof(int)) ? INTSXP : REALSXP;
		return 0;
	    case H5T_FLOAT:
		*ans_type = as_int ? INTSXP : REALSXP;
		return 0;
	    case H5T_STRING:
		if (as_int)
			warning("'as.integer' is ignored when "
				"dataset class is H5T_STRING");
		*ans_type = STRSXP;
		return 0;
	    case H5T_TIME: classname = "H5T_TIME"; break;
	    case H5T_BITFIELD: classname = "H5T_BITFIELD"; break;
	    case H5T_OPAQUE: classname = "H5T_OPAQUE"; break;
	    case H5T_COMPOUND: classname = "H5T_COMPOUND"; break;
	    case H5T_REFERENCE: classname = "H5T_REFERENCE"; break;
	    case H5T_ENUM: classname = "H5T_ENUM"; break;
	    case H5T_VLEN: classname = "H5T_VLEN"; break;
	    case H5T_ARRAY: classname = "H5T_ARRAY"; break;
	    default:
		PRINT_TO_ERRMSG_BUF("unknown dataset class identifier: %d",
				    class);
	    return -1;
	}
	PRINT_TO_ERRMSG_BUF("unsupported dataset class: %s", classname);
	return -1;
}

static size_t get_ans_elt_size_from_ans_type(SEXPTYPE ans_type, size_t size)
{
	switch (ans_type) {
	    case INTSXP:  return sizeof(int);
	    case REALSXP: return sizeof(double);
	    case STRSXP:  return size;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(ans_type)));
	return 0;
}

static hid_t get_mem_type_id_from_ans_type(SEXPTYPE ans_type, hid_t dtype_id)
{
	switch (ans_type) {
	    case INTSXP:  return H5T_NATIVE_INT;
	    case REALSXP: return H5T_NATIVE_DOUBLE;
	    case STRSXP:  return dtype_id;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(ans_type)));
	return -1;
}

void _close_DSet(DSet *dset)
{
	if (dset->h5nchunk != NULL)
		free(dset->h5nchunk);
	if (dset->h5dim != NULL)
		free(dset->h5dim);
	if (dset->plist_id != -1)
		H5Pclose(dset->plist_id);
	if (dset->space_id != -1)
		H5Sclose(dset->space_id);
	if (dset->dtype_id != -1)
		H5Tclose(dset->dtype_id);
	H5Dclose(dset->dset_id);
}

int _get_DSet(hid_t dset_id, int ndim, int as_int, int ans_type_only,
	      DSet *dset)
{
	hid_t space_id, plist_id, dtype_id, mem_type_id;
	int dset_ndim, *h5nchunk, h5along;
	hsize_t *h5dim, *h5chunk_spacings, d, spacing, nchunk;
	H5T_class_t class;
	size_t size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE ans_type;

	dset->dset_id = dset_id;
	dset->dtype_id = -1;
	dset->space_id = -1;
	dset->plist_id = -1;
	dset->h5dim = NULL;
	dset->h5nchunk = NULL;

	/* Set 'dtype_id'. */
	dtype_id = H5Dget_type(dset_id);
	if (dtype_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		return -1;
	}
	dset->dtype_id = dtype_id;

	/* Set 'class'. */
	class = H5Tget_class(dtype_id);
	if (class == H5T_NO_CLASS) {
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		return -1;
	}
	dset->class = class;

	/* Set 'size'. */
	size = H5Tget_size(dtype_id);
	if (size == 0) {
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF("H5Tget_size() returned 0");
		return -1;
	}
	dset->size = size;

	/* Set 'ans_type'. */
	if (get_ans_type_from_class(class, as_int, size, &ans_type) < 0) {
		H5Tclose(dtype_id);
		return -1;
	}
	if (ans_type == STRSXP && H5Tis_variable_str(dtype_id) != 0) {
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF("reading variable-length string data "
				    "is not supported at the moment");
		return -1;
	}
	dset->ans_type = ans_type;

	if (ans_type_only) {
		H5Tclose(dtype_id);
		return 0;
	}

	/* Set 'space_id'. */
	space_id = H5Dget_space(dset_id);
	if (space_id < 0) {
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		return -1;
	}
	dset->space_id = space_id;

	/* Set 'ndim'. */
	dset_ndim = H5Sget_simple_extent_ndims(space_id);
	if (dset_ndim < 0) {
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF(
			"H5Sget_simple_extent_ndims() returned an error");
		return -1;
	}
	if (ndim != dset_ndim) {
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF(
			"Dataset has %d dimensions but 'starts' has %d list "
			"element%s.\n  'starts' must have one list element "
			"per dimension in the dataset.",
			dset_ndim, ndim, ndim > 1 ? "s" : "");
		return -1;
	}
	dset->ndim = ndim;

	/* Set 'plist_id'. */
	plist_id = H5Dget_create_plist(dset_id);
	if (plist_id < 0) {
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		return -1;
	}
	dset->plist_id = plist_id;

	/* Allocate 'h5dim', 'h5chunk_spacings', and 'h5nchunk'. */
	h5dim = _alloc_hsize_t_buf(2 * ndim, 0,
				   "'h5dim' and 'h5chunk_spacings'");
	if (h5dim == NULL) {
		H5Pclose(plist_id);
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		return -1;
	}
	h5chunk_spacings = h5dim + ndim;
	h5nchunk = (int *) malloc(ndim * sizeof(int));
	if (h5nchunk == NULL) {
		free(h5dim);
		H5Pclose(plist_id);
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'h5nchunk'");
		return -1;
	}

	/* Set 'h5dim'. */
	if (H5Sget_simple_extent_dims(space_id, h5dim, NULL) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		goto on_error;
	}
	dset->h5dim = h5dim;

	/* Set 'h5chunk_spacings'. */
	if (H5Pget_chunk(plist_id, ndim, h5chunk_spacings) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
				    "an unexpected value");
		goto on_error;
	}
	dset->h5chunk_spacings = h5chunk_spacings;

	/* Set 'h5nchunk'. */
	for (h5along = 0; h5along < ndim; h5along++) {
		d = h5dim[h5along];
		if (d == 0) {
			h5nchunk[h5along] = 0;
			continue;
		}
		spacing = h5chunk_spacings[h5along];
		nchunk = d / spacing;
		if (d % spacing != 0)
			nchunk++;
		if (nchunk > INT_MAX) {
			PRINT_TO_ERRMSG_BUF("datasets with more than "
				"INT_MAX chunks along any dimension\n  "
				"are not supported at the moment");
			goto on_error;
		}
		h5nchunk[h5along] = nchunk;
	}
	dset->h5nchunk = h5nchunk;

	/* Set 'ans_elt_size'. */
	ans_elt_size = get_ans_elt_size_from_ans_type(ans_type, size);
	if (ans_elt_size == 0)
		goto on_error;
	dset->ans_elt_size = ans_elt_size;

	/* Set 'chunk_data_buf_size'. */
	chunk_data_buf_size = ans_elt_size;
	for (h5along = 0; h5along < ndim; h5along++)
		chunk_data_buf_size *= h5chunk_spacings[h5along];
	dset->chunk_data_buf_size = chunk_data_buf_size;

	/* Set 'mem_type_id'. */
	mem_type_id = get_mem_type_id_from_ans_type(ans_type, dtype_id);
	if (mem_type_id < 0)
		goto on_error;
	dset->mem_type_id = mem_type_id;
	return 0;

    on_error:
	_close_DSet(dset);
	return -1;
}

