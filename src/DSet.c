/****************************************************************************
 *                    Basic manipulation of a DSet struct                   *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for strcmp */
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


/****************************************************************************
 * _get_DSet() / _close_DSet()
 */

static char *get_storage_mode_attr(hid_t dset_id)
{
	hid_t attr_id, attr_type_id;
	H5T_class_t attr_H5class;
	hsize_t attr_size;
	char *storage_mode_attr;
	herr_t ret;

	attr_id = H5Aopen(dset_id, "storage.mode", H5P_DEFAULT);
	if (attr_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aopen() returned an error");
		return NULL;
	}
	attr_type_id = H5Aget_type(attr_id);
	if (attr_type_id < 0) {
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Aget_type() returned 0");
		return NULL;
	}
	attr_H5class = H5Tget_class(attr_type_id);
	if (attr_H5class == H5T_NO_CLASS) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		return NULL;
	}
	if (attr_H5class != H5T_STRING) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("attribute \"storage.mode\" is "
				    "not of expected type H5T_STRING");
		return NULL;
	}
	attr_size = H5Aget_storage_size(attr_id);
	if (attr_size == 0) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Aget_storage_size() returned 0");
		return NULL;
	}
	//printf("attr_size = %llu\n", attr_size);
	storage_mode_attr = (char *) malloc((size_t) attr_size);
	if (storage_mode_attr == NULL) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'storage_mode_attr'");
		return NULL;
	}
	ret = H5Aread(attr_id, attr_type_id, storage_mode_attr);
	H5Tclose(attr_type_id);
	H5Aclose(attr_id);
	if (ret < 0) {
		free(storage_mode_attr);
		PRINT_TO_ERRMSG_BUF("H5Aread() returned an error");
		return NULL;
	}
	//printf("storage_mode_attr: %s\n", storage_mode_attr);
	return storage_mode_attr;
}

static int map_storage_mode_to_Rtype(const char *storage_mode, int as_int,
				     SEXPTYPE *Rtype)
{
	if (strcmp(storage_mode, "logical") == 0) {
		*Rtype = as_int ? INTSXP : LGLSXP;
		return 0;
	}
	if (strcmp(storage_mode, "integer") == 0) {
		if (as_int)
			warning("'as.integer' is ignored when the dataset to "
				"read has a \"storage.mode\" attribute set "
				"to \"integer\"");
		*Rtype = INTSXP;
		return 0;
	}
	if (strcmp(storage_mode, "double") == 0 ||
	    strcmp(storage_mode, "numeric") == 0) {
		*Rtype = as_int ? INTSXP : REALSXP;
		return 0;
	}
	if (strcmp(storage_mode, "character") == 0) {
		if (as_int)
			warning("'as.integer' is ignored when the dataset to "
				"read has a \"storage.mode\" attribute set "
				"to \"character\"");
		*Rtype = STRSXP;
		return 0;
	}
	PRINT_TO_ERRMSG_BUF("the dataset to read has a \"storage.mode\" "
			    "attribute set to unsupported value: \"%s\"",
			    storage_mode);
	return -1;
}

/* See hdf5-1.10.3/src/H5Tpublic.h for the list of datatype classes. We only
   support H5T_INTEGER, H5T_FLOAT, and H5T_STRING for now. */
static int map_H5class_to_Rtype(H5T_class_t H5class, int as_int, size_t size,
				SEXPTYPE *Rtype)
{
	const char *classname;

	switch (H5class) {
	    case H5T_INTEGER:
		*Rtype = (as_int || size <= sizeof(int)) ? INTSXP : REALSXP;
		return 0;
	    case H5T_FLOAT:
		*Rtype = as_int ? INTSXP : REALSXP;
		return 0;
	    case H5T_STRING:
		if (as_int)
			warning("'as.integer' is ignored when "
				"dataset class is H5T_STRING");
		*Rtype = STRSXP;
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
				    H5class);
	    return -1;
	}
	PRINT_TO_ERRMSG_BUF("unsupported dataset class: %s", classname);
	return -1;
}

static size_t get_ans_elt_size_from_Rtype(SEXPTYPE Rtype, size_t size)
{
	switch (Rtype) {
	    case LGLSXP:
	    case INTSXP:  return sizeof(int);
	    case REALSXP: return sizeof(double);
	    case STRSXP:  return size;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return 0;
}

static hid_t get_mem_type_id_from_Rtype(SEXPTYPE Rtype, hid_t dtype_id)
{
	switch (Rtype) {
	    case LGLSXP:
	    case INTSXP:  return H5T_NATIVE_INT;
	    case REALSXP: return H5T_NATIVE_DOUBLE;
	    case STRSXP:  return dtype_id;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return -1;
}

void _close_DSet(DSet *dset)
{
	if (dset->h5nchunk != NULL)
		free(dset->h5nchunk);
	if (dset->h5chunk_spacings != NULL)
		free(dset->h5chunk_spacings);
	if (dset->h5dim != NULL)
		free(dset->h5dim);
	if (dset->plist_id != -1)
		H5Pclose(dset->plist_id);
	if (dset->space_id != -1)
		H5Sclose(dset->space_id);
	if (dset->dtype_id != -1)
		H5Tclose(dset->dtype_id);
	if (dset->storage_mode_attr != NULL)
		free(dset->storage_mode_attr);
	H5Dclose(dset->dset_id);
}

int _get_DSet(hid_t dset_id, int as_int, int get_Rtype_only, int ndim,
	      DSet *dset)
{
	char *storage_mode_attr;
	hid_t dtype_id, space_id, plist_id, mem_type_id;
	H5T_class_t H5class;
	size_t size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE Rtype;
	int dset_ndim, *h5nchunk, h5along;
	hsize_t *h5dim, *h5chunk_spacings, d, spacing, nchunk;
	htri_t ret;

	dset->dset_id = dset_id;

	/* Initialize the fields that _close_DSet() will look at. */
	dset->storage_mode_attr = NULL;
	dset->dtype_id = -1;
	dset->space_id = -1;
	dset->plist_id = -1;
	dset->h5dim = NULL;
	dset->h5chunk_spacings = NULL;
	dset->h5nchunk = NULL;

	/* Set 'storage_mode_attr'. */
	ret = H5Aexists(dset_id, "storage.mode");
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aexists() returned an error");
		goto on_error;
	}
	if (ret > 0) {
		storage_mode_attr = get_storage_mode_attr(dset_id);
		if (storage_mode_attr == NULL)
			goto on_error;
		dset->storage_mode_attr = storage_mode_attr;
	}

	/* Set 'dtype_id'. */
	dtype_id = H5Dget_type(dset_id);
	if (dtype_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		goto on_error;
	}
	dset->dtype_id = dtype_id;

	/* Set 'H5class'. */
	H5class = H5Tget_class(dtype_id);
	if (H5class == H5T_NO_CLASS) {
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		goto on_error;
	}
	dset->H5class = H5class;

	/* Set 'size'. */
	size = H5Tget_size(dtype_id);
	if (size == 0) {
		PRINT_TO_ERRMSG_BUF("H5Tget_size() returned 0");
		goto on_error;
	}
	dset->size = size;

	/* Set 'Rtype'. */
	if (dset->storage_mode_attr != NULL) {
		if (map_storage_mode_to_Rtype(storage_mode_attr, as_int,
					      &Rtype) < 0)
			goto on_error;
	} else {
		if (map_H5class_to_Rtype(H5class, as_int, size, &Rtype) < 0)
			goto on_error;
	}
	if (Rtype == STRSXP && H5Tis_variable_str(dtype_id) != 0) {
		PRINT_TO_ERRMSG_BUF("reading variable-length string data "
				    "is not supported at the moment");
		goto on_error;
	}
	dset->Rtype = Rtype;

	if (get_Rtype_only) {
		_close_DSet(dset);
		return 0;
	}

	/* Set 'space_id'. */
	space_id = H5Dget_space(dset_id);
	if (space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		goto on_error;
	}
	dset->space_id = space_id;

	/* Set 'ndim'. */
	dset_ndim = H5Sget_simple_extent_ndims(space_id);
	if (dset_ndim < 0) {
		PRINT_TO_ERRMSG_BUF(
			"H5Sget_simple_extent_ndims() returned an error");
		goto on_error;
	}
	if (ndim != dset_ndim) {
		PRINT_TO_ERRMSG_BUF(
			"Dataset has %d dimensions but 'starts' has %d list "
			"element%s.\n  'starts' must have one list element "
			"per dimension in the dataset.",
			dset_ndim, ndim, ndim > 1 ? "s" : "");
		goto on_error;
	}
	dset->ndim = ndim;

	/* Set 'plist_id'. */
	plist_id = H5Dget_create_plist(dset_id);
	if (plist_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		goto on_error;
	}
	dset->plist_id = plist_id;

	/* Set 'h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		goto on_error;
	if (H5Sget_simple_extent_dims(space_id, h5dim, NULL) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		goto on_error;
	}
	dset->h5dim = h5dim;

	/* Set 'h5chunk_spacings'. */
	if (H5Pget_layout(plist_id) == H5D_CHUNKED) {
		h5chunk_spacings = _alloc_hsize_t_buf(ndim, 0,
						      "'h5chunk_spacings'");
		if (h5chunk_spacings == NULL)
			goto on_error;
		if (H5Pget_chunk(plist_id, ndim, h5chunk_spacings) != ndim) {
			PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
					    "an unexpected value");
			goto on_error;
		}
		dset->h5chunk_spacings = h5chunk_spacings;
	}

	/* Set 'h5nchunk'. */
	if (dset->h5chunk_spacings != NULL) {
		h5nchunk = (int *) malloc(ndim * sizeof(int));
		if (h5nchunk == NULL) {
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for 'h5nchunk'");
			goto on_error;
		}
		for (h5along = 0; h5along < ndim; h5along++) {
			d = h5dim[h5along];
			if (d == 0) {
				h5nchunk[h5along] = 0;
				continue;
			}
			spacing = dset->h5chunk_spacings[h5along];
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
	}

	/* Set 'ans_elt_size'. */
	ans_elt_size = get_ans_elt_size_from_Rtype(Rtype, size);
	if (ans_elt_size == 0)
		goto on_error;
	dset->ans_elt_size = ans_elt_size;

	/* Set 'chunk_data_buf_size'. */
	if (dset->h5chunk_spacings != NULL) {
		chunk_data_buf_size = ans_elt_size;
		for (h5along = 0; h5along < ndim; h5along++)
			chunk_data_buf_size *= dset->h5chunk_spacings[h5along];
		dset->chunk_data_buf_size = chunk_data_buf_size;
	}

	/* Set 'mem_type_id'. */
	mem_type_id = get_mem_type_id_from_Rtype(Rtype, dtype_id);
	if (mem_type_id < 0)
		goto on_error;
	dset->mem_type_id = mem_type_id;
	return 0;

    on_error:
	_close_DSet(dset);
	return -1;
}


/****************************************************************************
 * Convenience wrappers to H5Fopen() and H5Dopen(), with argument checking
 *
 * These are called at the very beginning of the various .Call entry points
 * where they are used (and before any resource is allocated) so it's ok to
 * error() immediately in case of error.
 */

hid_t _get_file_id(SEXP filepath)
{
	SEXP filepath0;
	hid_t file_id;

	if (!(IS_CHARACTER(filepath) && LENGTH(filepath) == 1))
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' cannot be NA");
	file_id = H5Fopen(CHAR(filepath0), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		error("failed to open file '%s'", CHAR(filepath0));
	return file_id;
}

hid_t _get_dset_id(hid_t file_id, SEXP name, SEXP filepath)
{
	SEXP name0;
	hid_t dset_id;

	if (!(IS_CHARACTER(name) && LENGTH(name) == 1))
		error("'name' must be a single string");
	name0 = STRING_ELT(name, 0);
	if (name0 == NA_STRING)
		error("'name' cannot be NA");
	dset_id = H5Dopen(file_id, CHAR(name0), H5P_DEFAULT);
	if (dset_id < 0) {
		H5Fclose(file_id);
		error("failed to open dataset '%s' from file '%s'",
		      CHAR(name0), CHAR(STRING_ELT(filepath, 0)));
	}
	return dset_id;
}


/****************************************************************************
 * C_get_h5mread_returned_type()
 *
 * The R type returned by h5mread() is determined by arguments 'filepath',
 * 'name', and 'as_integer'.
 */

/* --- .Call ENTRY POINT --- */
SEXP C_get_h5mread_returned_type(SEXP filepath, SEXP name, SEXP as_integer)
{
	int as_int;
	hid_t file_id, dset_id;
	DSet dset;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	file_id = _get_file_id(filepath);
	dset_id = _get_dset_id(file_id, name, filepath);
	/* 4th arg ('ndim') is ignored when 3rd arg ('get_Rtype_only') is
	   set to 1. _get_DSet() will do H5Dclose(dset_id). */
	if (_get_DSet(dset_id, as_int, 1, 0, &dset) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf);
	}
	H5Fclose(file_id);
	return ScalarString(type2str(dset.Rtype));
}

