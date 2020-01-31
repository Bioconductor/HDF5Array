/****************************************************************************
 *                 Basic manipulation of a DSetHandle struct                *
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
 * _get_DSetHandle() / _close_DSetHandle()
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

static const char *H5class2str(H5T_class_t H5class)
{
	static char s[32];

	switch (H5class) {
	    case H5T_INTEGER:   return "H5T_INTEGER";
	    case H5T_FLOAT:     return "H5T_FLOAT";
	    case H5T_STRING:    return "H5T_STRING";
	    case H5T_TIME:      return "H5T_TIME";
	    case H5T_BITFIELD:  return "H5T_BITFIELD";
	    case H5T_OPAQUE:    return "H5T_OPAQUE";
	    case H5T_COMPOUND:  return "H5T_COMPOUND";
	    case H5T_REFERENCE: return "H5T_REFERENCE";
	    case H5T_ENUM:      return "H5T_ENUM";
	    case H5T_VLEN:      return "H5T_VLEN";
	    case H5T_ARRAY:     return "H5T_ARRAY";
	    default: break;
	}
	sprintf(s, "unknown (%d)", H5class);
	return s;
}

/* See hdf5-1.10.3/src/H5Tpublic.h for the list of datatype classes. We only
   support H5T_INTEGER, H5T_FLOAT, and H5T_STRING for now. */
static int map_H5class_to_Rtype(H5T_class_t H5class, int as_int, size_t size,
				SEXPTYPE *Rtype)
{
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
	    default: break;
	}
	PRINT_TO_ERRMSG_BUF("unsupported dataset class: %s",
			    H5class2str(H5class));
	return -1;
}

static const char *layout2str(H5D_layout_t layout)
{
	static char s[32];

	switch (layout) {
	    case H5D_COMPACT:    return "H5D_COMPACT";
	    case H5D_CONTIGUOUS: return "H5D_CONTIGUOUS";
	    case H5D_CHUNKED:    return "H5D_CHUNKED";
	    case H5D_VIRTUAL:    return "H5D_VIRTUAL";
	    default: break;
	}
	sprintf(s, "unknown (%d)", layout);
	return s;
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

void _close_DSetHandle(DSetHandle *dset_handle)
{
	if (dset_handle->h5nchunk != NULL)
		free(dset_handle->h5nchunk);
	if (dset_handle->h5chunkdim != NULL &&
	    dset_handle->h5chunkdim != dset_handle->h5dim)
		free(dset_handle->h5chunkdim);
	if (dset_handle->h5dim != NULL)
		free(dset_handle->h5dim);
	if (dset_handle->plist_id != -1)
		H5Pclose(dset_handle->plist_id);
	if (dset_handle->space_id != -1)
		H5Sclose(dset_handle->space_id);
	if (dset_handle->dtype_id != -1)
		H5Tclose(dset_handle->dtype_id);
	if (dset_handle->storage_mode_attr != NULL)
		free(dset_handle->storage_mode_attr);
	H5Dclose(dset_handle->dset_id);
}

int _get_DSetHandle(hid_t dset_id, int as_int, int get_Rtype_only,
		    DSetHandle *dset_handle)
{
	char *storage_mode_attr;
	hid_t dtype_id, space_id, plist_id, mem_type_id;
	H5T_class_t H5class;
	size_t size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE Rtype;
	int ndim, *h5nchunk, h5along;
	hsize_t *h5dim, *h5chunkdim, d, chunkd, nchunk;
	htri_t ret;

	dset_handle->dset_id = dset_id;

	/* Initialize the fields that _close_DSetHandle() will look at. */
	dset_handle->storage_mode_attr = NULL;
	dset_handle->dtype_id = -1;
	dset_handle->space_id = -1;
	dset_handle->plist_id = -1;
	dset_handle->h5dim = NULL;
	dset_handle->h5chunkdim = NULL;
	dset_handle->h5nchunk = NULL;

	/* Set 'dset_handle->storage_mode_attr'. */
	ret = H5Aexists(dset_id, "storage.mode");
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aexists() returned an error");
		goto on_error;
	}
	if (ret > 0) {
		storage_mode_attr = get_storage_mode_attr(dset_id);
		if (storage_mode_attr == NULL)
			goto on_error;
		dset_handle->storage_mode_attr = storage_mode_attr;
	}

	/* Set 'dset_handle->dtype_id'. */
	dtype_id = H5Dget_type(dset_id);
	if (dtype_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		goto on_error;
	}
	dset_handle->dtype_id = dtype_id;

	/* Set 'dset_handle->H5class'. */
	H5class = H5Tget_class(dtype_id);
	if (H5class == H5T_NO_CLASS) {
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		goto on_error;
	}
	dset_handle->H5class = H5class;

	/* Set 'dset_handle->size'. */
	size = H5Tget_size(dtype_id);
	if (size == 0) {
		PRINT_TO_ERRMSG_BUF("H5Tget_size() returned 0");
		goto on_error;
	}
	dset_handle->size = size;

	/* Set 'dset_handle->Rtype'. */
	if (dset_handle->storage_mode_attr != NULL) {
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
	dset_handle->Rtype = Rtype;

	if (get_Rtype_only) {
		_close_DSetHandle(dset_handle);
		return 0;
	}

	/* Set 'dset_handle->space_id'. */
	space_id = H5Dget_space(dset_id);
	if (space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		goto on_error;
	}
	dset_handle->space_id = space_id;

	/* Set 'dset_handle->ndim'. */
	ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim < 0) {
		PRINT_TO_ERRMSG_BUF(
			"H5Sget_simple_extent_ndims() returned an error");
		goto on_error;
	}
	dset_handle->ndim = ndim;

	/* Set 'dset_handle->plist_id'. */
	plist_id = H5Dget_create_plist(dset_id);
	if (plist_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		goto on_error;
	}
	dset_handle->plist_id = plist_id;

	/* Set 'dset_handle->h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		goto on_error;
	if (H5Sget_simple_extent_dims(space_id, h5dim, NULL) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		goto on_error;
	}
	dset_handle->h5dim = h5dim;

	/* Set 'dset_handle->layout'. */
	dset_handle->layout = H5Pget_layout(plist_id);

	/* Set 'dset_handle->h5chunkdim'. */
	if (dset_handle->layout == H5D_CHUNKED) {
		h5chunkdim = _alloc_hsize_t_buf(ndim, 0, "'h5chunkdim'");
		if (h5chunkdim == NULL)
			goto on_error;
		if (H5Pget_chunk(plist_id, ndim, h5chunkdim) != ndim) {
			PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
					    "an unexpected value");
			goto on_error;
		}
		dset_handle->h5chunkdim = h5chunkdim;
	} else if (dset_handle->Rtype == STRSXP) {
		/* Even though the dataset is contiguous, we treat it as
		   if it was made of a single chunk. This is so we can
		   use h5mread() methods 4 or 5 on it, which work only
		   on chunked data and are the only methods that know
		   how to handle string data. */
		dset_handle->h5chunkdim = dset_handle->h5dim;
	}

	/* Set 'dset_handle->h5nchunk'. */
	if (dset_handle->h5chunkdim != NULL) {
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
			chunkd = dset_handle->h5chunkdim[h5along];
			nchunk = d / chunkd;
			if (d % chunkd != 0)
				nchunk++;
			if (nchunk > INT_MAX) {
				PRINT_TO_ERRMSG_BUF("datasets with more than "
					"INT_MAX chunks along any dimension\n  "
					"are not supported at the moment");
				goto on_error;
			}
			h5nchunk[h5along] = nchunk;
		}
		dset_handle->h5nchunk = h5nchunk;
	}

	/* Set 'dset_handle->ans_elt_size'. */
	ans_elt_size = get_ans_elt_size_from_Rtype(Rtype, size);
	if (ans_elt_size == 0)
		goto on_error;
	dset_handle->ans_elt_size = ans_elt_size;

	/* Set 'dset_handle->chunk_data_buf_size'. */
	if (dset_handle->h5chunkdim != NULL) {
		chunk_data_buf_size = ans_elt_size;
		for (h5along = 0; h5along < ndim; h5along++)
			chunk_data_buf_size *=
				dset_handle->h5chunkdim[h5along];
		dset_handle->chunk_data_buf_size = chunk_data_buf_size;
	}

	/* Set 'dset_handle->mem_type_id'. */
	mem_type_id = get_mem_type_id_from_Rtype(Rtype, dtype_id);
	if (mem_type_id < 0)
		goto on_error;
	dset_handle->mem_type_id = mem_type_id;
	return 0;

    on_error:
	_close_DSetHandle(dset_handle);
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
 * Used in R/DSetHandle-class.R
 */

/* --- .Call ENTRY POINT --- */
SEXP C_destroy_DSetHandle_xp(SEXP xp)
{
	DSetHandle *dset_handle;

	dset_handle = R_ExternalPtrAddr(xp);
	if (dset_handle != NULL) {
		//printf("Destroying DSetHandle struct at address %p ... ",
		//       dset_handle);
		_close_DSetHandle(dset_handle);
		free(dset_handle);
		R_SetExternalPtrAddr(xp, NULL);
		//printf("OK\n");
	}

	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_create_DSetHandle_xp(SEXP filepath, SEXP name, SEXP as_integer)
{
	int as_int;
	hid_t file_id, dset_id;
	DSetHandle *dset_handle;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	file_id = _get_file_id(filepath);
	dset_id = _get_dset_id(file_id, name, filepath);

	dset_handle = (DSetHandle *) malloc(sizeof(DSetHandle));
	if (dset_handle == NULL) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error("C_create_DSetHandle_xp(): malloc() failed");
	}

	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, as_int, 0, dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}
	H5Fclose(file_id);
	//printf("DSetHandle struct created at address %p\n", dset_handle);

	return R_MakeExternalPtr(dset_handle, R_NilValue, R_NilValue);
}

/* --- .Call ENTRY POINT --- */
SEXP C_show_DSetHandle_xp(SEXP xp)
{
	const DSetHandle *dset_handle;
	int h5along;

	dset_handle = R_ExternalPtrAddr(xp);
	if (dset_handle == NULL) {
		Rprintf("Expired DSetHandle\n");
		return R_NilValue;
	}

	Rprintf("DSetHandle:\n");
	Rprintf("- dset_id = %lu\n", dset_handle->dset_id);

	Rprintf("- storage_mode_attr = ");
	if (dset_handle->storage_mode_attr == NULL) {
		Rprintf("NULL");
	} else {
		Rprintf("\"%s\"", dset_handle->storage_mode_attr);
	}
	Rprintf("\n");

	Rprintf("- dtype_id = %lu\n", dset_handle->dtype_id);

	Rprintf("- H5class = %s\n", H5class2str(dset_handle->H5class));

	Rprintf("- size = %lu\n", dset_handle->size);

	Rprintf("- Rtype = \"%s\"\n", CHAR(type2str(dset_handle->Rtype)));

	Rprintf("- space_id = %lu\n", dset_handle->space_id);

	Rprintf("- ndim = %d\n", dset_handle->ndim);

	Rprintf("- plist_id = %lu\n", dset_handle->plist_id);

	Rprintf("- h5dim =");
	for (h5along = 0; h5along < dset_handle->ndim; h5along++)
		Rprintf(" %llu", dset_handle->h5dim[h5along]);
	Rprintf("\n");

	Rprintf("- layout = %s\n", layout2str(dset_handle->layout));

	Rprintf("- h5chunkdim =");
	if (dset_handle->h5chunkdim == NULL) {
		Rprintf(" NULL\n");
	} else {
		for (h5along = 0; h5along < dset_handle->ndim; h5along++)
			Rprintf(" %llu",
				dset_handle->h5chunkdim[h5along]);
		if (dset_handle->layout != H5D_CHUNKED &&
		    dset_handle->h5chunkdim == dset_handle->h5dim)
			Rprintf(" (artificially set to h5dim)");
		Rprintf("\n");
		Rprintf("    h5nchunk =");
		for (h5along = 0; h5along < dset_handle->ndim; h5along++)
			Rprintf(" %d", dset_handle->h5nchunk[h5along]);
		Rprintf("\n");
		Rprintf("    chunk_data_buf_size = %lu\n",
			dset_handle->chunk_data_buf_size);
	}

	Rprintf("- ans_elt_size = %lu\n", dset_handle->ans_elt_size);

	Rprintf("- mem_type_id = %lu\n", dset_handle->mem_type_id);

	return R_NilValue;
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
	DSetHandle dset_handle;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	file_id = _get_file_id(filepath);
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, as_int, 1, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}
	H5Fclose(file_id);
	return ScalarString(type2str(dset_handle.Rtype));
}

