/****************************************************************************
 *              Basic manipulation of a H5DSetDescriptor struct             *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for strcmp */
#include <limits.h>  /* for INT_MAX */


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

static const char *H5layout2str(H5D_layout_t H5layout)
{
	static char s[32];

	switch (H5layout) {
	    case H5D_COMPACT:    return "H5D_COMPACT";
	    case H5D_CONTIGUOUS: return "H5D_CONTIGUOUS";
	    case H5D_CHUNKED:    return "H5D_CHUNKED";
	    case H5D_VIRTUAL:    return "H5D_VIRTUAL";
	    default: break;
	}
	sprintf(s, "unknown (%d)", H5layout);
	return s;
}

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

/* Does not distinguish between no name and an empty name (i.e. a name
   that is the empty string). */
static char *get_h5name(hid_t obj_id)
{
	ssize_t name_size;
	char *h5name;

	name_size = H5Iget_name(obj_id, NULL, 0);
	if (name_size < 0) {
		PRINT_TO_ERRMSG_BUF("H5Iget_name() returned an error");
		return NULL;
	}
	name_size++;
	h5name = (char *) malloc((size_t) name_size);
	if (h5name == NULL) {
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'h5name'");
		return NULL;
	}
	name_size = H5Iget_name(obj_id, h5name, (size_t) name_size);
	if (name_size < 0) {
		PRINT_TO_ERRMSG_BUF("H5Iget_name() returned an error");
		return NULL;
	}
	return h5name;
}

/* Get the value (expected to be a string) of a given attribute of class
   H5T_STRING. Return:
    -1: if an error occurs;
     0: if dataset 'dset_id' has no attribute with the name specified
        in 'attr_name';
     1: if dataset 'dset_id' has an attribute with the name specified
        in 'attr_name' but the attribute is not of class H5T_STRING;
     2: if dataset 'dset_id' has an attribute with the name specified
        in 'attr_name' and the attribute is of class H5T_STRING (in which
        case, and only in this case, the value of the attribute is copied
        to 'val').
 */
int _get_h5attrib_strval(hid_t dset_id, const char *attr_name, CharAE *val)
{
	int ret;
	hid_t attr_id, attr_type_id;
	H5T_class_t attr_H5class;
	hsize_t attr_size;

	ret = H5Aexists(dset_id, attr_name);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aexists() returned an error");
		return -1;
	}
	if (ret == 0)
		return 0;
	attr_id = H5Aopen(dset_id, attr_name, H5P_DEFAULT);
	if (attr_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aopen() returned an error");
		return -1;
	}
	attr_type_id = H5Aget_type(attr_id);
	if (attr_type_id < 0) {
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Aget_type() returned an error");
		return -1;
	}
	attr_H5class = H5Tget_class(attr_type_id);
	if (attr_H5class == H5T_NO_CLASS) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		return -1;
	}
	if (attr_H5class != H5T_STRING) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		return 1;
	}
	attr_size = H5Aget_storage_size(attr_id);
	if (attr_size == 0) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Aget_storage_size() returned 0");
		return -1;
	}
	if ((size_t) attr_size > val->_buflength)
		CharAE_extend(val, (size_t) attr_size);
	CharAE_set_nelt(val, (size_t) attr_size);
	ret = H5Aread(attr_id, attr_type_id, val->elts);
	H5Tclose(attr_type_id);
	H5Aclose(attr_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aread() returned an error");
		return -1;
	}
	return 2;
}

/* Get the value (expected to be a single int) of a given attribute of
   class H5T_INTEGER. Return:
    -1: if an error occurs;
     0: if dataset 'dset_id' has no attribute with the name specified
        in 'attr_name';
     1: if dataset 'dset_id' has an attribute with the name specified
        in 'attr_name' but the attribute is not of class H5T_INTEGER
        or its value is not a single int;
     2: if dataset 'dset_id' has an attribute with the name specified
        in 'attr_name' and the attribute is of class H5T_INTEGER and
        its value is a single int (in which case, and only in this case,
        the value of the attribute is copied to 'val').
 */
static int get_h5attrib_intval(hid_t dset_id, const char *attr_name, int *val)
{
	int ret;
	hid_t attr_id, attr_type_id;
	H5T_class_t attr_H5class;
	hsize_t attr_size;

	ret = H5Aexists(dset_id, attr_name);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aexists() returned an error");
		return -1;
	}
	if (ret == 0)
		return 0;
	attr_id = H5Aopen(dset_id, attr_name, H5P_DEFAULT);
	if (attr_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aopen() returned an error");
		return -1;
	}
	attr_type_id = H5Aget_type(attr_id);
	if (attr_type_id < 0) {
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Aget_type() returned an error");
		return -1;
	}
	attr_H5class = H5Tget_class(attr_type_id);
	if (attr_H5class == H5T_NO_CLASS) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		return -1;
	}
	if (attr_H5class != H5T_INTEGER) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		return 1;
	}
	attr_size = H5Aget_storage_size(attr_id);
	if (attr_size != sizeof(int)) {
		H5Tclose(attr_type_id);
		H5Aclose(attr_id);
		return 1;
	}
	ret = H5Aread(attr_id, attr_type_id, val);
	H5Tclose(attr_type_id);
	H5Aclose(attr_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Aread() returned an error");
		return -1;
	}
	return 2;
}


/****************************************************************************
 * _init_H5DSetDescriptor() / _destroy_H5DSetDescriptor()
 */

static int map_storage_mode_to_Rtype(const char *storage_mode, int as_int,
				     SEXPTYPE *Rtype)
{
	if (strcmp(storage_mode, "logical") == 0) {
		*Rtype = as_int ? INTSXP : LGLSXP;
		return 0;
	}
	if (strcmp(storage_mode, "integer") == 0) {
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
	if (strcmp(storage_mode, "raw") == 0) {
		*Rtype = as_int ? INTSXP : RAWSXP;
		return 0;
	}
	PRINT_TO_ERRMSG_BUF("the dataset to read has a \"storage.mode\" "
			    "attribute set to unsupported value: \"%s\"",
			    storage_mode);
	return -1;
}

/* See hdf5-1.10.3/src/H5Tpublic.h for the list of datatype classes. We only
   support H5T_INTEGER, H5T_FLOAT, and H5T_STRING for now. */
static int map_H5class_to_Rtype(H5T_class_t H5class, int as_int, size_t H5size,
				SEXPTYPE *Rtype)
{
	switch (H5class) {
	    case H5T_INTEGER:
		if (as_int) {
			*Rtype = INTSXP;
		} else if (H5size <= sizeof(char)) {
			*Rtype = RAWSXP;
		} else if (H5size <= sizeof(int)) {
			*Rtype = INTSXP;
		} else {
			*Rtype = REALSXP;
		}
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

static size_t get_ans_elt_size_from_Rtype(SEXPTYPE Rtype, size_t H5size)
{
	switch (Rtype) {
	    case LGLSXP:
	    case INTSXP:  return sizeof(int);
	    case REALSXP: return sizeof(double);
	    case STRSXP:  return H5size;
	    case RAWSXP:  return sizeof(char);
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
	    case RAWSXP:  return H5T_NATIVE_UCHAR;
	}
	PRINT_TO_ERRMSG_BUF("unsupported type: %s", CHAR(type2str(Rtype)));
	return -1;
}

void _destroy_H5DSetDescriptor(H5DSetDescriptor *h5dset)
{
	if (h5dset->h5nchunk != NULL)
		free(h5dset->h5nchunk);
	if (h5dset->h5chunkdim != NULL &&
	    h5dset->h5chunkdim != h5dset->h5dim)
		free(h5dset->h5chunkdim);
	if (h5dset->h5dim != NULL)
		free(h5dset->h5dim);
	if (h5dset->plist_id != -1)
		H5Pclose(h5dset->plist_id);
	if (h5dset->space_id != -1)
		H5Sclose(h5dset->space_id);
	if (h5dset->dtype_id != -1)
		H5Tclose(h5dset->dtype_id);
	if (h5dset->storage_mode_attr != NULL)
		free(h5dset->storage_mode_attr);
	if (h5dset->h5name != NULL)
		free(h5dset->h5name);
	return;
}

int _init_H5DSetDescriptor(H5DSetDescriptor *h5dset, hid_t dset_id,
		     int as_int, int get_Rtype_only)
{
	char *h5name, *storage_mode_attr;
	hid_t dtype_id, space_id, plist_id, mem_type_id;
	H5T_class_t H5class;
	size_t H5size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE Rtype;
	int as_na_attr, ndim, *h5nchunk, h5along;
	hsize_t *h5dim, *h5chunkdim, d, chunkd, nchunk;
	htri_t ret;
	CharAE *buf;

	h5dset->dset_id = dset_id;

	/* Initialize the fields that _destroy_H5DSetDescriptor() will free
	   or close. */
	h5dset->h5name = NULL;
	h5dset->storage_mode_attr = NULL;
	h5dset->dtype_id = -1;
	h5dset->space_id = -1;
	h5dset->plist_id = -1;
	h5dset->h5dim = NULL;
	h5dset->h5chunkdim = NULL;
	h5dset->h5nchunk = NULL;

	/* Set 'h5dset->h5name'. */
	h5name = get_h5name(dset_id);
	if (h5name == NULL)
		goto on_error;
	h5dset->h5name = h5name;

	/* Set 'h5dset->storage_mode_attr'. */
	buf = new_CharAE(0);
	ret = _get_h5attrib_strval(dset_id, "storage.mode", buf);
	if (ret < 0)
		goto on_error;
	if (ret == 1) {
		PRINT_TO_ERRMSG_BUF("attribute \"storage.mode\" is "
				    "not of expected class H5T_STRING");
		goto on_error;
	}
	if (ret == 2) {
		storage_mode_attr = (char *) malloc(CharAE_get_nelt(buf));
		if (storage_mode_attr == NULL) {
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for 'storage_mode_attr'");
			goto on_error;
		}
		strcpy(storage_mode_attr, buf->elts);
		h5dset->storage_mode_attr = storage_mode_attr;
	}

	/* Set 'h5dset->dtype_id'. */
	dtype_id = H5Dget_type(dset_id);
	if (dtype_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		goto on_error;
	}
	h5dset->dtype_id = dtype_id;

	/* Set 'h5dset->H5class'. */
	H5class = H5Tget_class(dtype_id);
	if (H5class == H5T_NO_CLASS) {
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		goto on_error;
	}
	h5dset->H5class = H5class;

	/* Set 'h5dset->H5size'. */
	H5size = H5Tget_size(dtype_id);
	if (H5size == 0) {
		PRINT_TO_ERRMSG_BUF("H5Tget_size() returned 0");
		goto on_error;
	}
	h5dset->H5size = H5size;

	/* Set 'h5dset->Rtype'. */
	if (h5dset->storage_mode_attr != NULL) {
		if (map_storage_mode_to_Rtype(storage_mode_attr, as_int,
					      &Rtype) < 0)
			goto on_error;
	} else {
		if (map_H5class_to_Rtype(H5class, as_int, H5size, &Rtype) < 0)
			goto on_error;
	}
	if (Rtype == STRSXP && H5Tis_variable_str(dtype_id) != 0) {
		PRINT_TO_ERRMSG_BUF("reading variable-length string data "
				    "is not supported at the moment");
		goto on_error;
	}
	h5dset->Rtype = Rtype;

	if (get_Rtype_only)
		return 0;

	/* Set 'h5dset->as_na_attr'. */
	ret = get_h5attrib_intval(dset_id, "as.na", &as_na_attr);
	if (ret < 0)
		goto on_error;
	if (ret == 1) {
		PRINT_TO_ERRMSG_BUF("attribute \"as.na\" is "
				    "not of expected class H5T_INTEGER"
				    "or its value is not a single int");
		goto on_error;
	}
	h5dset->as_na_attr = ret == 2 ? as_na_attr : 0;

	/* Set 'h5dset->space_id'. */
	space_id = H5Dget_space(dset_id);
	if (space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		goto on_error;
	}
	h5dset->space_id = space_id;

	/* Set 'h5dset->ndim'. */
	ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim < 0) {
		PRINT_TO_ERRMSG_BUF(
			"H5Sget_simple_extent_ndims() returned an error");
		goto on_error;
	}
	h5dset->ndim = ndim;

	/* Set 'h5dset->plist_id'. */
	plist_id = H5Dget_create_plist(dset_id);
	if (plist_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		goto on_error;
	}
	h5dset->plist_id = plist_id;

	/* Set 'h5dset->h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		goto on_error;
	if (H5Sget_simple_extent_dims(space_id, h5dim, NULL) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		goto on_error;
	}
	h5dset->h5dim = h5dim;

	/* Set 'h5dset->H5layout'. */
	h5dset->H5layout = H5Pget_layout(plist_id);

	/* Set 'h5dset->h5chunkdim'. */
	if (h5dset->H5layout == H5D_CHUNKED) {
		h5chunkdim = _alloc_hsize_t_buf(ndim, 0, "'h5chunkdim'");
		if (h5chunkdim == NULL)
			goto on_error;
		if (H5Pget_chunk(plist_id, ndim, h5chunkdim) != ndim) {
			PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
					    "an unexpected value");
			goto on_error;
		}
		h5dset->h5chunkdim = h5chunkdim;
	} else if (h5dset->Rtype == STRSXP) {
		/* Even though the dataset is contiguous, we treat it as
		   if it was made of a single chunk. This is so we can
		   use h5mread() methods 4 or 5 on it, which work only
		   on chunked data and are the only methods that know
		   how to handle string data. */
		h5dset->h5chunkdim = h5dset->h5dim;
	}

	/* Set 'h5dset->h5nchunk'. */
	if (h5dset->h5chunkdim != NULL) {
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
			chunkd = h5dset->h5chunkdim[h5along];
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
		h5dset->h5nchunk = h5nchunk;
	}

	/* Set 'h5dset->ans_elt_size'. */
	ans_elt_size = get_ans_elt_size_from_Rtype(Rtype, H5size);
	if (ans_elt_size == 0)
		goto on_error;
	h5dset->ans_elt_size = ans_elt_size;

	/* Set 'h5dset->chunk_data_buf_size'. */
	if (h5dset->h5chunkdim != NULL) {
		chunk_data_buf_size = ans_elt_size;
		for (h5along = 0; h5along < ndim; h5along++)
			chunk_data_buf_size *=
				h5dset->h5chunkdim[h5along];
		h5dset->chunk_data_buf_size = chunk_data_buf_size;
	}

	/* Set 'h5dset->mem_type_id'. */
	mem_type_id = get_mem_type_id_from_Rtype(Rtype, dtype_id);
	if (mem_type_id < 0)
		goto on_error;
	h5dset->mem_type_id = mem_type_id;
	return 0;

    on_error:
	_destroy_H5DSetDescriptor(h5dset);
	return -1;
}


/****************************************************************************
 * Convenience wrappers to H5Fopen() and H5Dopen(), with argument checking
 *
 * These are called at the very beginning of the various .Call entry points
 * where they are used (and before any resource is allocated) so it's ok to
 * error() immediately in case of error.
 */

hid_t _get_file_id(SEXP filepath, int readonly)
{
	SEXP filepath0;
	herr_t ret;
	unsigned int flags;
	hid_t file_id;

	if (!(IS_CHARACTER(filepath) && LENGTH(filepath) == 1))
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' cannot be NA");
	ret = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	if (ret < 0)
		error("H5Eset_auto() returned an error");
	flags = readonly ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
	file_id = H5Fopen(CHAR(filepath0), flags, H5P_DEFAULT);
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
 * Used in R/H5DSetDescriptor-class.R
 */

/* --- .Call ENTRY POINT --- */
SEXP C_destroy_H5DSetDescriptor_xp(SEXP xp)
{
	H5DSetDescriptor *h5dset;

	h5dset = R_ExternalPtrAddr(xp);
	if (h5dset != NULL) {
		//printf("Destroying H5DSetDescriptor struct at address %p ... ",
		//       h5dset);
		_destroy_H5DSetDescriptor(h5dset);
		H5Dclose(h5dset->dset_id);
		free(h5dset);
		R_SetExternalPtrAddr(xp, NULL);
		//printf("OK\n");
	}

	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_new_H5DSetDescriptor_xp(SEXP filepath, SEXP name, SEXP as_integer)
{
	int as_int;
	hid_t file_id, dset_id;
	H5DSetDescriptor *h5dset;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);

	h5dset = (H5DSetDescriptor *) malloc(sizeof(H5DSetDescriptor));
	if (h5dset == NULL) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error("C_new_H5DSetDescriptor_xp(): malloc() failed");
	}

	if (_init_H5DSetDescriptor(h5dset, dset_id, as_int, 0) < 0) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}
	H5Fclose(file_id);
	//printf("H5DSetDescriptor struct created at address %p\n", h5dset);

	return R_MakeExternalPtr(h5dset, R_NilValue, R_NilValue);
}

/* --- .Call ENTRY POINT --- */
SEXP C_show_H5DSetDescriptor_xp(SEXP xp)
{
	const H5DSetDescriptor *h5dset;
	int h5along;

	h5dset = R_ExternalPtrAddr(xp);
	if (h5dset == NULL) {
		Rprintf("Expired H5DSetDescriptor\n");
		return R_NilValue;
	}

	Rprintf("H5DSetDescriptor:\n");
	Rprintf("- dset_id = %lu\n", h5dset->dset_id);

	Rprintf("- h5name = \"%s\"\n", h5dset->h5name);

	Rprintf("- storage_mode_attr = ");
	if (h5dset->storage_mode_attr == NULL) {
		Rprintf("NULL");
	} else {
		Rprintf("\"%s\"", h5dset->storage_mode_attr);
	}
	Rprintf("\n");

	Rprintf("- dtype_id = %lu\n", h5dset->dtype_id);

	Rprintf("- H5class = %s\n", H5class2str(h5dset->H5class));

	Rprintf("- H5size = %lu\n", h5dset->H5size);

	Rprintf("- Rtype = \"%s\"\n", CHAR(type2str(h5dset->Rtype)));

	Rprintf("- as_na_attr = %d\n", h5dset->as_na_attr);

	Rprintf("- space_id = %lu\n", h5dset->space_id);

	Rprintf("- ndim = %d\n", h5dset->ndim);

	Rprintf("- plist_id = %lu\n", h5dset->plist_id);

	Rprintf("- h5dim =");
	for (h5along = 0; h5along < h5dset->ndim; h5along++)
		Rprintf(" %llu", h5dset->h5dim[h5along]);
	Rprintf("\n");

	Rprintf("- H5layout = %s\n", H5layout2str(h5dset->H5layout));

	Rprintf("- h5chunkdim =");
	if (h5dset->h5chunkdim == NULL) {
		Rprintf(" NULL\n");
	} else {
		for (h5along = 0; h5along < h5dset->ndim; h5along++)
			Rprintf(" %llu",
				h5dset->h5chunkdim[h5along]);
		if (h5dset->H5layout != H5D_CHUNKED &&
		    h5dset->h5chunkdim == h5dset->h5dim)
			Rprintf(" (artificially set to h5dim)");
		Rprintf("\n");
		Rprintf("    h5nchunk =");
		for (h5along = 0; h5along < h5dset->ndim; h5along++)
			Rprintf(" %d", h5dset->h5nchunk[h5along]);
		Rprintf("\n");
		Rprintf("    chunk_data_buf_size = %lu\n",
			h5dset->chunk_data_buf_size);
	}

	Rprintf("- ans_elt_size = %lu\n", h5dset->ans_elt_size);

	Rprintf("- mem_type_id = %lu\n", h5dset->mem_type_id);

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
	int as_int, ret;
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	ret = _init_H5DSetDescriptor(&h5dset, dset_id, as_int, 1);
	/* It's ok to close 'dset_id' **before** destroying its descriptor. */
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf());
	_destroy_H5DSetDescriptor(&h5dset);
	return ScalarString(type2str(h5dset.Rtype));
}

