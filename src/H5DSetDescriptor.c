/****************************************************************************
 *              Basic manipulation of a H5DSetDescriptor struct             *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "H5DSetDescriptor.h"

#include "global_errmsg_buf.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for strcmp */
#include <limits.h>  /* for INT_MAX */


/* Make sure to keep map_native_type_to_predef_type() and
   predef_native_type_as_string() in sync! */
#define	RETURN_PREDEF_TYPE_IF(predef_type_id) \
	if (H5Tequal(native_type_id, predef_type_id) > 0) \
		return predef_type_id;
static hid_t map_native_type_to_predef_type(hid_t native_type_id)
{
	/* See H5Tpublic.h for authoritative list of native types.
	   Note that the list published at
	     https://portal.hdfgroup.org/display/HDF5/Predefined+Datatypes
	   might be incomplete or out-of-sync with H5Tpublic.h. */
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_CHAR)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_SCHAR)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UCHAR)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_SHORT)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_USHORT)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_LONG)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_ULONG)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_LLONG)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_ULLONG)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_FLOAT)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_DOUBLE)
#ifdef H5T_NATIVE_LDOUBLE
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_LDOUBLE)
#endif
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_B8)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_B16)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_B32)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_B64)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_OPAQUE)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_HADDR)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_HSIZE)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_HSSIZE)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_HERR)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_HBOOL)

	/* C9x integer types */

	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT8)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT8)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_LEAST8)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_LEAST8)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_FAST8)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_FAST8)

	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT16)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT16)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_LEAST16)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_LEAST16)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_FAST16)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_FAST16)

	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT32)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT32)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_LEAST32)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_LEAST32)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_FAST32)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_FAST32)

	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT64)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT64)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_LEAST64)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_LEAST64)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_INT_FAST64)
	RETURN_PREDEF_TYPE_IF(H5T_NATIVE_UINT_FAST64)

	PRINT_TO_ERRMSG_BUF("failed to map native type id %ld "
			    "to predef type id", native_type_id);
	return -1;
}

/* Make sure to keep map_native_type_to_predef_type() and
   predef_native_type_as_string() in sync! */
#define	RETURN_TYPE_STRING_IF(predef_type_id) \
	if (native_type_id == predef_type_id) \
		return #predef_type_id;  /* stringize 'predef_type_id' */
static const char *predef_native_type_as_string(hid_t native_type_id)
{
	static char s[50];

	RETURN_TYPE_STRING_IF(H5T_NATIVE_CHAR)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_SCHAR)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UCHAR)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_SHORT)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_USHORT)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_LONG)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_ULONG)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_LLONG)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_ULLONG)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_FLOAT)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_DOUBLE)
#ifdef H5T_NATIVE_LDOUBLE
	RETURN_TYPE_STRING_IF(H5T_NATIVE_LDOUBLE)
#endif
	RETURN_TYPE_STRING_IF(H5T_NATIVE_B8)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_B16)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_B32)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_B64)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_OPAQUE)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_HADDR)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_HSIZE)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_HSSIZE)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_HERR)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_HBOOL)

	/* C9x integer types */

	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT8)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT8)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_LEAST8)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_LEAST8)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_FAST8)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_FAST8)

	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT16)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT16)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_LEAST16)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_LEAST16)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_FAST16)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_FAST16)

	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT32)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT32)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_LEAST32)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_LEAST32)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_FAST32)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_FAST32)

	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT64)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT64)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_LEAST64)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_LEAST64)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_INT_FAST64)
	RETURN_TYPE_STRING_IF(H5T_NATIVE_UINT_FAST64)

	/* Should never happen if predef_native_type_as_string()
	   is kept in sync with map_native_type_to_predef_type(). */
	sprintf(s, "unknown native type (%ld)", native_type_id);
	return s;
}

static const char *H5class2str(H5T_class_t h5class)
{
	static char s[32];

	switch (h5class) {
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
	/* Only to avoid gcc warning about enumeration value not handled
	   in switch. */
	    default: break;
	}
	sprintf(s, "unknown (%d)", h5class);
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
	/* Only to avoid gcc warning about enumeration value not handled
	   in switch. */
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
 * Handle Rtype
 */

/* See hdf5-1.10.3/src/H5Tpublic.h for the list of datatype classes. We only
   support H5T_INTEGER, H5T_FLOAT, and H5T_STRING for now. */
static int map_numeric_h5class_to_Rtype(H5T_class_t h5class, size_t h5size,
		SEXPTYPE *Rtype)
{
	switch (h5class) {
	    case H5T_INTEGER:
		if (h5size <= sizeof(char)) {
			*Rtype = RAWSXP;
		} else if (h5size <= sizeof(int)) {
			*Rtype = INTSXP;
		} else {
			*Rtype = REALSXP;
		}
		return 0;
	    case H5T_FLOAT:
		*Rtype = REALSXP;
		return 0;
	/* Only to avoid gcc warning about enumeration value not handled
	   in switch. */
	    default: break;
	}
	PRINT_TO_ERRMSG_BUF("unsupported dataset class: %s",
			    H5class2str(h5class));
	return -1;
}

static int map_storage_mode_to_Rtype(const char *storage_mode, SEXPTYPE *Rtype)
{
	if (strcmp(storage_mode, "logical") == 0) {
		*Rtype = LGLSXP;
		return 0;
	}
	if (strcmp(storage_mode, "integer") == 0) {
		*Rtype = INTSXP;
		return 0;
	}
	if (strcmp(storage_mode, "double") == 0 ||
	    strcmp(storage_mode, "numeric") == 0) {
		*Rtype = REALSXP;
		return 0;
	}
	if (strcmp(storage_mode, "character") == 0) {
		warning("dataset has a \"storage.mode\" attribute "
			"set to \"character\" -- we ignore it");
		return 0;
	}
	if (strcmp(storage_mode, "raw") == 0) {
		*Rtype = RAWSXP;
		return 0;
	}
	PRINT_TO_ERRMSG_BUF("dataset has a \"storage.mode\" "
			    "attribute set to unsupported value: \"%s\"",
			    storage_mode);
	return -1;
}

static int set_Rtype(H5T_class_t h5class, size_t h5type_size,
		int as_int, const char *storage_mode, SEXPTYPE *Rtype)
{
	int ret;

	if (h5class == H5T_STRING) {
		if (as_int)
			warning("'as.integer' is ignored when "
				"dataset class is H5T_STRING");
		*Rtype = STRSXP;
		if (storage_mode == NULL
		 || strcmp(storage_mode, "character") == 0)
			return 0;
		warning("dataset class is H5T_STRING but it has a "
			"\"storage.mode\" attribute set to %s -- we ignore it",
			storage_mode);
		return 0;
	}
	ret = map_numeric_h5class_to_Rtype(h5class, h5type_size, Rtype);
	if (ret < 0)
		return -1;
	if (as_int) {
		/* We override the *Rtype value set earlier by
		   map_numeric_h5class_to_Rtype() but calling
		   map_numeric_h5class_to_Rtype() still had the merit
		   of checking that 'h5class' is a supported class. */
		*Rtype = INTSXP;
		return 0;
	}
	if (storage_mode == NULL)
		return 0;
	return map_storage_mode_to_Rtype(storage_mode, Rtype);
}

static size_t get_Rtype_size(SEXPTYPE Rtype, size_t h5size)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return sizeof(int);
	    case REALSXP:             return sizeof(double);
	    case STRSXP:              return h5size;
	    case RAWSXP:              return sizeof(char);
	}
	PRINT_TO_ERRMSG_BUF("unsupported Rtype: %s", CHAR(type2str(Rtype)));
	return 0;
}

/* Must NOT be called on STRSXP! */
static hid_t map_Rtype_to_native_type(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return H5T_NATIVE_INT;
	    case REALSXP:             return H5T_NATIVE_DOUBLE;
	    case RAWSXP:              return H5T_NATIVE_UCHAR;
	}
	/* Should never happen. */
	PRINT_TO_ERRMSG_BUF("failed to map Rtype %s to a native type",
			    CHAR(type2str(Rtype)));
	return -1;
}


/****************************************************************************
 * _init_H5DSetDescriptor() / _destroy_H5DSetDescriptor()
 */

void _destroy_H5DSetDescriptor(H5DSetDescriptor *h5dset)
{
	if (h5dset->h5nchunk != NULL)
		free(h5dset->h5nchunk);
	if (h5dset->h5chunkdim != NULL &&
	    h5dset->h5chunkdim != h5dset->h5dim)
		free(h5dset->h5chunkdim);
	if (h5dset->h5dim != NULL)
		free(h5dset->h5dim);
	if (h5dset->h5plist_id != -1)
		H5Pclose(h5dset->h5plist_id);
	if (h5dset->h5space_id != -1)
		H5Sclose(h5dset->h5space_id);
	if (h5dset->h5type_id != -1)
		H5Tclose(h5dset->h5type_id);
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
	hid_t h5type_id, native_type_id, native_type_id_for_Rtype,
	      h5space_id, h5plist_id, tmp;
	H5T_class_t h5class;
	size_t h5type_size, Rtype_size, native_type_size;
	SEXPTYPE Rtype;
	int as_na_attr, ndim, *h5nchunk, h5along;
	hsize_t *h5dim, *h5chunkdim, d, chunkd, nchunk;
	H5D_layout_t h5layout;
	htri_t ret;
	CharAE *buf;

	h5dset->dset_id = dset_id;

	/* Initialize the members that _destroy_H5DSetDescriptor() will free
	   or close. */
	h5dset->h5name = NULL;
	h5dset->storage_mode_attr = NULL;
	h5dset->h5type_id = -1;
	h5dset->h5space_id = -1;
	h5dset->h5plist_id = -1;
	h5dset->h5dim = NULL;
	h5dset->h5chunkdim = NULL;
	h5dset->h5nchunk = NULL;

	/* Set member 'h5name'. */
	h5name = get_h5name(dset_id);
	if (h5name == NULL)
		goto on_error;
	h5dset->h5name = h5name;

	/* Set member 'storage_mode_attr'. */
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

	/* Set member 'h5type_id'. */
	h5type_id = H5Dget_type(dset_id);
	if (h5type_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_type() returned an error");
		goto on_error;
	}
	h5dset->h5type_id = h5type_id;

	/* Set member 'h5class'. */
	h5class = H5Tget_class(h5dset->h5type_id);
	if (h5class == H5T_NO_CLASS) {
		PRINT_TO_ERRMSG_BUF("H5Tget_class() returned an error");
		goto on_error;
	}
	if (h5class == H5T_STRING &&
	    H5Tis_variable_str(h5dset->h5type_id) != 0) {
		PRINT_TO_ERRMSG_BUF("reading variable-length string data "
				    "is not supported at the moment");
		goto on_error;
	}
	h5dset->h5class = h5class;

	/* Set member 'h5type_size'. */
	h5type_size = H5Tget_size(h5dset->h5type_id);
	if (h5type_size == 0) {
		PRINT_TO_ERRMSG_BUF("H5Tget_size(h5type_id) returned 0");
		goto on_error;
	}
	h5dset->h5type_size = h5type_size;

	/* Set member 'Rtype'. */
	ret = set_Rtype(h5dset->h5class, h5dset->h5type_size,
			as_int, h5dset->storage_mode_attr, &Rtype);
	if (ret < 0)
		goto on_error;
	h5dset->Rtype = Rtype;

	if (get_Rtype_only)
		return 0;

	/* Set member 'Rtype_size'. */
	Rtype_size = get_Rtype_size(h5dset->Rtype, h5dset->h5type_size);
	if (Rtype_size == 0)
		goto on_error;
	h5dset->Rtype_size = Rtype_size;

	if (h5dset->h5class != H5T_STRING) {
		/* Set member 'native_type_id'.
		   Gosh, it really sucks that H5Tget_native_type() doesn't
		   actually return the id of a predefined native type like
		   H5T_NATIVE_INT. Instead it returns a stupid random type
		   id that cannot be compared directly with the predefined
		   ones using something like the == operator. Nope, you need
		   to compare with H5Tequal()! And to make the matter even
		   worse, the documentation for H5Tget_native_type() doesn't
		   say anything about this! */
		tmp = H5Tget_native_type(h5dset->h5type_id, H5T_DIR_DEFAULT);
		if (tmp < 0) {
			PRINT_TO_ERRMSG_BUF("H5Tget_native_type() "
					    "returned an error");
			goto on_error;
		}
		native_type_id = map_native_type_to_predef_type(tmp);
		H5Tclose(tmp);
		if (native_type_id < 0)
			goto on_error;
		h5dset->native_type_id = native_type_id;

		/* Set member 'native_type_size'. */
		native_type_size = H5Tget_size(h5dset->native_type_id);
		if (native_type_size == 0) {
			PRINT_TO_ERRMSG_BUF("H5Tget_size(native_type_id) "
					    "returned 0");
			goto on_error;
		}
		h5dset->native_type_size = native_type_size;

		/* Set member 'native_type_id_for_Rtype'. */
		native_type_id_for_Rtype =
			map_Rtype_to_native_type(h5dset->Rtype);
		if (native_type_id_for_Rtype < 0)
			goto on_error;
		h5dset->native_type_id_for_Rtype = native_type_id_for_Rtype;
	}

	/* Set member 'as_na_attr'. */
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

	/* Set member 'h5space_id'. */
	h5space_id = H5Dget_space(dset_id);
	if (h5space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_space() returned an error");
		goto on_error;
	}
	h5dset->h5space_id = h5space_id;

	/* Set member 'ndim'. */
	ndim = H5Sget_simple_extent_ndims(h5dset->h5space_id);
	if (ndim < 0) {
		PRINT_TO_ERRMSG_BUF(
			"H5Sget_simple_extent_ndims() returned an error");
		goto on_error;
	}
	h5dset->ndim = ndim;

	/* Set member 'h5plist_id'. */
	h5plist_id = H5Dget_create_plist(dset_id);
	if (h5plist_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_create_plist() returned an error");
		goto on_error;
	}
	h5dset->h5plist_id = h5plist_id;

	/* Set member 'h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		goto on_error;
	if (H5Sget_simple_extent_dims(h5space_id, h5dim, NULL) != ndim) {
		PRINT_TO_ERRMSG_BUF("H5Sget_simple_extent_dims() returned "
				    "an unexpected value");
		goto on_error;
	}
	h5dset->h5dim = h5dim;

	/* Set member 'h5layout'. */
	h5layout = H5Pget_layout(h5dset->h5plist_id);
	if (h5layout < 0) {
		PRINT_TO_ERRMSG_BUF("H5Pget_layout() returned an error");
		goto on_error;
	}
	h5dset->h5layout = h5layout;

	/* Set member 'h5chunkdim'. */
	if (h5dset->h5layout == H5D_CHUNKED) {
		h5chunkdim = _alloc_hsize_t_buf(ndim, 0, "'h5chunkdim'");
		if (h5chunkdim == NULL)
			goto on_error;
		if (H5Pget_chunk(h5plist_id, ndim, h5chunkdim) != ndim) {
			PRINT_TO_ERRMSG_BUF("H5Pget_chunk() returned "
					    "an unexpected value");
			goto on_error;
		}
		h5dset->h5chunkdim = h5chunkdim;
	} else if (h5dset->Rtype == STRSXP) {
		/* Even though the dataset is contiguous, we treat it as
		   if it was made of a single chunk. This is so that we
		   can use h5mread() method 4 on it, which works only on
		   chunked data and is the only method that knows how to
		   handle string data. */
		h5dset->h5chunkdim = h5dset->h5dim;
	}

	/* Set member 'h5nchunk'. */
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
		error(_HDF5Array_global_errmsg_buf());
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
	Rprintf("- dset_id = %ld\n", h5dset->dset_id);

	Rprintf("- h5name = \"%s\"\n", h5dset->h5name);

	Rprintf("- storage_mode_attr = ");
	if (h5dset->storage_mode_attr == NULL) {
		Rprintf("NULL");
	} else {
		Rprintf("\"%s\"", h5dset->storage_mode_attr);
	}
	Rprintf("\n");

	Rprintf("- h5type_id = %ld\n", h5dset->h5type_id);

	Rprintf("- h5class = %s\n", H5class2str(h5dset->h5class));

	Rprintf("- h5type_size = %lu\n", h5dset->h5type_size);

	Rprintf("- Rtype = \"%s\"\n", CHAR(type2str(h5dset->Rtype)));

	Rprintf("- Rtype_size = %lu\n", h5dset->Rtype_size);

	if (h5dset->h5class != H5T_STRING) {
		Rprintf("- native_type_id = %s\n",
			predef_native_type_as_string(h5dset->native_type_id));
		Rprintf("- native_type_size = %lu\n", h5dset->native_type_size);
		Rprintf("- native_type_id_for_Rtype = %s\n",
			predef_native_type_as_string(
				h5dset->native_type_id_for_Rtype));
	} else {
		Rprintf("- native_type_id = %s\n", "none          "
			" (not set when h5class is H5T_STRING)");
		Rprintf("- native_type_size = %s\n", "none        "
			" (not set when h5class is H5T_STRING)");
		Rprintf("- native_type_id_for_Rtype = %s\n", "none"
			" (not set when h5class is H5T_STRING)");
	}

	Rprintf("- as_na_attr = %d\n", h5dset->as_na_attr);

	Rprintf("- h5space_id = %ld\n", h5dset->h5space_id);

	Rprintf("- ndim = %d\n", h5dset->ndim);

	Rprintf("- h5plist_id = %ld\n", h5dset->h5plist_id);

	Rprintf("- h5dim =");
	for (h5along = 0; h5along < h5dset->ndim; h5along++)
		Rprintf(" %llu", h5dset->h5dim[h5along]);
	Rprintf("\n");

	Rprintf("- h5layout = %s\n", H5layout2str(h5dset->h5layout));

	Rprintf("- h5chunkdim =");
	if (h5dset->h5chunkdim == NULL) {
		Rprintf(" NULL\n");
	} else {
		for (h5along = 0; h5along < h5dset->ndim; h5along++)
			Rprintf(" %llu",
				h5dset->h5chunkdim[h5along]);
		if (h5dset->h5layout != H5D_CHUNKED &&
		    h5dset->h5chunkdim == h5dset->h5dim)
			Rprintf(" (artificially set to h5dim)");
		Rprintf("\n");
		Rprintf("    h5nchunk =");
		for (h5along = 0; h5along < h5dset->ndim; h5along++)
			Rprintf(" %d", h5dset->h5nchunk[h5along]);
		Rprintf("\n");
	}

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
		error(_HDF5Array_global_errmsg_buf());
	_destroy_H5DSetDescriptor(&h5dset);
	return ScalarString(type2str(h5dset.Rtype));
}

