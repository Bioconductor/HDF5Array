#ifndef _H5DSETDESCRIPTOR_H_
#define _H5DSETDESCRIPTOR_H_

#include <Rdefines.h>
#include "S4Vectors_interface.h"
#include "hdf5.h"

/* A data structure for describing an H5 type. */
typedef struct h5type_descriptor {
	/* Core struct members (always set). */
	hid_t h5type_id;
	H5T_class_t h5class;
	size_t h5type_size;
	int Rtype_is_set;
	int num_h5tmembers;
	struct h5tmember_descriptor **h5tmembers;

	/* Struct members below will be set only if h5class is **not**
	   H5T_COMPOUND or H5T_ENUM. */
	int signedness;			// set only if h5class == H5T_INTEGER
	SEXPTYPE Rtype;

	/* Struct members below will be set only if h5class is **not**
	   H5T_COMPOUND or H5T_ENUM, and when new_H5TypeDescriptor() is
	   called with 'get_Rtype_only' set to 0. */

	int is_variable_str;		// set to 0 if h5class != H5T_STRING

	size_t Rtype_size;		// set only if Rtype is set and
					// is_variable_str == 0

	/* Struct members below will be set only if h5class is H5T_INTEGER
	   or H5T_FLOAT, and when new_H5TypeDescriptor() is called
	   with 'get_Rtype_only' set to 0. */
	hid_t native_type_id;
	size_t native_type_size;
	hid_t native_type_id_for_Rtype; // set only if Rtype is set
} H5TypeDescriptor;

typedef struct h5tmember_descriptor {
	char *name;

	/* Struct members below will be set for a compound H5 type member, but
	   **not** for an enumeration H5 type member. */
	H5T_class_t h5class;
	struct h5type_descriptor *h5type;
} H5TMemberDescriptor;

/* A data structure for handling an HDF5 dataset. Collect various information
   about the dataset. What is collected is basically the union of the things
   needed by functions C_h5mread(), C_h5getdimscales(), and C_h5setdimscales().
 */
typedef struct h5dset_descriptor_t {
	hid_t dset_id;
	char *h5name;  // canonical name as retrieved by H5Iget_name()
	char *storage_mode_attr;
	struct h5type_descriptor *h5type;

	/* Struct members below will be set only when _init_H5DSetDescriptor()
	   is called with 'get_Rtype_only' set to 0. */
	int as_na_attr;
	hid_t h5space_id;
	int ndim;
	hid_t h5plist_id;
	hsize_t *h5dim;
	H5D_layout_t h5layout;
	hsize_t *h5chunkdim;
	int *h5nchunk;
} H5DSetDescriptor;


hsize_t *_alloc_hsize_t_buf(
	size_t buflength,
	int zeroes,
	const char *what
);

int _get_h5attrib_strval(
	hid_t dset_id,
	const char *attr_name,
	CharAE *buf
);

void _destroy_H5DSetDescriptor(
	H5DSetDescriptor *h5dset
);

int _init_H5DSetDescriptor(
	H5DSetDescriptor *h5dset,
	hid_t dset_id,
	int as_int,
	int Rtype_only
);

hid_t _get_file_id(
	SEXP filepath,
	int readonly
);

hid_t _get_dset_id(
	hid_t file_id,
	SEXP name,
	SEXP filepath
);

SEXP C_destroy_H5DSetDescriptor_xp(SEXP xp);

SEXP C_new_H5DSetDescriptor_xp(
	SEXP filepath,
	SEXP name,
	SEXP as_integer
);

SEXP C_show_H5DSetDescriptor_xp(SEXP xp);

SEXP C_get_h5mread_returned_type(
	SEXP filepath,
	SEXP name,
	SEXP as_integer
);

#endif  /* _H5DSETDESCRIPTOR_H_ */

