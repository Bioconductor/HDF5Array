#ifndef _H5DSETDESCRIPTOR_H_
#define _H5DSETDESCRIPTOR_H_

#include <Rdefines.h>
#include "S4Vectors_interface.h"
#include "hdf5.h"

/* A data structure for handling an HDF5 dataset. Collect various information
   about the dataset. What is collected is basically the union of the things
   needed by functions C_h5mread(), C_h5getdimscales(), and C_h5setdimscales().
 */
typedef struct {
	hid_t dset_id, dtype_id, space_id, plist_id, mem_type_id;
	char *h5name;  /* canonical name as retrieved by H5Iget_name() */
	char *storage_mode_attr;
	H5T_class_t H5class;
	size_t H5size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE Rtype;
	int as_na_attr, ndim, *h5nchunk;
	hsize_t *h5dim, *h5chunkdim;
	H5D_layout_t H5layout;
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

