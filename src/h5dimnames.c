/****************************************************************************
 *           Set/find the HDF5 datasets representing the dimnames           *
 *                         of a given HDF5 dataset                          *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include "hdf5_hl.h"

#include <stdlib.h>  /* for malloc, free */
#include <limits.h>  /* for INT_MAX */


/****************************************************************************
 * C_h5setdimnames() and C_h5getdimnames()
 */

static int set_names_along(hid_t file_id, hid_t dset_id, int along,
			   const char *dn)
{
	hid_t dsset_id;
	int ret;

	dsset_id = H5Dopen(file_id, dn, H5P_DEFAULT);
	if (dsset_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dn);
		return -1;
	}
	ret = H5DSset_scale(dsset_id, "names");
	if (ret < 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("H5DSset_scale() failed");
		return -1;
	}
	ret = H5DSattach_scale(dset_id, dsset_id, (unsigned int) along);
	if (ret < 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("H5DSattach_scale() failed");
		return -1;
	}
	H5Dclose(dsset_id);
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimnames(SEXP filepath, SEXP name, SEXP dimnames)
{
	hid_t file_id, dset_id;
	int ndim, along, ret;
	SEXP dn;

	file_id = _get_file_id(filepath, 0);
	dset_id = _get_dset_id(file_id, name, filepath);

	ndim = LENGTH(dimnames);
	for (along = 0; along < ndim; along++) {
		dn = STRING_ELT(dimnames, along);
		if (dn != NA_STRING) {
			ret = set_names_along(file_id, dset_id, along,
					      CHAR(dn));
			if (ret < 0) {
				H5Dclose(dset_id);
				H5Fclose(file_id);
				error(_HDF5Array_errmsg_buf());
			}
		}
	}

	H5Dclose(dset_id);
	H5Fclose(file_id);
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimnames(SEXP filepath, SEXP name)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	int along, ret;

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, 0, 0, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	for (along = 0; along < dset_handle.ndim; along++) {
		ret = H5DSget_num_scales(dset_id, (unsigned int) along);
		printf("num_scales = %d\n", ret);
	}

	_close_DSetHandle(&dset_handle);
	H5Fclose(file_id);
	return R_NilValue;
}


/****************************************************************************
 * C_h5setdimlabels() and C_h5getdimlabels()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimlabels(SEXP filepath, SEXP name, SEXP labels)
{
	hid_t file_id, dset_id;
	int ndim, along, ret;
	SEXP labels_elt;
	const char *label;

	if (labels == R_NilValue)
		return R_NilValue;

	file_id = _get_file_id(filepath, 0);
	dset_id = _get_dset_id(file_id, name, filepath);

	ndim = LENGTH(labels);
	for (along = 0; along < ndim; along++) {
		labels_elt = STRING_ELT(labels, along);
		if (labels_elt == NA_STRING || LENGTH(labels_elt) == 0) {
			label = "";
		} else {
			label = CHAR(labels_elt);
		}
		ret = H5DSset_label(dset_id, (unsigned int) along, label);
		if (ret < 0) {
			H5Dclose(dset_id);
			H5Fclose(file_id);
			PRINT_TO_ERRMSG_BUF("H5DSset_label() failed on "
					    "label %d", along + 1);
			error(_HDF5Array_errmsg_buf());
		}
	}

	H5Dclose(dset_id);
	H5Fclose(file_id);
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimlabels(SEXP filepath, SEXP name)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	int along;
	ssize_t max_label_size, label_size;
	char *label_buf;
	SEXP labels, label;

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, 0, 0, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	/* First pass */
	max_label_size = 0;
	for (along = 0; along < dset_handle.ndim; along++) {
		label_size = H5DSget_label(dset_id, (unsigned int) along,
					   NULL, 0);
		if (label_size < 0) {
			_close_DSetHandle(&dset_handle);
			H5Fclose(file_id);
			error("H5DSget_label() failed");
		}
		//printf("label_size = %ld\n", label_size);
		if (label_size > max_label_size)
			max_label_size = label_size;
	}

	if (max_label_size == 0) {
		_close_DSetHandle(&dset_handle);
		H5Fclose(file_id);
		return R_NilValue;
	}

	/* Second pass */
	if (max_label_size > INT_MAX) {
		max_label_size = INT_MAX;
		warning("some dimension labels were too big "
			"so have been truncated");
	}
	label_buf = (char *) malloc((size_t) max_label_size + 1);
	if (label_buf == NULL) {
		_close_DSetHandle(&dset_handle);
		H5Fclose(file_id);
		error("failed to allocate memory for 'label_buf'");
	}
	labels = PROTECT(NEW_CHARACTER(dset_handle.ndim));
	for (along = 0; along < dset_handle.ndim; along++) {
		label_size = H5DSget_label(dset_id, (unsigned int) along,
					   label_buf, max_label_size + 1);
		/* Should never happen. */
		if (label_size < 0) {
			free(label_buf);
			_close_DSetHandle(&dset_handle);
			H5Fclose(file_id);
			error("H5DSget_label() failed");
		}
		if (label_size > INT_MAX)
			label_size = INT_MAX;
		label = PROTECT(mkCharLen(label_buf, (int) label_size));
		SET_STRING_ELT(labels, along, label);
		UNPROTECT(1);
	}

	free(label_buf);
	_close_DSetHandle(&dset_handle);
	H5Fclose(file_id);
	UNPROTECT(1);
	return labels;
}

