/****************************************************************************
 *         Low-level manipulation of HDF5 Dimension Scale datasets          *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include "hdf5_hl.h"

#include <stdlib.h>  /* for malloc, free */
#include <limits.h>  /* for INT_MAX */


/****************************************************************************
 * C_h5setdimscales() and C_h5getdimscales()
 */

static int check_NAME_attribute(hid_t dsset_id, const char *dsset_name,
				int is_scale,
				const char *scalename, CharAE *NAME_buf)
{
	int ret;

	ret = _get_h5_attrib_str(dsset_id, "NAME", NAME_buf);
	if (ret < 0)
		return -1;
	if (ret == 0) {
		if (!is_scale || scalename == NULL)
			return 0;
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "dimension scale named \"%s\":\n  "
				    "dataset is already a dimension scale "
				    "and is an **unnamed** one", dsset_name,
				    scalename);
		return -1;
	}
	if (scalename == NULL) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to "
				    "an unnamed dimension scale:\n  "
				    "there is a \"NAME\" attribute "
				    "on this dataset", dsset_name);
		return -1;
	}
	if (ret == 1) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "dimension scale named \"%s\":\n  "
				    "there is already a \"NAME\" attribute "
				    "on this dataset but it's not of class "
				    "H5T_STRING\n  (converting the dataset "
				    "to a dimension scale would replace this "
				    "attribute)", dsset_name, scalename);
		return -1;
	}
	if (strcmp(NAME_buf->elts, scalename) != 0) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "dimension scale named \"%s\":\n  "
				    "there is already a \"NAME\" attribute "
				    "set to \"%s\" on this dataset",
				    dsset_name, scalename, NAME_buf->elts);
		return -1;
	}
	return 0;
}

/* If 'NAME_buf' != NULL then only check feasibility but don't do it (dry-run
   mode). */
static int set_scale_along(hid_t file_id, hid_t dset_id, int along,
			   const char *dsset_name, const char *scalename,
			   CharAE *NAME_buf)
{
	hid_t dsset_id;
	int is_scale, ret;

	dsset_id = H5Dopen(file_id, dsset_name, H5P_DEFAULT);
	if (dsset_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dsset_name);
		return -1;
	}
	is_scale = H5DSis_scale(dsset_id);
	if (is_scale < 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("H5DSis_scale() returned an error");
		return -1;
	}
	if (NAME_buf != NULL) {
		ret = check_NAME_attribute(dsset_id, dsset_name, is_scale,
					   scalename, NAME_buf);
		H5Dclose(dsset_id);
		return ret;
	}
	if (is_scale == 0) {
		ret = H5DSset_scale(dsset_id, scalename);
		if (ret < 0) {
			H5Dclose(dsset_id);
			PRINT_TO_ERRMSG_BUF("H5DSset_scale() "
					    "returned an error");
			return -1;
		}
	}
	ret = H5DSis_attached(dset_id, dsset_id, (unsigned int) along);
	if (ret < 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("H5DSis_attached() returned an error");
		return -1;
	}
	if (ret == 0) {
		ret = H5DSattach_scale(dset_id, dsset_id, (unsigned int) along);
		if (ret < 0) {
			H5Dclose(dsset_id);
			PRINT_TO_ERRMSG_BUF("H5DSattach_scale() "
					    "returned an error");
			return -1;
		}
	}
	H5Dclose(dsset_id);
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimscales(SEXP filepath, SEXP name, SEXP scalename, SEXP dsnames)
{
	hid_t file_id, dset_id;
	const char *scalename0;
	int ndim, along, ret = 0;
	SEXP dsset_name;
	CharAE *NAME_buf;

	file_id = _get_file_id(filepath, 0);
	dset_id = _get_dset_id(file_id, name, filepath);
	/* We treat an empty scalename as no scalename at all. */
	if (scalename == R_NilValue || LENGTH(STRING_ELT(scalename, 0)) == 0) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}
	ndim = LENGTH(dsnames);

        /* We use 2 passes to make the change atomic. */

	/* 1st pass: dry-run */
	NAME_buf = new_CharAE(0);
	for (along = 0; along < ndim; along++) {
		dsset_name = STRING_ELT(dsnames, along);
		if (dsset_name != NA_STRING) {
			ret = set_scale_along(file_id, dset_id, along,
					      CHAR(dsset_name),
					      scalename0,
					      NAME_buf);
			if (ret < 0)
				goto on_error;
		}
	}

	/* 2nd pass: do it */
	for (along = 0; along < ndim; along++) {
		dsset_name = STRING_ELT(dsnames, along);
		if (dsset_name != NA_STRING) {
			ret = set_scale_along(file_id, dset_id, along,
					      CHAR(dsset_name),
					      scalename0,
					      NULL);
			if (ret < 0)  /* should never happen */
				goto on_error;
		}
	}

    on_error:
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf());
	return R_NilValue;
}

static herr_t scale_visitor1(hid_t dset_id, unsigned int along, hid_t dsset_id,
			     void *visitor_data)
{
	ssize_t max_name_size, name_size;
	char buffer[100];

	printf("scale_visitor1: along = %u -- dsset_id = %lu\n",
		along, dsset_id);
	name_size = H5Iget_name(dsset_id, NULL, 0);
	printf("name_size = %ld\n", name_size);
	H5Iget_name(dsset_id, buffer, name_size + 1);
	printf("name: %s\n", buffer);
	return 0;
}

static herr_t scale_visitor2(hid_t dset_id, unsigned int along, hid_t dsset_id,
			     void *visitor_data)
{
	const char *scalename;

	scalename = visitor_data;

	//printf("along = %u -- dsset_id = %lu\n", along, dsset_id);
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimscales(SEXP filepath, SEXP name, SEXP scalename)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	const char *scalename0;
	int along, ret, *idx;

	file_id = _get_file_id(filepath, 1);
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, 0, 0, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}
	/* We treat an empty scalename as no scalename at all. */
	if (scalename == R_NilValue || LENGTH(STRING_ELT(scalename, 0)) == 0) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}

	idx = NULL;
	for (along = 0; along < dset_handle.ndim; along++) {
		ret = H5DSget_num_scales(dset_id, (unsigned int) along);
		printf("num_scales = %d\n", ret);
		ret = H5DSiterate_scales(dset_id, (unsigned int) along,
				idx, scale_visitor1, NULL);
		if (ret < 0) {
			_close_DSetHandle(&dset_handle);
			H5Fclose(file_id);
			error(_HDF5Array_errmsg_buf());
		}
		ret = H5DSiterate_scales(dset_id, (unsigned int) along,
				idx, scale_visitor2, (void *) scalename0);
		if (ret < 0) {
			_close_DSetHandle(&dset_handle);
			H5Fclose(file_id);
			error(_HDF5Array_errmsg_buf());
		}
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
			error("H5DSget_label() returned an error");
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
			error("H5DSget_label() returned an error");
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

