/****************************************************************************
 *         Low-level manipulation of HDF5 Dimension Scale datasets          *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include "hdf5_hl.h"

#include <stdlib.h>  /* for malloc, free */
#include <limits.h>  /* for INT_MAX */


static int get_dsname(hid_t dset_id, CharAE *name_buf)
{
	ssize_t name_size;

	name_size = H5Iget_name(dset_id, NULL, 0);
	if (name_size < 0) {
		PRINT_TO_ERRMSG_BUF("H5Iget_name() returned an error");
		return -1;
	}
	name_size++;
	if ((size_t) name_size > name_buf->_buflength)
		CharAE_extend(name_buf, (size_t) name_size);
	CharAE_set_nelt(name_buf, (size_t) name_size);
	name_size = H5Iget_name(dset_id, name_buf->elts, (size_t) name_size);
	if (name_size < 0) {
		PRINT_TO_ERRMSG_BUF("H5Iget_name() returned an error");
		return -1;
	}
	return 0;
}

typedef struct {
	const char *dset_name;
	const char *dsset_name;
	const char *scalename;
	CharAE *name_buf;
	CharAE *NAME_buf;
} VisitorData;

static herr_t scale_visitor1(hid_t dset_id, unsigned int along, hid_t dsset_id,
			     void *data)
{
	VisitorData *visitor_data;
	int ret;

	visitor_data = (VisitorData *) data;
	ret = get_dsname(dsset_id, visitor_data->name_buf);
	if (ret < 0)
		return -1;
	ret = _get_h5_attrib_str(dsset_id, "NAME", visitor_data->NAME_buf);
	if (ret < 0)
		return -1;
	if (ret == 0) {
		if (visitor_data->scalename != NULL)
			return 0;
		return 1;
	}
	if (visitor_data->scalename == NULL ||
	    ret == 1 ||
	    strcmp(visitor_data->NAME_buf->elts, visitor_data->scalename) != 0)
		return 0;
	return 1;
}

static herr_t scale_visitor2(hid_t dset_id, unsigned int along, hid_t dsset_id,
			     void *data)
{
	VisitorData *visitor_data;
	int ret;

	visitor_data = (VisitorData *) data;
	ret = get_dsname(dsset_id, visitor_data->name_buf);
	if (ret < 0)
		return -1;
	if (strcmp(visitor_data->name_buf->elts, visitor_data->dsset_name) == 0)
		return 0;
	ret = _get_h5_attrib_str(dsset_id, "NAME", visitor_data->NAME_buf);
	if (ret < 0)
		return -1;
	if (ret == 0) {
		if (visitor_data->scalename != NULL)
			return 0;
		PRINT_TO_ERRMSG_BUF("dimension %u of dataset '%s' is already "
				    "linked to unnamed\n  Dimension Scale "
				    "dataset '%s'", along + 1,
				    visitor_data->dset_name,
				    visitor_data->name_buf->elts);
		return -1;
	}
	if (visitor_data->scalename == NULL ||
	    ret == 1 ||
	    strcmp(visitor_data->NAME_buf->elts, visitor_data->scalename) != 0)
		return 0;
	PRINT_TO_ERRMSG_BUF("dimension %u of dataset '%s' is already "
			    "linked to Dimension Scale\n  dataset '%s' "
			    "named \"%s\" (via \"NAME\" attribute)", along + 1,
			    visitor_data->dset_name,
			    visitor_data->name_buf->elts,
			    visitor_data->scalename);
	return -1;
}

static int get_scale_along(hid_t dset_id, int along, const char *scalename,
		CharAE *name_buf, CharAE *NAME_buf)
{
	int *idx;
	VisitorData visitor_data;

	visitor_data.scalename = scalename;
	visitor_data.name_buf = name_buf;
	visitor_data.NAME_buf = NAME_buf;

	idx = NULL;
	return H5DSiterate_scales(dset_id, (unsigned int) along, idx,
				  scale_visitor1, &visitor_data);
}

static int check_no_other_links_with_this_NAME(
		hid_t dset_id, const char *dset_name,
		int along, hid_t dsset_id, const char *scalename,
		CharAE *buf1, CharAE *buf2, CharAE *buf3)
{
	int ret, *idx;
	VisitorData visitor_data;

	ret = get_dsname(dsset_id, buf1);
	if (ret < 0)
		return -1;
	visitor_data.dset_name = dset_name;
	visitor_data.dsset_name = buf1->elts;
	visitor_data.scalename = scalename;
	visitor_data.name_buf = buf2;
	visitor_data.NAME_buf = buf3;

	idx = NULL;
	return H5DSiterate_scales(dset_id, (unsigned int) along, idx,
				  scale_visitor2, &visitor_data);
}


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
				    "Dimension Scale named \"%s\":\n  "
				    "dataset is already a Dimension Scale "
				    "and is an **unnamed** one", dsset_name,
				    scalename);
		return -1;
	}
	if (scalename == NULL) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to "
				    "an unnamed Dimension Scale:\n  "
				    "there is a \"NAME\" attribute "
				    "on this dataset", dsset_name);
		return -1;
	}
	if (ret == 1) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "Dimension Scale named \"%s\":\n  "
				    "there is already a \"NAME\" attribute "
				    "on this dataset but it's not of class "
				    "H5T_STRING\n  (converting the dataset "
				    "to a Dimension Scale would replace this "
				    "attribute)", dsset_name, scalename);
		return -1;
	}
	if (strcmp(NAME_buf->elts, scalename) != 0) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "Dimension Scale named \"%s\":\n  "
				    "there is already a \"NAME\" attribute "
				    "set to \"%s\" on this dataset",
				    dsset_name, scalename, NAME_buf->elts);
		return -1;
	}
	return 0;
}

/* If 'buf1' != NULL then only check feasibility but don't do it (dry-run
   mode). */
static int set_scale_along(hid_t file_id,
		hid_t dset_id, const char *dset_name,
		int along, const char *dsset_name, const char *scalename,
		CharAE *buf1, CharAE *buf2, CharAE *buf3)
{
	hid_t dsset_id;
	int ret, is_scale, is_attached;

	dsset_id = H5Dopen(file_id, dsset_name, H5P_DEFAULT);
	if (dsset_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dsset_name);
		return -1;
	}
	if (buf1 != NULL) {
		/* Retrieve **canonical name** of dataset 'dsset_id'. */
		ret = get_dsname(dsset_id, buf1);
		if (ret < 0) {
			H5Dclose(dsset_id);
			return -1;
		}
		if (strcmp(buf1->elts, dset_name) == 0) {
			H5Dclose(dsset_id);
			PRINT_TO_ERRMSG_BUF("cannot set dataset '%s' as a "
					    "Dimension Scale on itself",
					    dset_name);
			return -1;
		}
	}
	is_scale = H5DSis_scale(dsset_id);
	if (is_scale < 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("H5DSis_scale() returned an error");
		return -1;
	}
	if (is_scale) {
		is_attached = H5DSis_attached(dset_id, dsset_id,
					      (unsigned int) along);
		if (is_attached < 0) {
			H5Dclose(dsset_id);
			PRINT_TO_ERRMSG_BUF("H5DSis_attached() "
					    "returned an error");
			return -1;
		}
	} else {
		is_attached = 0;
	}
	if (buf1 != NULL) {
		ret = check_NAME_attribute(dsset_id, dsset_name, is_scale,
					   scalename, buf1);
		if (ret < 0) {
			H5Dclose(dsset_id);
			return -1;
		}
		ret = check_no_other_links_with_this_NAME(dset_id, dset_name,
							  along, dsset_id,
							  scalename,
							  buf1, buf2, buf3);
		H5Dclose(dsset_id);
		return ret;
	}
	if (!is_scale) {
		ret = H5DSset_scale(dsset_id, scalename);
		if (ret < 0) {
			H5Dclose(dsset_id);
			PRINT_TO_ERRMSG_BUF("H5DSset_scale() "
					    "returned an error");
			return -1;
		}
	}
	if (!is_attached) {
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
	CharAE *buf0, *buf1, *buf2, *buf3;
	int ret, ndim, along;
	const char *scalename0;
	SEXP dsset_name;

	file_id = _get_file_id(filepath, 0);
	dset_id = _get_dset_id(file_id, name, filepath);

	/* Retrieve **canonical name** of dataset 'dset_id' (we cannot use
	   'CHAR(STRING_ELT(name, 0))' for that). */
	buf0 = new_CharAE(0);
	ret = get_dsname(dset_id, buf0);
	if (ret < 0)
		goto on_error;

	ret = H5DSis_scale(dset_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5DSis_scale() returned an error");
		goto on_error;
	}
	if (ret > 0) {
		PRINT_TO_ERRMSG_BUF("dataset '%s' is a Dimension Scale "
				    "(a Dimension Scale dataset\n  cannot "
				    "have Dimension Scales)", buf0->elts);
		ret = -1;
		goto on_error;
	}

	/* We treat an empty scalename as no scalename at all. */
	if (scalename == R_NilValue || LENGTH(STRING_ELT(scalename, 0)) == 0) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}
	ndim = LENGTH(dsnames);

	/* We use 2 passes to make the change atomic. */

	/* 1st pass: dry-run */
	buf1 = new_CharAE(0);
	buf2 = new_CharAE(0);
	buf3 = new_CharAE(0);
	for (along = 0; along < ndim; along++) {
		dsset_name = STRING_ELT(dsnames, along);
		if (dsset_name != NA_STRING) {
			ret = set_scale_along(file_id, dset_id, buf0->elts,
					      along, CHAR(dsset_name),
					      scalename0, buf1, buf2, buf3);
			if (ret < 0)
				goto on_error;
		}
	}

	/* 2nd pass: do it */
	for (along = 0; along < ndim; along++) {
		dsset_name = STRING_ELT(dsnames, along);
		if (dsset_name != NA_STRING) {
			ret = set_scale_along(file_id, dset_id, buf0->elts,
					      along, CHAR(dsset_name),
					      scalename0, NULL, NULL, NULL);
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

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimscales(SEXP filepath, SEXP name, SEXP scalename)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	const char *scalename0;
	CharAE *name_buf, *NAME_buf;
	int along, ret;

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

	name_buf = new_CharAE(0);
	NAME_buf = new_CharAE(0);
	for (along = 0; along < dset_handle.ndim; along++) {
		printf("along = %d: ", along);
		ret = get_scale_along(dset_id, along, scalename0,
				      name_buf, NAME_buf);
		if (ret < 0) {
			_close_DSetHandle(&dset_handle);
			H5Fclose(file_id);
			error(_HDF5Array_errmsg_buf());
		}
		if (ret == 0) {
			printf("NA\n");
		} else {
			printf("%s\n", name_buf->elts);
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

