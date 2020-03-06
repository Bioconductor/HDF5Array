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

static herr_t scale_visitor(hid_t dset_id, unsigned int along, hid_t dsset_id,
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
				  scale_visitor, &visitor_data);
}


/****************************************************************************
 * C_h5setdimscales() and C_h5getdimscales()
 */

static int check_NAME_attribute(hid_t dsset_id, const char *dsname,
				const char *scalename, CharAE *NAME_buf)
{
	int ret;

	ret = _get_h5_attrib_str(dsset_id, "NAME", NAME_buf);
	if (ret < 0)
		return -1;
	if (ret == 0) {
		if (scalename == NULL)
			return 0;
		ret = H5DSis_scale(dsset_id);
		if (ret < 0) {
			PRINT_TO_ERRMSG_BUF("H5DSis_scale() returned an error");
			return -1;
		}
		if (ret > 0) {
			PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
					    "Dimension Scale named \"%s\"\n  "
					    "because dataset is already a "
					    "Dimension Scale but is "
					    "an **unnamed** one",
					    dsname, scalename);
			return -1;
		}
		return 0;
	}
	if (scalename == NULL) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to an "
				    "**unnamed** Dimension Scale\n  "
				    "because dataset has a \"NAME\" "
				    "attribute", dsname);
		return -1;
	}
	if (ret == 1) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "Dimension Scale named \"%s\"\n  "
				    "because dataset already has a \"NAME\" "
				    "attribute that is **not** of class "
				    "H5T_STRING\n  (converting the dataset "
				    "to a Dimension Scale would replace this "
				    "attribute)", dsname, scalename);
		return -1;
	}
	if (strcmp(NAME_buf->elts, scalename) != 0) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "Dimension Scale named \"%s\"\n  "
				    "because dataset already has a \"NAME\" "
				    "attribute set to \"%s\"",
				    dsname, scalename, NAME_buf->elts);
		return -1;
	}
	return 0;
}

/* Check whether or not dataset 'dsname' is already set as a Dimension Scale
   named 'scalename' (via its "NAME" attribute) on dimension 'along' of
   dataset 'dset_id'. Return: 0 if no, 1 if yes, -1 if error.
 */
static int check_scale_along(hid_t file_id,
		hid_t dset_id, const char *dset_name,
		int along, const char *dsname, const char *scalename,
		CharAE *name_buf1, CharAE *name_buf2, CharAE *buf)
{
	hid_t dsset_id;
	int ret;

	dsset_id = H5Dopen(file_id, dsname, H5P_DEFAULT);
	if (dsset_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dsname);
		return -1;
	}
	/* Retrieve **canonical name** of dataset 'dsset_id'. */
	ret = get_dsname(dsset_id, name_buf1);
	if (ret < 0) {
		H5Dclose(dsset_id);
		return -1;
	}
	if (strcmp(name_buf1->elts, dset_name) == 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("cannot set dataset '%s' as a "
				    "Dimension Scale on itself",
				    dset_name);
		return -1;
	}
	/* Check if 'dsset_id' has scales on it, in which case we won't be
	   able to turn it into a Dimension Scale. */
	ret = H5Aexists(dsset_id, "DIMENSION_LIST");
	if (ret < 0) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("H5Aexists() returned an error");
		return -1;
	}
	if (ret) {
		H5Dclose(dsset_id);
		PRINT_TO_ERRMSG_BUF("dataset '%s' already has Dimension "
				    "Scales on it so cannot be turned\n  "
				    "itself into a Dimension Scale (a "
                                    "Dimension Scale dataset cannot have\n  "
                                    "Dimension Scales)", name_buf1->elts);
		return -1;
	}
	ret = get_scale_along(dset_id, along, scalename,
			      name_buf2, buf);
	if (ret < 0) {
		H5Dclose(dsset_id);
		return -1;
	}
	if (ret == 0) {
		/* Now that we know that 'dsset_id' is not already set as
		   a Dimension Scale named 'scalename' on dimension 'along'
		   of 'dset_id', we also want to make sure that it **would**
		   be ok to create that link. We return -1 if not. */
		ret = check_NAME_attribute(dsset_id, dsname, scalename, buf);
		H5Dclose(dsset_id);
		return ret < 0 ? -1 : 0;
	}
	if (strcmp(name_buf2->elts, name_buf1->elts) == 0) {
		H5Dclose(dsset_id);
		return 1;
	}
	if (scalename == NULL) {
		PRINT_TO_ERRMSG_BUF("dimension %d of dataset '%s' is already "
				    "linked to unnamed\n  Dimension Scale "
				    "dataset '%s'", along + 1,
				    dset_name, name_buf2->elts);
	} else {
		PRINT_TO_ERRMSG_BUF("dimension %d of dataset '%s' is already "
				    "linked to Dimension Scale\n  "
				    "dataset '%s' named \"%s\" (via \"NAME\" "
				    "attribute)", along + 1,
				    dset_name, name_buf2->elts, scalename);
	}
	H5Dclose(dsset_id);
	return -1;
}

static int set_scale_along(hid_t file_id, hid_t dset_id,
		int along, const char *dsname, const char *scalename)
{
	hid_t dsset_id;
	int is_scale, is_attached, ret;

	dsset_id = H5Dopen(file_id, dsname, H5P_DEFAULT);
	if (dsset_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dsname);
		return -1;
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
		ret = H5DSset_scale(dsset_id, scalename);
		if (ret < 0) {
			H5Dclose(dsset_id);
			PRINT_TO_ERRMSG_BUF("H5DSset_scale() "
					    "returned an error");
			return -1;
		}
		is_attached = 0;
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

static SEXP check_scales(SEXP filepath, SEXP name, SEXP dsnames,
			 const char *scalename, CharAE *dset_name_buf)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	CharAE *buf1, *buf2, *buf3;
	int ret, along;
	SEXP ans, dsname;

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, 0, 0, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	/* Retrieve **canonical name** of dataset 'dset_id' ('name' contains
	   a user-supplied name that is not necessarily canonical). */
	ret = get_dsname(dset_id, dset_name_buf);
	if (ret < 0)
		goto on_error;

	if (LENGTH(dsnames) > dset_handle.ndim) {
		PRINT_TO_ERRMSG_BUF("'dsnames' cannot be longer than the "
				    "nb of dimensions of dataset '%s' (%d)",
				    dset_name_buf->elts, dset_handle.ndim);
		ret = -1;
		goto on_error;
	}

	ret = H5DSis_scale(dset_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5DSis_scale() returned an error");
		goto on_error;
	}
	if (ret > 0) {
		PRINT_TO_ERRMSG_BUF("dataset '%s' is a Dimension Scale "
				    "(cannot set Dimension Scales\n  "
				    "on a Dimension Scale dataset)",
				    dset_name_buf->elts);
		ret = -1;
		goto on_error;
	}

	ans = PROTECT(NEW_LOGICAL(dset_handle.ndim));
	buf1 = new_CharAE(0);
	buf2 = new_CharAE(0);
	buf3 = new_CharAE(0);
	for (along = 0; along < LENGTH(dsnames); along++) {
		dsname = STRING_ELT(dsnames, along);
		if (dsname == NA_STRING) {
			LOGICAL(ans)[along] = 0;
			continue;
		}
		ret = check_scale_along(file_id,
					dset_id, dset_name_buf->elts,
					along, CHAR(dsname), scalename,
					buf1, buf2, buf3);
		if (ret < 0)
			goto on_error;
		LOGICAL(ans)[along] = !ret;
	}

    on_error:
	_close_DSetHandle(&dset_handle);
	H5Fclose(file_id);
	UNPROTECT(1);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf());
	return ans;
}

static void set_scales(SEXP filepath, SEXP name, SEXP dsnames,
		       const char *scalename, const char *dset_name)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	int ret, along;
	SEXP dsname;

	file_id = _get_file_id(filepath, 0);  /* read/write */
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, 0, 0, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	for (along = 0; along < LENGTH(dsnames); along++) {
		dsname = STRING_ELT(dsnames, along);
		if (dsname == NA_STRING)
			continue;
		ret = set_scale_along(file_id, dset_id,
				      along, CHAR(dsname), scalename);
		if (ret < 0)  /* should never happen */
			goto on_error;
	}

    on_error:
	_close_DSetHandle(&dset_handle);
	H5Fclose(file_id);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf());
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimscales(SEXP filepath, SEXP name, SEXP dsnames, SEXP scalename,
		      SEXP dry_run)
{
	const char *scalename0;
	CharAE *dset_name_buf;
	SEXP ans;
	int along, do_it;

	if (!IS_CHARACTER(dsnames))
		error("'dsnames' must be a character vector");

	if (STRING_ELT(scalename, 0) == NA_STRING) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}

	dset_name_buf = new_CharAE(0);

	/* We use 2 passes to make the changes to the file atomic. */

	/* 1st pass: dry-run */
	ans = PROTECT(check_scales(filepath, name, dsnames, scalename0,
				   dset_name_buf));

	/* 2nd pass: do it */
	if (!LOGICAL(dry_run)[0]) {
		/* We could just call set_scales() and it will do the right
		   thing even if there is nothing to do (it acts as a no-op
		   along the dimensions that already have the requested scales).
		   Except that it starts by opening the file in read/write
		   mode so will fail on a read-only file!
		   By calling set_scales() only if there is something to do,
		   we make h5setdimscales() work on a read-only file if all
		   the requested scales are already set on dataset 'name'.
		   It also makes the thing twice faster which is nice. */
		do_it = 0;
		for (along = 0; along < LENGTH(dsnames); along++) {
			if (LOGICAL(ans)[along]) {
				do_it = 1;
				break;
			}
		}
		if (do_it)
			set_scales(filepath, name, dsnames, scalename0,
				   dset_name_buf->elts);
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimscales(SEXP filepath, SEXP name, SEXP scalename)
{
	hid_t file_id, dset_id;
	DSetHandle dset_handle;
	const char *scalename0;
	CharAE *name_buf, *NAME_buf;
	SEXP ans, ans_elt;
	int along, ret;

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, 0, 0, &dset_handle) < 0) {
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}
	if (STRING_ELT(scalename, 0) == NA_STRING) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}

	ans = PROTECT(NEW_CHARACTER(dset_handle.ndim));

	name_buf = new_CharAE(0);
	NAME_buf = new_CharAE(0);
	for (along = 0; along < dset_handle.ndim; along++) {
		ret = get_scale_along(dset_id, along, scalename0,
				      name_buf, NAME_buf);
		if (ret < 0) {
			_close_DSetHandle(&dset_handle);
			H5Fclose(file_id);
			error(_HDF5Array_errmsg_buf());
		}
		if (ret == 0) {
			SET_STRING_ELT(ans, along, NA_STRING);
		} else {
			ans_elt = PROTECT(mkChar(name_buf->elts));
			SET_STRING_ELT(ans, along, ans_elt);
			UNPROTECT(1);
		}
	}

	_close_DSetHandle(&dset_handle);
	H5Fclose(file_id);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_h5setdimlabels() and C_h5getdimlabels()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimlabels(SEXP filepath, SEXP name, SEXP labels)
{
	hid_t file_id, dset_id;
	int ndim, along, ret;
	SEXP label;

	if (labels == R_NilValue)
		return R_NilValue;

	file_id = _get_file_id(filepath, 0);
	dset_id = _get_dset_id(file_id, name, filepath);
	ndim = LENGTH(labels);

	for (along = 0; along < ndim; along++) {
		label = STRING_ELT(labels, along);
		if (label == NA_STRING)
			continue;
		ret = H5DSset_label(dset_id, (unsigned int) along, CHAR(label));
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

	file_id = _get_file_id(filepath, 1);  /* read-only */
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

