/****************************************************************************
 *         Low-level manipulation of HDF5 Dimension Scale datasets          *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"

#include "hdf5_hl.h"

#include <stdlib.h>  /* for malloc, free */
#include <limits.h>  /* for INT_MAX */


/****************************************************************************
 * C_h5isdimscale()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_h5isdimscale(SEXP filepath, SEXP name)
{
	hid_t file_id, dset_id;
	int is_scale;

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	is_scale = H5DSis_scale(dset_id);
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (is_scale < 0)
		error("H5DSis_scale() returned an error");
	return ScalarLogical(is_scale);
}


/****************************************************************************
 * C_h5getdimscales()
 */

typedef struct {
	const char *scalename;
	H5DSetDescriptor *h5dimscale;
	CharAE *NAME_buf;
} VisitorData;

static herr_t visitor(hid_t dset_id, unsigned int along, hid_t dimscale_id,
		      void *data)
{
	VisitorData *visitor_data;
	H5DSetDescriptor *h5dimscale;
	int ret;

	visitor_data = (VisitorData *) data;
	h5dimscale = visitor_data->h5dimscale;
	ret = _init_H5DSetDescriptor(h5dimscale, dimscale_id, 0, 0);
	if (ret < 0)
		return -1;
	ret = _get_h5attrib_strval(dimscale_id, "NAME", visitor_data->NAME_buf);
	if (ret < 0) {
		_destroy_H5DSetDescriptor(h5dimscale);
		return -1;
	}
	if (ret == 0) {
		if (visitor_data->scalename == NULL)
			return 1;
		_destroy_H5DSetDescriptor(h5dimscale);
		return 0;
	}
	if (ret == 2 &&
	    visitor_data->scalename != NULL &&
	    strcmp(visitor_data->NAME_buf->elts, visitor_data->scalename) == 0)
		return 1;
	_destroy_H5DSetDescriptor(h5dimscale);
	return 0;
}

static int get_scale_along(const H5DSetDescriptor *h5dset,
		int along, const char *scalename,
		H5DSetDescriptor *h5dimscale, CharAE *NAME_buf)
{
	int *idx;
	VisitorData visitor_data;

	visitor_data.scalename = scalename;
	visitor_data.h5dimscale = h5dimscale;
	visitor_data.NAME_buf = NAME_buf;
	idx = NULL;
	return H5DSiterate_scales(h5dset->dset_id,
				  (unsigned int) along, idx,
				  visitor, &visitor_data);
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimscales(SEXP filepath, SEXP name, SEXP scalename)
{
	const char *scalename0;
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset, h5dimscale;
	SEXP ans, ans_elt;
	CharAE *NAME_buf;
	int along, ret;

	if (STRING_ELT(scalename, 0) == NA_STRING) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	if (_init_H5DSetDescriptor(&h5dset, dset_id, 0, 0) < 0) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	ans = PROTECT(NEW_CHARACTER(h5dset.ndim));

	NAME_buf = new_CharAE(0);
	for (along = 0; along < h5dset.ndim; along++) {
		ret = get_scale_along(&h5dset, along, scalename0,
				      &h5dimscale, NAME_buf);
		if (ret < 0) {
			_destroy_H5DSetDescriptor(&h5dset);
			H5Dclose(dset_id);
			H5Fclose(file_id);
			error(_HDF5Array_errmsg_buf());
		}
		if (ret == 0) {
			SET_STRING_ELT(ans, along, NA_STRING);
		} else {
			ans_elt = PROTECT(mkChar(h5dimscale.h5name));
			_destroy_H5DSetDescriptor(&h5dimscale);
			SET_STRING_ELT(ans, along, ans_elt);
			UNPROTECT(1);
		}
	}

	_destroy_H5DSetDescriptor(&h5dset);
	H5Dclose(dset_id);
	H5Fclose(file_id);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_h5setdimscales()
 */

static int check_NAME_attribute(hid_t dimscale_id, const char *dimscale,
				const char *scalename, CharAE *NAME_buf)
{
	int ret;

	ret = _get_h5attrib_strval(dimscale_id, "NAME", NAME_buf);
	if (ret < 0)
		return -1;
	if (ret == 0) {
		if (scalename == NULL)
			return 0;
		ret = H5DSis_scale(dimscale_id);
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
					    dimscale, scalename);
			return -1;
		}
		return 0;
	}
	if (scalename == NULL) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to an "
				    "**unnamed** Dimension Scale\n  "
				    "because dataset has a \"NAME\" "
				    "attribute", dimscale);
		return -1;
	}
	if (ret == 1) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "Dimension Scale named \"%s\"\n  "
				    "because dataset already has a \"NAME\" "
				    "attribute that is **not** of class "
				    "H5T_STRING\n  (converting the dataset "
				    "to a Dimension Scale would replace this "
				    "attribute)", dimscale, scalename);
		return -1;
	}
	if (strcmp(NAME_buf->elts, scalename) != 0) {
		PRINT_TO_ERRMSG_BUF("won't convert dataset '%s' to a "
				    "Dimension Scale named \"%s\"\n  "
				    "because dataset already has a \"NAME\" "
				    "attribute set to \"%s\"",
				    dimscale, scalename, NAME_buf->elts);
		return -1;
	}
	return 0;
}

/* Check whether or not dataset 'dimscale' is already set as a Dimension
   Scale named 'scalename' (via its "NAME" attribute) on dimension 'along'
   of dataset 'h5dset'. Return: 0 if no, 1 if yes, -1 if error.
 */
static int check_scale_along(hid_t file_id, const H5DSetDescriptor *h5dset,
		int along, const char *dimscale, const char *scalename,
		CharAE *NAME_buf)
{
	hid_t dimscale_id;
	H5DSetDescriptor h5dimscale, h5dimscale2;
	int ret;

	dimscale_id = H5Dopen(file_id, dimscale, H5P_DEFAULT);
	if (dimscale_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dimscale);
		return -1;
	}
	ret = _init_H5DSetDescriptor(&h5dimscale, dimscale_id, 0, 0);
	if (ret < 0)
		return -1;
	if (strcmp(h5dimscale.h5name, h5dset->h5name) == 0) {
		_destroy_H5DSetDescriptor(&h5dimscale);
		H5Dclose(dimscale_id);
		PRINT_TO_ERRMSG_BUF("cannot set dataset '%s' as a "
				    "Dimension Scale on itself",
				    h5dset->h5name);
		return -1;
	}
	/* Check if 'dimscale_id' has scales on it, in which case we won't be
	   able to turn it into a Dimension Scale. */
	ret = H5Aexists(dimscale_id, "DIMENSION_LIST");
	if (ret < 0) {
		_destroy_H5DSetDescriptor(&h5dimscale);
		H5Dclose(dimscale_id);
		PRINT_TO_ERRMSG_BUF("H5Aexists() returned an error");
		return -1;
	}
	if (ret) {
		_destroy_H5DSetDescriptor(&h5dimscale);
		H5Dclose(dimscale_id);
		PRINT_TO_ERRMSG_BUF("dataset '%s' already has Dimension "
				    "Scales on it so cannot be turned\n  "
				    "itself into a Dimension Scale (a "
				    "Dimension Scale dataset cannot have\n  "
				    "Dimension Scales)",
				    h5dimscale.h5name);
		return -1;
	}
	ret = get_scale_along(h5dset, along, scalename,
			      &h5dimscale2, NAME_buf);
	if (ret < 0) {
		_destroy_H5DSetDescriptor(&h5dimscale);
		H5Dclose(dimscale_id);
		return -1;
	}
	if (ret == 0) {
		/* Now that we know that 'dimscale_id' is not already set as
		   a Dimension Scale named 'scalename' on dimension 'along'
		   of 'h5dset', we also want to make sure that it
		   **would** be ok to create that link. We return -1 if not. */
		ret = check_NAME_attribute(dimscale_id, dimscale,
					   scalename, NAME_buf);
		_destroy_H5DSetDescriptor(&h5dimscale);
		H5Dclose(dimscale_id);
		return ret < 0 ? -1 : 0;
	}
	if (strcmp(h5dimscale2.h5name, h5dimscale.h5name) == 0) {
		_destroy_H5DSetDescriptor(&h5dimscale2);
		_destroy_H5DSetDescriptor(&h5dimscale);
		H5Dclose(dimscale_id);
		return 1;
	}
	if (scalename == NULL) {
		PRINT_TO_ERRMSG_BUF("dimension %d of dataset '%s' is already "
				    "linked to unnamed\n  Dimension Scale "
				    "dataset '%s'", along + 1,
				    h5dset->h5name,
				    h5dimscale2.h5name);
	} else {
		PRINT_TO_ERRMSG_BUF("dimension %d of dataset '%s' is already "
				    "linked to Dimension Scale\n  "
				    "dataset '%s' named \"%s\" (via \"NAME\" "
				    "attribute)", along + 1,
				    h5dset->h5name,
				    h5dimscale2.h5name,
				    scalename);
	}
	_destroy_H5DSetDescriptor(&h5dimscale2);
	_destroy_H5DSetDescriptor(&h5dimscale);
	H5Dclose(dimscale_id);
	return -1;
}

static SEXP check_scales(SEXP filepath, SEXP name, SEXP dimscales,
			 const char *scalename)
{
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;
	int ret, along;
	SEXP ans, dimscale;
	CharAE *NAME_buf;

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	if (_init_H5DSetDescriptor(&h5dset, dset_id, 0, 0) < 0) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	if (LENGTH(dimscales) > h5dset.ndim) {
		PRINT_TO_ERRMSG_BUF("'dimscales' cannot be longer than the "
				    "nb of dimensions of dataset '%s' (%d)",
				    h5dset.h5name, h5dset.ndim);
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
				    h5dset.h5name);
		ret = -1;
		goto on_error;
	}

	ans = PROTECT(NEW_LOGICAL(h5dset.ndim));
	NAME_buf = new_CharAE(0);
	for (along = 0; along < LENGTH(dimscales); along++) {
		dimscale = STRING_ELT(dimscales, along);
		if (dimscale == NA_STRING) {
			LOGICAL(ans)[along] = 0;
			continue;
		}
		ret = check_scale_along(file_id, &h5dset,
					along, CHAR(dimscale), scalename,
					NAME_buf);
		if (ret < 0)
			goto on_error;
		LOGICAL(ans)[along] = !ret;
	}

    on_error:
	_destroy_H5DSetDescriptor(&h5dset);
	H5Dclose(dset_id);
	H5Fclose(file_id);
	UNPROTECT(1);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf());
	return ans;
}

static int set_scale_along(hid_t file_id, hid_t dset_id,
		int along, const char *dimscale, const char *scalename)
{
	hid_t dimscale_id;
	int is_scale, is_attached, ret;

	dimscale_id = H5Dopen(file_id, dimscale, H5P_DEFAULT);
	if (dimscale_id < 0) {
		PRINT_TO_ERRMSG_BUF("failed to open dataset '%s'", dimscale);
		return -1;
	}
	is_scale = H5DSis_scale(dimscale_id);
	if (is_scale < 0) {
		H5Dclose(dimscale_id);
		PRINT_TO_ERRMSG_BUF("H5DSis_scale() returned an error");
		return -1;
	}
	if (is_scale) {
		is_attached = H5DSis_attached(dset_id, dimscale_id,
					      (unsigned int) along);
		if (is_attached < 0) {
			H5Dclose(dimscale_id);
			PRINT_TO_ERRMSG_BUF("H5DSis_attached() "
					    "returned an error");
			return -1;
		}
	} else {
		ret = H5DSset_scale(dimscale_id, scalename);
		if (ret < 0) {
			H5Dclose(dimscale_id);
			PRINT_TO_ERRMSG_BUF("H5DSset_scale() "
					    "returned an error");
			return -1;
		}
		is_attached = 0;
	}
	if (!is_attached) {
		ret = H5DSattach_scale(dset_id, dimscale_id,
				       (unsigned int) along);
		if (ret < 0) {
			H5Dclose(dimscale_id);
			PRINT_TO_ERRMSG_BUF("H5DSattach_scale() "
					    "returned an error");
			return -1;
		}
	}
	H5Dclose(dimscale_id);
	return 0;
}

static void set_scales(SEXP filepath, SEXP name, SEXP dimscales,
		       const char *scalename)
{
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;
	int ret, along;
	SEXP dimscale;

	file_id = _get_file_id(filepath, 0);  /* read/write */
	dset_id = _get_dset_id(file_id, name, filepath);
	if (_init_H5DSetDescriptor(&h5dset, dset_id, 0, 0) < 0) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	for (along = 0; along < LENGTH(dimscales); along++) {
		dimscale = STRING_ELT(dimscales, along);
		if (dimscale == NA_STRING)
			continue;
		ret = set_scale_along(file_id, dset_id,
				      along, CHAR(dimscale), scalename);
		if (ret < 0)  /* should never happen */
			goto on_error;
	}

    on_error:
	_destroy_H5DSetDescriptor(&h5dset);
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ret < 0)
		error(_HDF5Array_errmsg_buf());
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimscales(SEXP filepath, SEXP name, SEXP dimscales, SEXP scalename,
		      SEXP dry_run)
{
	const char *scalename0;
	SEXP ans;
	int along, do_it;

	if (!IS_CHARACTER(dimscales))
		error("'dimscales' must be a character vector");

	if (STRING_ELT(scalename, 0) == NA_STRING) {
		scalename0 = NULL;
	} else {
		scalename0 = CHAR(STRING_ELT(scalename, 0));
	}

	/* We use 2 passes to make the changes to the file atomic. */

	/* 1st pass: dry-run */
	ans = PROTECT(check_scales(filepath, name, dimscales, scalename0));

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
		for (along = 0; along < LENGTH(dimscales); along++) {
			if (LOGICAL(ans)[along]) {
				do_it = 1;
				break;
			}
		}
		if (do_it)
			set_scales(filepath, name, dimscales, scalename0);
	}

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_h5getdimlabels() and C_h5setdimlabels()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_h5getdimlabels(SEXP filepath, SEXP name)
{
	hid_t file_id, dset_id;
	H5DSetDescriptor h5dset;
	int along;
	ssize_t max_label_size, label_size;
	char *label_buf;
	SEXP ans, ans_elt;

	file_id = _get_file_id(filepath, 1);  /* read-only */
	dset_id = _get_dset_id(file_id, name, filepath);
	if (_init_H5DSetDescriptor(&h5dset, dset_id, 0, 0) < 0) {
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error(_HDF5Array_errmsg_buf());
	}

	/* First pass */
	max_label_size = 0;
	for (along = 0; along < h5dset.ndim; along++) {
		label_size = H5DSget_label(dset_id, (unsigned int) along,
					   NULL, 0);
		if (label_size < 0) {
			_destroy_H5DSetDescriptor(&h5dset);
			H5Dclose(dset_id);
			H5Fclose(file_id);
			error("H5DSget_label() returned an error");
		}
		//printf("label_size = %ld\n", label_size);
		if (label_size > max_label_size)
			max_label_size = label_size;
	}

	if (max_label_size == 0) {
		_destroy_H5DSetDescriptor(&h5dset);
		H5Dclose(dset_id);
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
		_destroy_H5DSetDescriptor(&h5dset);
		H5Dclose(dset_id);
		H5Fclose(file_id);
		error("failed to allocate memory for 'label_buf'");
	}
	ans = PROTECT(NEW_CHARACTER(h5dset.ndim));
	for (along = 0; along < h5dset.ndim; along++) {
		label_size = H5DSget_label(dset_id, (unsigned int) along,
					   label_buf, max_label_size + 1);
		/* Should never happen. */
		if (label_size < 0) {
			free(label_buf);
			_destroy_H5DSetDescriptor(&h5dset);
			H5Dclose(dset_id);
			H5Fclose(file_id);
			error("H5DSget_label() returned an error");
		}
		if (label_size > INT_MAX)
			label_size = INT_MAX;
		ans_elt = PROTECT(mkCharLen(label_buf, (int) label_size));
		SET_STRING_ELT(ans, along, ans_elt);
		UNPROTECT(1);
	}

	free(label_buf);
	_destroy_H5DSetDescriptor(&h5dset);
	H5Dclose(dset_id);
	H5Fclose(file_id);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5setdimlabels(SEXP filepath, SEXP name, SEXP dimlabels)
{
	hid_t file_id, dset_id;
	int ndim, along, ret;
	SEXP label;

	if (dimlabels == R_NilValue)
		return R_NilValue;

	file_id = _get_file_id(filepath, 0);
	dset_id = _get_dset_id(file_id, name, filepath);
	ndim = LENGTH(dimlabels);

	for (along = 0; along < ndim; along++) {
		label = STRING_ELT(dimlabels, along);
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

