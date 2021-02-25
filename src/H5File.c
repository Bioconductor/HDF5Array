/****************************************************************************
 *                              H5File objects                              *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "H5File.h"

#include <stdlib.h>  /* for strtoll */
#include <limits.h>  /* for LLONG_MAX */


static hid_t string_to_hid(const char *s)
{
	long long hid;
	char *endptr;

	if (s[0] == '\0')
		return -1;
	hid = strtoll(s, &endptr, 10);
	if (endptr[0] != '\0')
		return -1;
	if (hid < 0 || hid == LLONG_MAX)
		return -1;
	return (hid_t) hid;
}

static const char *hid_to_string(hid_t hid)
{
	static char buf[21];

	sprintf(buf, "%lld", (long long) hid);
	return buf;
}

/* 'ID' must be a string representing a plausible h5 id (i.e. non-negative
   hid_t value), or NA_character_. The function returns:
     0: if 'ID' is invalid;
     1: if 'ID' is NA_character_;
     2: if 'ID' is a string representing a plausible h5 id, in which
        case the ID converted to hid_t is stored in '*file_id' (but
        only if 'file_id' is not NULL).
*/
static int valid_ID(SEXP ID, hid_t *file_id)
{
	SEXP ID0;
	hid_t hid;

	if (!(IS_CHARACTER(ID) && LENGTH(ID) == 1))
		return 0;
	ID0 = STRING_ELT(ID, 0);
	if (ID0 == NA_STRING)
		return 1;
	hid = string_to_hid(CHAR(ID0));
	if (hid < 0)
		return 0;
	if (file_id != NULL)
		*file_id = hid;
	return 2;
}


/****************************************************************************
 * C_h5openlocalfile()
 */

hid_t _h5openlocalfile(SEXP filepath, int readonly)
{
	SEXP filepath0;
	herr_t ret;
	unsigned int flags;
	hid_t file_id;

	/* Check 'filepath'. */
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

SEXP C_h5openlocalfile(SEXP filepath, SEXP readonly)
{
	int readonly0;
	hid_t file_id;

	/* Check 'readonly'. */
	if (!(IS_LOGICAL(readonly) && LENGTH(readonly) == 1))
		error("'readonly' must be TRUE or FALSE");
	readonly0 = LOGICAL(readonly)[0];

	file_id = _h5openlocalfile(filepath, readonly0);
	return mkString(hid_to_string(file_id));
}


/****************************************************************************
 * C_h5openS3file()
 *
 * This is just an all-C transcription of what
 *
 *     hdf5::h5read(..., s3=TRUE, s3credentials=s3credentials)
 *
 * does (with R code) to obtain a file ID for the remote HDF5 file. The R
 * code that hdf5::h5read() uses for this boils down to:
 *
 *     fapl <- H5Pcreate("H5P_FILE_ACCESS")
 *     on.exit(H5Pclose(fapl))
 *     H5Pset_fapl_ros3(fapl, s3credentials)
 *     loc <- rhdf5:::h5checktypeOrOpenLocS3(filepath, readonly=TRUE,
 *                                                     fapl=fapl,
 *                                                     native=FALSE)
 *     ID <- loc$H5Identifier@ID
 */

#ifdef H5_HAVE_ROS3_VFD
static int set_fapl_ros3(hid_t fapl_id, int auth,
					const char *aws_region,
					const char *secret_id,
					const char *secret_key)
{
	H5FD_ros3_fapl_t fapl;
	int n;

	/* See Rhdf5lib/src/hdf5/src/H5FDros3.h for the definition of
	   the H5FD_ros3_fapl_t struct. */
	fapl.version = H5FD_CURR_ROS3_FAPL_T_VERSION;
	fapl.authenticate = (hbool_t) auth;
	n = snprintf(fapl.aws_region, H5FD_ROS3_MAX_REGION_LEN + 1,
		     "%s", aws_region);
	if (n < 0 || n > H5FD_ROS3_MAX_REGION_LEN)
		return -1;
	n = snprintf(fapl.secret_id, H5FD_ROS3_MAX_SECRET_ID_LEN + 1,
		     "%s", secret_id);
	if (n < 0 || n > H5FD_ROS3_MAX_SECRET_ID_LEN)
		return -1;
	n = snprintf(fapl.secret_key, H5FD_ROS3_MAX_SECRET_KEY_LEN + 1,
		     "%s", secret_key);
	if (n < 0 || n > H5FD_ROS3_MAX_SECRET_KEY_LEN)
		return -1;
	return H5Pset_fapl_ros3(fapl_id, &fapl);
}
#endif

static hid_t h5openS3file(const char *url, int auth,
					   const char *aws_region,
					   const char *secret_id,
					   const char *secret_key)
{
#ifdef H5_HAVE_ROS3_VFD
	int ret;
	hid_t fapl_id, file_id;

	ret = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	if (ret < 0)
		error("H5Eset_auto() returned an error");

	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	if (fapl_id < 0)
		error("H5Pcreate() returned an error");
	ret = set_fapl_ros3(fapl_id, auth, aws_region, secret_id, secret_key);
	if (ret < 0) {
		H5Pclose(fapl_id);
		error("set_fapl_ros3() returned an error");
	}
	file_id = H5Fopen(url, H5F_ACC_RDONLY, fapl_id);
	H5Pclose(fapl_id);
	if (file_id < 0)
		error("failed to open file '%s'", url);
	return file_id;
#else
	error("Rhdf5lib was not compiled with support for the S3 VFD");
#endif
	return 0;
}

SEXP C_h5openS3file(SEXP filepath, SEXP auth,
		    SEXP aws_region, SEXP secret_id, SEXP secret_key)
{
	SEXP filepath0, aws_region0, secret_id0, secret_key0;
	int auth0;
	hid_t file_id;

	/* Check 'filepath'. */
	if (!(IS_CHARACTER(filepath) && LENGTH(filepath) == 1))
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' cannot be NA");

	/* Check 'auth'. */
	if (!(IS_LOGICAL(auth) && LENGTH(auth) == 1))
		error("'auth' must be TRUE or FALSE");
	auth0 = LOGICAL(auth)[0];

	/* Check 'aws_region'. */
	if (!(IS_CHARACTER(aws_region) && LENGTH(aws_region) == 1))
		error("'aws_region' must be a single string");
	aws_region0 = STRING_ELT(aws_region, 0);
	if (aws_region0 == NA_STRING)
		error("'aws_region' cannot be NA");

	/* Check 'secret_id'. */
	if (!(IS_CHARACTER(secret_id) && LENGTH(secret_id) == 1))
		error("'secret_id' must be a single string");
	secret_id0 = STRING_ELT(secret_id, 0);
	if (secret_id0 == NA_STRING)
		error("'secret_id' cannot be NA");

	/* Check 'secret_key'. */
	if (!(IS_CHARACTER(secret_key) && LENGTH(secret_key) == 1))
		error("'secret_key' must be a single string");
	secret_key0 = STRING_ELT(secret_key, 0);
	if (secret_key0 == NA_STRING)
		error("'secret_key' cannot be NA");

	file_id = h5openS3file(CHAR(filepath0), auth0,
						CHAR(aws_region0),
						CHAR(secret_id0),
						CHAR(secret_key0));
	return mkString(hid_to_string(file_id));
}


/****************************************************************************
 * C_h5closefile()
 */

SEXP C_h5closefile(SEXP ID)
{
	int ret;
	hid_t file_id;

	ret = valid_ID(ID, &file_id);
	if (ret == 0)
		error("invalid H5FileID object (invalid 'ID')");
	if (ret != 1)
		ret = H5Fclose(file_id);
	return R_NilValue;
}


/****************************************************************************
 * H5FileID objects
 */

/* Strangely R_SetExternalPtrTag() and R_ExternalPtrTag() work on any type of
   SEXP and manage to return something. These wrappers check that their first
   argument is an external pointer and return an error is they are not. */
static SEXP set_xp_tag(SEXP xp, SEXP tag)
{
	if (TYPEOF(xp) != EXTPTRSXP)
		error("invalid H5FileID object "
		      "('xp' slot not an external pointer)");
	R_SetExternalPtrTag(xp, tag);
	return xp;
}

static SEXP get_xp_tag(SEXP xp)
{
	if (TYPEOF(xp) != EXTPTRSXP)
		error("invalid H5FileID object "
		      "('xp' slot not an external pointer)");
	return R_ExternalPtrTag(xp);
}

/* --- .Call ENTRY POINT --- */
SEXP C_set_H5FileID_xp_ID(SEXP xp, SEXP ID)
{
	if (!valid_ID(ID, NULL))
		error("supplied 'ID' is invalid");
	return set_xp_tag(xp, ID);
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_H5FileID_xp_ID(SEXP xp)
{
	return get_xp_tag(xp);
}

/* --- .Call ENTRY POINT --- */
SEXP C_new_H5FileID_xp(SEXP ID)
{
	if (!valid_ID(ID, NULL))
		error("supplied 'ID' is invalid");
	return R_MakeExternalPtr(NULL, ID, R_NilValue);
}


/****************************************************************************
 * _get_file_id()
 */

static hid_t ID_to_file_id(SEXP ID)
{
	int ret;
	hid_t file_id;

	ret = valid_ID(ID, &file_id);
	if (ret == 0)
		error("invalid H5File object (invalid 'ID')");
	if (ret == 1)
		error("H5File object is closed");
	return file_id;
}

static hid_t get_H5File_file_id(SEXP x)
{
	SEXP HDF5Array_h5id, xp, ID;

	HDF5Array_h5id = GET_SLOT(x, install("HDF5Array_h5id"));
	xp = GET_SLOT(HDF5Array_h5id, install("xp"));
	ID = get_xp_tag(xp);
	return ID_to_file_id(ID);
}

static SEXP get_H5File_filepath(SEXP x)
{
	return GET_SLOT(x, install("filepath"));
}

/* 'filepath' can be a single string or an H5File object:
     - If a single string, it's considered to be the path to a local HDF5
       file. The file is opened and its h5 id returned.
     - If an H5File object, the h5 id is extracted from the object.
*/
hid_t _get_file_id(SEXP filepath, int readonly)
{
	const char *class;

	if (!isObject(filepath))
		return _h5openlocalfile(filepath, readonly);
	class = CHAR(STRING_ELT(GET_CLASS(filepath), 0));
	if (strcmp(class, "H5File") != 0)
		error("'filepath' must be a single string "
		      "or an H5File object");
	if (!readonly)
		error("H5File objects are read-only and cannot be "
		      "accessed in read/write mode at the moment");
	return get_H5File_file_id(filepath);
}

/* 'filepath' should have been already validated by _get_file_id() above
   so can be trusted. */
const char *_get_file_string(SEXP filepath)
{
	if (isObject(filepath))
		filepath = get_H5File_filepath(filepath);
	return CHAR(STRING_ELT(filepath, 0));
}

