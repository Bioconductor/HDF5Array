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

static hid_t ID_to_file_id(SEXP ID)
{
	int ret;
	hid_t file_id;

	ret = valid_ID(ID, &file_id);
	if (ret == 0)
		error("invalid H5FileID object (invalid 'ID')");
	if (ret == 1)
		error("H5FileID object is closed");
	return file_id;
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

SEXP C_h5openS3file(SEXP filepath, SEXP s3credentials)
{
	error("'s3=TRUE' is not supported yet");
	return R_NilValue;
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

static hid_t get_H5File_file_id(SEXP x)
{
	SEXP h5fid, xp, ID;

	h5fid = GET_SLOT(x, install("h5fid"));
	xp = GET_SLOT(h5fid, install("xp"));
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

