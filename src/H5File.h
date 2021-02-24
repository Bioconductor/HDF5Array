#ifndef _H5FILE_H_
#define _H5FILE_H_

#include <Rdefines.h>
#include "hdf5.h"

hid_t _h5openlocalfile(SEXP filepath, int readonly);

SEXP C_h5openlocalfile(SEXP filepath, SEXP readonly);

SEXP C_h5openS3file(SEXP filepath, SEXP s3credentials);

SEXP C_h5closefile(SEXP ID);

SEXP C_set_H5FileID_xp_ID(SEXP xp, SEXP ID);

SEXP C_get_H5FileID_xp_ID(SEXP xp);

SEXP C_new_H5FileID_xp(SEXP ID);

hid_t _get_file_id(SEXP filepath, int readonly);

const char *_get_file_string(SEXP filepath);

#endif  /* _H5FILE_H_ */

