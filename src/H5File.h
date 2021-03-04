#ifndef _H5FILE_H_
#define _H5FILE_H_

/* Avoid spurious MinGW warnings about %lld format specifier. For example:
     H5File.c: In function 'hid_to_string':
     H5File.c:30:18: warning: unknown conversion type character 'l'
                              in format [-Wformat=]
     sprintf(buf, "%lld", (long long) hid);
                     ^
   See https://sourceforge.net/p/mingw-w64/wiki2/gnu%20printf/ for more info. */
#define __USE_MINGW_ANSI_STDIO 1

#include <Rdefines.h>
#include "hdf5.h"

hid_t _h5openlocalfile(SEXP filepath, int readonly);

SEXP C_h5openlocalfile(SEXP filepath, SEXP readonly);

SEXP C_h5openS3file(
	SEXP filepath,
	SEXP auth,
	SEXP aws_region,
	SEXP secret_id,
	SEXP secret_key
);

SEXP C_h5closefile(SEXP ID);

SEXP C_set_H5FileID_xp_ID(SEXP xp, SEXP ID);

SEXP C_get_H5FileID_xp_ID(SEXP xp);

SEXP C_new_H5FileID_xp(SEXP ID);

hid_t _get_file_id(SEXP filepath, int readonly);

const char *_get_file_string(SEXP filepath);

#endif  /* _H5FILE_H_ */

