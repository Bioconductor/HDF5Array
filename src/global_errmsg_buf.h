#ifndef _GLOBAL_ERRMSG_BUF_H_
#define _GLOBAL_ERRMSG_BUF_H_

char * _HDF5Array_global_errmsg_buf();

#define ERRMSG_BUF_LENGTH 256

#define PRINT_TO_ERRMSG_BUF(...) \
	snprintf(_HDF5Array_global_errmsg_buf(), ERRMSG_BUF_LENGTH, __VA_ARGS__)

#endif  /* _GLOBAL_ERRMSG_BUF_H_ */

