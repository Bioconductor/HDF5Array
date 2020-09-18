/****************************************************************************
 *                  HDF5Array global error message buffer                   *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "global_errmsg_buf.h"


char * _HDF5Array_global_errmsg_buf()
{
	static char buf[ERRMSG_BUF_LENGTH];

	return buf;
}

