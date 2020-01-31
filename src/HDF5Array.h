#include <Rdefines.h>
#include "S4Vectors_interface.h"

#include "hdf5.h"

#define ERRMSG_BUF_LENGTH 256

#define PRINT_TO_ERRMSG_BUF(...) \
	snprintf(_HDF5Array_errmsg_buf(), ERRMSG_BUF_LENGTH, __VA_ARGS__)

/* A data structure for handling an HDF5 dataset. Collect various information
   about the dataset. What is collected are basically the things needed by the
   C_h5mread() function. */
typedef struct {
	hid_t dset_id, dtype_id, space_id, plist_id, mem_type_id;
	char *storage_mode_attr;
	H5T_class_t H5class;
	size_t size, ans_elt_size, chunk_data_buf_size;
	SEXPTYPE Rtype;
	int ndim, *h5nchunk;
	hsize_t *h5dim, *h5chunkdim;
	H5D_layout_t layout;
} DSetHandle;

/* Like VECTOR_ELT(x, i) except that 'x' can be R_NilValue. */
#define	GET_LIST_ELT(x, i) ((x) != R_NilValue ? VECTOR_ELT(x, i) : R_NilValue)

static inline long long int _get_trusted_elt(SEXP x, int i)
{
	return IS_INTEGER(x) ? (long long int) INTEGER(x)[i] :
			       (long long int) REAL(x)[i];
}


/* array_selection.c */

char * _HDF5Array_errmsg_buf();

int _shallow_check_selection(
	int ndim,
	SEXP starts,
	SEXP counts
);

long long int _check_selection(
	int ndim,
	const long long int *dim,
	SEXP starts,
	SEXP counts,
	int *selection_dim_buf
);

SEXP C_check_selection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

long long int _check_ordered_selection(
	int ndim,
	const long long int *dim,
	SEXP starts,
	SEXP counts,
	int *selection_dim_buf,
	int *nstart_buf,
	int *nblock_buf,
	long long int *last_block_start_buf
);

SEXP C_check_ordered_selection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _selection_can_be_reduced(
	int ndim,
	const int *nstart,
	const int *nblock
);

SEXP _reduce_selection(
	int ndim,
	SEXP starts, SEXP counts,
	const int *selection_dim,
	const int *nblock,
	const long long int *last_block_start
);

SEXP C_reduce_selection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _map_starts_to_chunks(
	int ndim,
	const long long int *dim,
	const long long int *chunkdim,
	SEXP starts,
	int *nstart_buf,
	IntAEAE *breakpoint_bufs,
	LLongAEAE *chunkidx_bufs
);

SEXP C_map_starts_to_chunks(
	SEXP starts,
	SEXP dim,
	SEXP chunkdim
);


/* DSetHandle.c */

hsize_t *_alloc_hsize_t_buf(
	size_t buflength,
	int zeroes,
	const char *what
);

void _close_DSetHandle(DSetHandle *dset_handle);

int _get_DSetHandle(
	hid_t dset_id,
	int as_int,
	int Rtype_only,
	DSetHandle *dset_handle
);

hid_t _get_file_id(SEXP filepath);

hid_t _get_dset_id(
	hid_t file_id,
	SEXP name,
	SEXP filepath
);

SEXP C_destroy_DSetHandle_xp(SEXP xp);

SEXP C_create_DSetHandle_xp(
	SEXP filepath,
	SEXP name,
	SEXP as_integer
);

SEXP C_show_DSetHandle_xp(SEXP xp);

SEXP C_get_h5mread_returned_type(
	SEXP filepath,
	SEXP name,
	SEXP as_integer
);


/* h5mread.c */

SEXP C_h5mread(
	SEXP filepath,
	SEXP name,
	SEXP starts,
	SEXP counts,
	SEXP noreduce,
	SEXP as_integer,
	SEXP method
);

