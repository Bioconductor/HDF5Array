/****************************************************************************
 *  Some low-level helper functions to support the various h5mread methods  *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_helpers.h"

/*
Some useful links:
- Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
    https://support.hdfgroup.org/HDF5/doc/RM/RM_H5S.html
- Documentation of H5Dread():
    https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
- An H5Dread() example:
    https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
*/

#include "global_errmsg_buf.h"
#include "uaselection.h"
#include "H5DSetDescriptor.h"

#include <stdlib.h>  /* for malloc, free */
#include <zlib.h>  /* for uncompress(), Z_OK, Z_MEM_ERROR, etc.. */


static void print_chunk_data(const H5DSetDescriptor *h5dset,
			     void *chunk_data_buf)
{
	printf("chunk_data_buf:");
/*
	for (size_t i = 0; i < h5dset->chunk_data_buf_size; i++) {
		if (i % 12 == 0)
			printf("\n ");
		printf(" '%c'", ((char *) chunk_data_buf)[i]);
	}
*/
	size_t nval = h5dset->chunk_data_buf_size / h5dset->ans_elt_size;
	for (size_t i = 0; i < nval; i++) {
		if (i % 12 == 0)
			printf("\n ");
		printf(" %4d", ((int *) chunk_data_buf)[i]);
	}
	printf("\n");
	return;
}


/****************************************************************************
 * Memory management of H5Viewport structs
 */

int _alloc_H5Viewport(H5Viewport *vp, int ndim, int mode)
{
	vp->h5off = NULL;
	vp->off = NULL;
	if (mode != ALLOC_OFF_AND_DIM) {
		/* Allocate memory for the 'h5off' and 'h5dim' fields. */
		vp->h5off = _alloc_hsize_t_buf(2 * ndim, 0,
					       "H5Viewport fields");
		if (vp->h5off == NULL)
			return -1;
		vp->h5dim = vp->h5off + ndim;
	}
	if (mode != ALLOC_H5OFF_AND_H5DIM) {
		/* Allocate memory for the 'off' and 'dim' fields. */
		vp->off = (int *) malloc(2 * ndim * sizeof(int));
		if (vp->off == NULL) {
			if (mode != ALLOC_OFF_AND_DIM)
				free(vp->h5off);
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for H5Viewport fields");
			return -1;
		}
		vp->dim = vp->off + ndim;
	}
	return 0;
}

void _free_H5Viewport(H5Viewport *vp)
{
	if (vp->h5off != NULL)
		free(vp->h5off);
	if (vp->off != NULL)
		free(vp->off);
	return;
}

/* Used in read_data_4_5(), read_data_7(), and read_data_8(). */
int _alloc_tchunk_vp_middle_vp_dest_vp(int ndim,
		H5Viewport *tchunk_vp,
		H5Viewport *middle_vp,
		H5Viewport *dest_vp, int dest_vp_mode)
{
	if (_alloc_H5Viewport(tchunk_vp, ndim, ALLOC_H5OFF_AND_H5DIM) < 0)
		return -1;
	middle_vp->h5off = _alloc_hsize_t_buf(ndim, 1, "'middle_vp->h5off'");
	if (middle_vp->h5off == NULL) {
		_free_H5Viewport(tchunk_vp);
		return -1;
	}
	middle_vp->h5dim = tchunk_vp->h5dim;
	if (_alloc_H5Viewport(dest_vp, ndim, dest_vp_mode) < 0) {
		free(middle_vp->h5off);
		_free_H5Viewport(tchunk_vp);
		return -1;
	}
	return 0;
}

/* Used in read_data_4_5(), read_data_7(), and read_data_8(). */
void _free_tchunk_vp_middle_vp_dest_vp(
		H5Viewport *tchunk_vp,
		H5Viewport *middle_vp,
		H5Viewport *dest_vp)
{
	_free_H5Viewport(dest_vp);
	free(middle_vp->h5off);
	_free_H5Viewport(tchunk_vp);
	return;
}

/* Used in read_data_6(). */
int _alloc_tchunk_vp_inner_vp_dest_vp(int ndim,
		H5Viewport *tchunk_vp,
		H5Viewport *inner_vp,
		H5Viewport *dest_vp)
{
	if (_alloc_H5Viewport(tchunk_vp, ndim, ALLOC_H5OFF_AND_H5DIM) < 0)
		return -1;
	if (_alloc_H5Viewport(inner_vp, ndim, ALLOC_H5OFF_AND_H5DIM) < 0) {
		_free_H5Viewport(tchunk_vp);
		return -1;
	}
	if (_alloc_H5Viewport(dest_vp, ndim, ALLOC_ALL_FIELDS) < 0) {
		_free_H5Viewport(inner_vp);
		_free_H5Viewport(tchunk_vp);
		return -1;
	}
	return 0;
}

/* Used in read_data_6(). */
void _free_tchunk_vp_inner_vp_dest_vp(
		H5Viewport *tchunk_vp,
		H5Viewport *inner_vp,
		H5Viewport *dest_vp)
{
	_free_H5Viewport(dest_vp);
	_free_H5Viewport(inner_vp);
	_free_H5Viewport(tchunk_vp);
	return;
}


/****************************************************************************
 * Other helpers
 */

int _map_starts_to_h5chunks(const H5DSetDescriptor *h5dset,
		SEXP starts,
		int *nstart_buf,
		IntAEAE *breakpoint_bufs, LLongAEAE *tchunkidx_bufs)
{
	int ndim, along, h5along;
	LLongAE *dim_buf, *chunkdim_buf;

	ndim = h5dset->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	chunkdim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		dim_buf->elts[along] =
			(long long int) h5dset->h5dim[h5along];
		chunkdim_buf->elts[along] =
			(long long int) h5dset->h5chunkdim[h5along];
	}
	return _map_starts_to_chunks(ndim, dim_buf->elts, chunkdim_buf->elts,
				     starts,
				     nstart_buf,
				     breakpoint_bufs, tchunkidx_bufs);
}

long long int _set_num_tchunks(const H5DSetDescriptor *h5dset,
		const SEXP starts,
		const LLongAEAE *tchunkidx_bufs,
		int *num_tchunks_buf)
{
	int ndim, along, h5along, n;
	long long int total_num_tchunks;  /* total nb of touched chunks */
	SEXP start;

	ndim = h5dset->ndim;
	total_num_tchunks = 1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			n = LLongAE_get_nelt(tchunkidx_bufs->elts[along]);
		} else {
			n = h5dset->h5nchunk[h5along];
		}
		total_num_tchunks *= num_tchunks_buf[along] = n;
	}
	return total_num_tchunks;
}

hid_t _create_mem_space(int ndim, const int *dim)
{
	hsize_t *h5dim;
	int along, h5along;
	hid_t mem_space_id;

	/* Allocate and set 'h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		return -1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		h5dim[h5along] = dim[along];
	mem_space_id = H5Screate_simple(ndim, h5dim, NULL);
	if (mem_space_id < 0)
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
	free(h5dim);
	return mem_space_id;
}

int _select_H5Viewport(hid_t space_id, const H5Viewport *vp)
{
	int ret;

	ret = H5Sselect_hyperslab(space_id, H5S_SELECT_SET,
				  vp->h5off, NULL, vp->h5dim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

int _add_H5Viewport_to_h5selection(hid_t space_id, const H5Viewport *vp)
{
	int ret;

	ret = H5Sselect_hyperslab(space_id, H5S_SELECT_OR,
				  vp->h5off, NULL, vp->h5dim, NULL);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_hyperslab() returned an error");
		return -1;
	}
	return 0;
}

int _read_H5Viewport(const H5DSetDescriptor *h5dset,
		const H5Viewport *h5dset_vp,
		const H5Viewport *mem_vp,
		void *mem, hid_t mem_space_id)
{
	int ret;

	ret = _select_H5Viewport(h5dset->space_id, h5dset_vp);
	if (ret < 0)
		return -1;
	ret = _select_H5Viewport(mem_space_id, mem_vp);
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	//print_chunk_data(h5dset, mem);
	return ret;
}

int _read_h5selection(const H5DSetDescriptor *h5dset,
		const H5Viewport *mem_vp,
		void *mem, hid_t mem_space_id)
{
	int ret;

	if (mem_vp == NULL) {
		ret = H5Sselect_all(mem_space_id);
		if (ret < 0)
		    PRINT_TO_ERRMSG_BUF("H5Sselect_all() returned an error");
	} else {
		ret = _select_H5Viewport(mem_space_id, mem_vp);
	}
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static void update_tchunk_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const LLongAEAE *tchunkidx_bufs,
		H5Viewport *tchunk_vp)
{
	int ndim, along, h5along, i;
	SEXP start;
	long long int tchunkidx;
	hsize_t chunkd, off, d;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = tchunk_midx[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue) {
			tchunkidx = tchunkidx_bufs->elts[along]->elts[i];
		} else {
			tchunkidx = i;
		}
		chunkd = h5dset->h5chunkdim[h5along];
		off = tchunkidx * chunkd;
		d = h5dset->h5dim[h5along] - off;
		if (d > chunkd)
			d = chunkd;
		tchunk_vp->h5off[h5along] = off;
		tchunk_vp->h5dim[h5along] = d;
	}
	//printf("# tchunk_vp->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", tchunk_vp->h5off[h5along]);
	//printf("\n");
	//printf("# tchunk_vp->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", tchunk_vp->h5dim[h5along]);
	//printf("\n");
	return;
}

static void update_dest_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const IntAEAE *breakpoint_bufs,
		const H5Viewport *tchunk_vp, H5Viewport *dest_vp)
{
	int ndim, along, h5along, i, off, d;
	SEXP start;
	const int *breakpoint;

	ndim = h5dset->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = tchunk_midx[along];
		start = GET_LIST_ELT(starts, along);
		if (start != R_NilValue ) {
			breakpoint = breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = tchunk_vp->h5off[h5along];
			d = tchunk_vp->h5dim[h5along];
		}
		if (dest_vp->h5off != NULL) {
			dest_vp->h5off[h5along] = off;
			dest_vp->h5dim[h5along] = d;
		}
		dest_vp->off[along] = off;
		dest_vp->dim[along] = d;
	}
	//printf("# dest_vp (offsets):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", dest_vp->off[along]);
	//printf("\n");
	//printf("# dest_vp (dims):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", dest_vp->dim[along]);
	//printf("\n");
	return;
}

void _update_tchunk_vp_dest_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		H5Viewport *tchunk_vp, H5Viewport *dest_vp)
{
	update_tchunk_vp(h5dset,
			tchunk_midx, moved_along,
			starts, tchunkidx_bufs,
			tchunk_vp);
	update_dest_vp(h5dset,
			tchunk_midx, moved_along,
			starts, breakpoint_bufs,
			tchunk_vp, dest_vp);
	return;
}

/* Return 1 if the chunk that 'tchunk_vp' is pointing at is "truncated"
   (a.k.a. "partial edge chunk" in HDF5's terminology), and 0 otherwise
   (i.e. if the new chunk is a full-size chunk). */
int _tchunk_is_truncated(const H5DSetDescriptor *h5dset,
			 const H5Viewport *tchunk_vp)
{
	int ndim, h5along;
	hsize_t chunkd, d;

	ndim = h5dset->ndim;
	for (h5along = 0; h5along < ndim; h5along++) {
		chunkd = h5dset->h5chunkdim[h5along];
		d = tchunk_vp->h5dim[h5along];
		if (d != chunkd)
			return 1;
	}
	return 0;
}

int _tchunk_is_fully_selected(int ndim,
		const H5Viewport *tchunk_vp,
		const H5Viewport *dest_vp)
{
	int along, h5along, not_fully;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		not_fully = tchunk_vp->h5dim[h5along] !=
			    (hsize_t) dest_vp->dim[along];
		if (not_fully)
			return 0;
	}
	return 1;
}


/****************************************************************************
 * _read_h5chunk()
 */

static int uncompress_chunk_data(const void *compressed_chunk_data,
				 hsize_t compressed_size,
				 void *uncompressed_chunk_data,
				 size_t uncompressed_size)
{
	int ret;
	uLong destLen;

	destLen = (uLong) uncompressed_size;
	ret = uncompress((Bytef *) uncompressed_chunk_data, &destLen,
			 compressed_chunk_data, (uLong) compressed_size);
	if (ret == Z_OK) {
		if (destLen == uncompressed_size)
			return 0;
		PRINT_TO_ERRMSG_BUF("error in uncompress_chunk_data(): "
				    "chunk data smaller than expected "
				    "after decompression");
		return -1;
	}
	switch (ret) {
	    case Z_MEM_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "not enough memory to uncompress chunk");
	    break;
	    case Z_BUF_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "not enough room in output buffer");
	    break;
	    case Z_DATA_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "chunk data corrupted or incomplete");
	    break;
	    default:
		PRINT_TO_ERRMSG_BUF("unknown error in uncompress()");
	}
	return -1;
}

static void transpose_bytes(const char *in, size_t nrow, size_t ncol, char *out)
{
	size_t i, j, in_offset;

	for (i = 0; i < nrow; i++) {
		in_offset = i;
		for (j = 0; j < ncol; j++) {
			*(out++) = *(in + in_offset);
			in_offset += nrow;
		}
	}
	return;
}

/*
 * Unfortunately H5Dread_chunk() is NOT listed here:
 *   https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html
 * Header file for declaration:
 *   hdf5-1.10.3/src/H5Dpublic.h
 * See hdf5-1.10.3/test/direct_chunk.c for plenty of examples.
 *
 * Call stack for H5Dread_chunk()
 *   H5Dread_chunk                (H5Dio.c)
 *     H5D__chunk_direct_read     (H5Dchunk.c)
 *       H5F_block_read           (H5Fio.c)
 *         H5PB_read              (H5PB.c)
 *           H5F__accum_read      (H5Faccum.c)
 *             or
 *           H5FD_read            (H5FDint.c)
 *            ??
 *
 * Call stack for H5Dread()
 *   H5Dread                      (H5Dio.c)
 *     H5D__read                  (H5Dio.c)
 *       H5D__chunk_read          (H5Dchunk.c)
 *         H5D__select_read
 *           or
 *         H5D__scatgath_read     (H5Dscatgath.c)
 *           H5D__gather_file     (H5Dscatgath.c)
 *       call ser_read member of a H5D_layout_ops_t object
 *            ??
 */
#define CHUNK_COMPRESSION_OVERHEAD 8  // empirical (increase if necessary)

int _read_h5chunk(const H5DSetDescriptor *h5dset,
		const H5Viewport *h5chunk_vp,
		void *compressed_chunk_data_buf,
		void *chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(h5dset->dset_id,
					h5chunk_vp->h5off,
					&chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > h5dset->chunk_data_buf_size +
				 CHUNK_COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%lu + %d)",
				    chunk_storage_size,
				    h5dset->chunk_data_buf_size,
				    CHUNK_COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(h5dset->dset_id, H5P_DEFAULT,
			    h5chunk_vp->h5off, &filters,
			    compressed_chunk_data_buf);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dread_chunk() returned an error");
		return -1;
	}

	//printf("filters = %u\n", filters);

	//FIXME: This will error if chunk data is not compressed!
	//TODO: Decompress only if chunk data is compressed. There should be
	//a bit in the returned 'filters' that indicates this.
	ret = uncompress_chunk_data(compressed_chunk_data_buf,
				    chunk_storage_size,
				    chunk_data_buf,
				    h5dset->chunk_data_buf_size);
	if (ret < 0)
		return -1;
	size_t nval = h5dset->chunk_data_buf_size / h5dset->ans_elt_size;
	transpose_bytes(chunk_data_buf, nval, h5dset->ans_elt_size,
			compressed_chunk_data_buf);
	//print_chunk_data(h5dset, compressed_chunk_data_buf);
	return 0;
}

