/****************************************************************************
 *  Some low-level helper functions to support the various h5mread methods  *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "h5mread_helpers.h"

#include "global_errmsg_buf.h"
#include "array_selection.h"
#include "H5DSetDescriptor.h"

#include <stdlib.h>  /* for malloc, free */
#include <zlib.h>  /* for uncompress(), Z_OK, Z_MEM_ERROR, etc.. */


/****************************************************************************
 * Memory management of H5Viewport structs
 */

/* 'mode' controls what fields we should allocate memory for:
     - mode == 0: 'h5off' and 'h5dim' fields only
     - mode == 1: 'off' and 'dim' fields only
     - mode == 2: for all the fields */
int _alloc_H5Viewport(H5Viewport *vp, int ndim, int mode)
{
	vp->h5off = NULL;
	vp->off = NULL;
	if (mode != 1) {
		/* Allocate memory for the 'h5off' and 'h5dim' fields. */
		vp->h5off = _alloc_hsize_t_buf(2 * ndim, 0,
					       "H5Viewport fields");
		if (vp->h5off == NULL)
			return -1;
		vp->h5dim = vp->h5off + ndim;
	}
	if (mode != 0) {
		/* Allocate memory for the 'off' and 'dim' fields. */
		vp->off = (int *) malloc(2 * ndim * sizeof(int));
		if (vp->off == NULL) {
			if (mode != 1)
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
int _alloc_h5chunkvp_middlevp_destvp_bufs(int ndim,
		H5Viewport *h5chunkvp_buf,
		H5Viewport *middlevp_buf,
		H5Viewport *destvp_buf, int destvp_mode)
{
	if (_alloc_H5Viewport(h5chunkvp_buf, ndim, 0) < 0)
		return -1;
	middlevp_buf->h5off =
		_alloc_hsize_t_buf(ndim, 1, "'middlevp_buf->h5off'");
	middlevp_buf->h5dim = h5chunkvp_buf->h5dim;
	if (middlevp_buf->h5off == NULL) {
		_free_H5Viewport(h5chunkvp_buf);
		return -1;
	}
	if (_alloc_H5Viewport(destvp_buf, ndim, destvp_mode) < 0) {
		free(middlevp_buf->h5off);
		_free_H5Viewport(h5chunkvp_buf);
		return -1;
	}
	return 0;
}

/* Used in read_data_4_5(), read_data_7(), and read_data_8(). */
void _free_h5chunkvp_middlevp_destvp_bufs(
		H5Viewport *h5chunkvp_buf,
		H5Viewport *middlevp_buf,
		H5Viewport *destvp_buf)
{
	_free_H5Viewport(destvp_buf);
	free(middlevp_buf->h5off);
	_free_H5Viewport(h5chunkvp_buf);
	return;
}

/* Used in read_data_6(). */
int _alloc_h5chunkvp_innervp_destvp_bufs(int ndim,
		H5Viewport *h5chunkvp_buf,
		H5Viewport *innervp_buf,
		H5Viewport *destvp_buf)
{
	if (_alloc_H5Viewport(h5chunkvp_buf, ndim, 0) < 0)
		return -1;
	if (_alloc_H5Viewport(innervp_buf, ndim, 0) < 0) {
		_free_H5Viewport(h5chunkvp_buf);
		return -1;
	}
	if (_alloc_H5Viewport(destvp_buf, ndim, 2) < 0) {
		_free_H5Viewport(innervp_buf);
		_free_H5Viewport(h5chunkvp_buf);
		return -1;
	}
	return 0;
}

/* Used in read_data_6(). */
void _free_h5chunkvp_innervp_destvp_bufs(
		H5Viewport *h5chunkvp_buf,
		H5Viewport *innervp_buf,
		H5Viewport *destvp_buf)
{
	_free_H5Viewport(destvp_buf);
	_free_H5Viewport(innervp_buf);
	_free_H5Viewport(h5chunkvp_buf);
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

int _add_H5Viewport_to_selection(hid_t space_id, const H5Viewport *vp)
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

int _read_H5Viewport(const H5DSetDescriptor *h5dset,
		const H5Viewport *dsetvp,
		const H5Viewport *memvp,
		void *mem, hid_t mem_space_id)
{
	int ret;

	ret = _select_H5Viewport(h5dset->space_id, dsetvp);
	if (ret < 0)
		return -1;
	ret = _select_H5Viewport(mem_space_id, memvp);
	if (ret < 0)
		return -1;
	ret = H5Dread(h5dset->dset_id,
		      h5dset->mem_type_id, mem_space_id,
		      h5dset->space_id, H5P_DEFAULT, mem);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static void update_h5chunkvp_buf(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const LLongAEAE *tchunkidx_bufs,
		H5Viewport *h5chunkvp_buf)
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
		h5chunkvp_buf->h5off[h5along] = off;
		h5chunkvp_buf->h5dim[h5along] = d;
	}
	//printf("# h5chunkvp_buf->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", h5chunkvp_buf->h5off[h5along]);
	//printf("\n");
	//printf("# h5chunkvp_buf->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", h5chunkvp_buf->h5dim[h5along]);
	//printf("\n");
	return;
}

static void update_destvp_buf(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const IntAEAE *breakpoint_bufs,
		const H5Viewport *h5chunkvp, H5Viewport *destvp_buf)
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
			off = h5chunkvp->h5off[h5along];
			d = h5chunkvp->h5dim[h5along];
		}
		if (destvp_buf->h5off != NULL) {
			destvp_buf->h5off[h5along] = off;
			destvp_buf->h5dim[h5along] = d;
		}
		destvp_buf->off[along] = off;
		destvp_buf->dim[along] = d;
	}
	//printf("# destvp_buf (offsets):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", destvp_buf->off[along]);
	//printf("\n");
	//printf("# destvp_buf (dims):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", destvp_buf->dim[along]);
	//printf("\n");
	return;
}

void _update_h5chunkvp_destvp_bufs(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		H5Viewport *h5chunkvp_buf, H5Viewport *destvp_buf)
{
	update_h5chunkvp_buf(h5dset,
			tchunk_midx, moved_along,
			starts, tchunkidx_bufs,
			h5chunkvp_buf);
	update_destvp_buf(h5dset,
			tchunk_midx, moved_along,
			starts, breakpoint_bufs,
			h5chunkvp_buf, destvp_buf);
	return;
}

int _tchunk_is_fully_selected(int ndim,
		const H5Viewport *h5chunkvp,
		const H5Viewport *destvp)
{
	int ok, along, h5along;

	ok = 1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		//printf("h5chunkvp[%d]=%llu destvp[%d]=%llu\n",
		//       along, h5chunkvp->h5dim[h5along],
		//       along, destvp->h5dim[h5along]);
		if (h5chunkvp->h5dim[h5along] != destvp->h5dim[h5along]) {
			ok = 0;
			break;
		}
	}
	//printf("ok=%d\n", ok);
	return ok;
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

	destLen = uncompressed_size;
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
		const H5Viewport *h5chunkvp,
		void *chunk_data_out,
		void *compressed_chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(h5dset->dset_id,
					h5chunkvp->h5off,
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
			    h5chunkvp->h5off, &filters,
			    compressed_chunk_data_buf);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dread_chunk() returned an error");
		return -1;
	}

	//printf("filters = %u\n", filters);

	//FIXME: This will error if chunk data is not compressed!
	//TODO: Decompress only if chunk data is compressed. There should be
	// a bit in the returned 'filters' that indicates this.
	return uncompress_chunk_data(compressed_chunk_data_buf,
				     chunk_storage_size,
				     chunk_data_out,
				     h5dset->chunk_data_buf_size);
}

