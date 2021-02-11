/****************************************************************************
 *               Iterating over the chunks of an HDF5 dataset               *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "ChunkIterator.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"

#include <stdlib.h>  /* for malloc, free */
#include <zlib.h>  /* for uncompress(), Z_OK, Z_MEM_ERROR, etc.. */


/****************************************************************************
 * Low-level helpers (non-exported)
 */

static int alloc_h5dset_vp_mem_vp(int ndim,
		H5Viewport *h5dset_vp,
		H5Viewport *mem_vp, int mem_vp_mode)
{
	if (_alloc_H5Viewport(h5dset_vp, ndim, ALLOC_H5OFF_AND_H5DIM) < 0)
		return -1;
	if (_alloc_H5Viewport(mem_vp, ndim, mem_vp_mode) < 0) {
		_free_H5Viewport(h5dset_vp);
		return -1;
	}
	return 0;
}

static void free_h5dset_vp_mem_vp(H5Viewport *h5dset_vp, H5Viewport *mem_vp)
{
	_free_H5Viewport(mem_vp);
	_free_H5Viewport(h5dset_vp);
	return;
}

static int map_starts_to_h5chunks(const H5DSetDescriptor *h5dset,
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

static long long int set_num_tchunks(const H5DSetDescriptor *h5dset,
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

static void update_h5dset_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const LLongAEAE *tchunkidx_bufs,
		H5Viewport *h5dset_vp)
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
		h5dset_vp->h5off[h5along] = off;
		h5dset_vp->h5dim[h5along] = d;
	}
	//printf("# h5dset_vp->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", h5dset_vp->h5off[h5along]);
	//printf("\n");
	//printf("# h5dset_vp->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      printf(" %llu", h5dset_vp->h5dim[h5along]);
	//printf("\n");
	return;
}

static void update_mem_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts, const IntAEAE *breakpoint_bufs,
		const H5Viewport *h5dset_vp, H5Viewport *mem_vp)
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
			off = h5dset_vp->h5off[h5along];
			d = h5dset_vp->h5dim[h5along];
		}
		if (mem_vp->h5off != NULL) {
			mem_vp->h5off[h5along] = off;
			mem_vp->h5dim[h5along] = d;
		}
		mem_vp->off[along] = off;
		mem_vp->dim[along] = d;
	}
	//printf("# mem_vp (offsets):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", mem_vp->off[along]);
	//printf("\n");
	//printf("# mem_vp (dims):");
	//for (along = 0; along < ndim; along++)
	//      printf(" %d", mem_vp->dim[along]);
	//printf("\n");
	return;
}

static void update_h5dset_vp_mem_vp(const H5DSetDescriptor *h5dset,
		const int *tchunk_midx, int moved_along,
		SEXP starts,
		const IntAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		H5Viewport *h5dset_vp, H5Viewport *mem_vp)
{
	update_h5dset_vp(h5dset,
			tchunk_midx, moved_along,
			starts, tchunkidx_bufs,
			h5dset_vp);
	update_mem_vp(h5dset,
			tchunk_midx, moved_along,
			starts, breakpoint_bufs,
			h5dset_vp, mem_vp);
	return;
}


/****************************************************************************
 * read_h5chunk()
 *
 * Based on H5Dread_chunk(), which is NOT listed here for some mysterious
 * reasons: https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html
 *
 * Header file for declaration: hdf5-1.10.3/src/H5Dpublic.h
 *
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
 *             ??
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
static void print_chunk_data(void *data, size_t data_length, size_t data_size)
{
	printf("chunk data:");
	//for (size_t i = 0; i < data_size; i++) {
	//	if (i % 12 == 0)
	//		printf("\n ");
	//	printf(" '%c'", ((char *) data)[i]);
	//}
	for (size_t i = 0; i < data_length; i++) {
		if (i % 12 == 0)
			printf("\n ");
		printf(" %4d", ((int *) data)[i]);
	}
	printf("\n");
	return;
}
*/

#define	CHUNK_COMPRESSION_OVERHEAD 8  // empirical (increase if necessary)

/* WARNING: read_h5chunk() is not ready yet! It is NOT working properly on
   some datasets:
      library(HDF5Array)
      library(ExperimentHub)
      hub <- ExperimentHub()
      fname0 <- hub[["EH1039"]]
      h5mread(fname0, "mm10/barcodes", list(1), method=4L)
      # [1] "AAACCTGAGATAGGAG-1"
      h5mread(fname0, "mm10/barcodes", list(1),
              method=4L, use.H5Dread_chunk=TRUE)
      # [1] "AAAAAAAAAAAAAAAAAAAA"
   Looks like the chunk data has been shuffled (transposed in that case)
   before being written to disk in order to improve compression.
   TODO: Investigate this further. I suspect we need to check whether a
   "Data shuffling filter" (H5Z_FILTER_SHUFFLE) was used at creation time.
   Check H5Pget_filter() for how to know whether this filter was used or not.
   There should be a way to retrieve information about how the data was
   shuffled. */
static int read_h5chunk(hid_t dset_id,
		const H5Viewport *h5chunk_vp,
		ChunkDataBuffer *chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(dset_id,
					h5chunk_vp->h5off,
					&chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > chunk_data_buf->data_size +
				 CHUNK_COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%lu + %d)",
				    chunk_storage_size,
				    chunk_data_buf->data_size,
				    CHUNK_COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(dset_id, H5P_DEFAULT,
			    h5chunk_vp->h5off, &filters,
			    chunk_data_buf->compressed_data);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dread_chunk() returned an error");
		return -1;
	}

	//printf("filters = %u\n", filters);

	//FIXME: This will error if chunk data is not compressed!
	//TODO: Decompress only if chunk data is compressed. There should be
	//a bit in the returned 'filters' that indicates this.
	ret = uncompress_chunk_data(chunk_data_buf->compressed_data,
				    chunk_storage_size,
				    chunk_data_buf->data,
				    chunk_data_buf->data_size);
	if (ret < 0)
		return -1;
	transpose_bytes(chunk_data_buf->data,
			chunk_data_buf->data_length,
			chunk_data_buf->data_type_size,
			chunk_data_buf->compressed_data);
	//print_chunk_data(chunk_data_buf->compressed_data,
	//		   chunk_data_buf->data_length,
	//		   chunk_data_buf->data_size);
	return 0;
}


/****************************************************************************
 * Exported functions
 */

void _destroy_ChunkIterator(ChunkIterator *chunk_iter)
{
	if (chunk_iter->h5dset_vp.h5off != NULL)
		free_h5dset_vp_mem_vp(&chunk_iter->h5dset_vp,
				      &chunk_iter->mem_vp);
	return;
}

int _init_ChunkIterator(ChunkIterator *chunk_iter,
		const H5DSetDescriptor *h5dset, SEXP index,
		int *selection_dim,
		int alloc_full_mem_vp)
{
	int ndim, ret;

	if (h5dset->h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("'h5dset->h5chunkdim' is NULL");
		return -1;
	}

	chunk_iter->h5dset = h5dset;
	chunk_iter->index = index;
	ndim = h5dset->ndim;

	/* Initialize ChunkIterator struct members that control
	   what _destroy_ChunkIterator() needs to free or close. */
	chunk_iter->h5dset_vp.h5off = NULL;

	/* Set struct members 'breakpoint_bufs' and 'tchunkidx_bufs'.
	   Also populate 'selection_dim' if not set to NULL. */
	chunk_iter->breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunk_iter->tchunkidx_bufs = new_LLongAEAE(ndim, ndim);
	ret = map_starts_to_h5chunks(h5dset, index, selection_dim,
				     chunk_iter->breakpoint_bufs,
				     chunk_iter->tchunkidx_bufs);
	if (ret < 0)
		goto on_error;

	/* Set struct members 'num_tchunks' and 'total_num_tchunks'. */
	chunk_iter->num_tchunks = new_IntAE(ndim, ndim, 0)->elts;
	chunk_iter->total_num_tchunks = set_num_tchunks(h5dset, index,
						chunk_iter->tchunkidx_bufs,
						chunk_iter->num_tchunks);

	/* Allocate struct members 'h5dset_vp' and 'mem_vp'. */
	ret = alloc_h5dset_vp_mem_vp(ndim,
				&chunk_iter->h5dset_vp,
				&chunk_iter->mem_vp,
				alloc_full_mem_vp ? ALLOC_ALL_FIELDS
						  : ALLOC_OFF_AND_DIM);
	if (ret < 0)
		goto on_error;

	/* Set struct member 'tchunk_midx_buf'. */
	chunk_iter->tchunk_midx_buf = new_IntAE(ndim, ndim, 0)->elts;

	/* Set struct member 'tchunk_rank'. */
	chunk_iter->tchunk_rank = -1;
	return 0;

    on_error:
	_destroy_ChunkIterator(chunk_iter);
	return -1;
}

/* Return:
 *     1 = if the chunk before the move was not the last chunk and the move to
 *         the next chunk was successful;
 *     0 = if the chunk before the move was the last chunk and so the move to
 *         the next chunk was not possible;
 *   < 0 = if error
 * Typical use:
 *     while (ret = _next_chunk(&chunk_iter)) {
 *         if (ret < 0) {
 *             an error occured
 *         }
 *         handle current chunk
 *     }
 */
int _next_chunk(ChunkIterator *chunk_iter)
{
	const H5DSetDescriptor *h5dset;

	chunk_iter->tchunk_rank++;
	if (chunk_iter->tchunk_rank == chunk_iter->total_num_tchunks)
		return 0;
	h5dset = chunk_iter->h5dset;
	chunk_iter->moved_along = chunk_iter->tchunk_rank == 0 ?
					h5dset->ndim :
					_next_midx(h5dset->ndim,
						chunk_iter->num_tchunks,
						chunk_iter->tchunk_midx_buf);
	update_h5dset_vp_mem_vp(h5dset,
			chunk_iter->tchunk_midx_buf,
			chunk_iter->moved_along,
			chunk_iter->index,
			chunk_iter->breakpoint_bufs, chunk_iter->tchunkidx_bufs,
			&chunk_iter->h5dset_vp, &chunk_iter->mem_vp);
	return 1;
}

void _print_tchunk_info(const ChunkIterator *chunk_iter)
{
	int ndim, along, h5along, i;
	long long int tchunkidx;

	Rprintf("processing chunk %lld/%lld: [",
		chunk_iter->tchunk_rank + 1, chunk_iter->total_num_tchunks);
	ndim = chunk_iter->h5dset->ndim;
	for (along = 0; along < ndim; along++) {
		i = chunk_iter->tchunk_midx_buf[along] + 1;
		if (along != 0)
			Rprintf(", ");
		Rprintf("%d/%d", i, chunk_iter->num_tchunks[along]);
	}
	Rprintf("] -- <<");
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		i = chunk_iter->tchunk_midx_buf[along];
		if (GET_LIST_ELT(chunk_iter->index, along) != R_NilValue) {
			tchunkidx =
			    chunk_iter->tchunkidx_bufs->elts[along]->elts[i];
		} else {
			tchunkidx = i;
		}
		if (along != 0)
			Rprintf(", ");
		Rprintf("#%lld=%llu:%llu", tchunkidx + 1,
			chunk_iter->h5dset_vp.h5off[h5along] + 1,
			chunk_iter->h5dset_vp.h5off[h5along] +
			    chunk_iter->h5dset_vp.h5dim[h5along]);
	}
	Rprintf(">>\n");
	return;
}

/* Return 1 if the chunk that 'h5dset_vp' is pointing at is "truncated"
   (a.k.a. "partial edge chunk" in HDF5's terminology), and 0 otherwise
   (i.e. if the new chunk is a full-size chunk). */
int _tchunk_is_truncated(const H5DSetDescriptor *h5dset,
			 const H5Viewport *h5dset_vp)
{
	int ndim, h5along;
	hsize_t chunkd, d;

	ndim = h5dset->ndim;
	for (h5along = 0; h5along < ndim; h5along++) {
		chunkd = h5dset->h5chunkdim[h5along];
		d = h5dset_vp->h5dim[h5along];
		if (d != chunkd)
			return 1;
	}
	return 0;
}

int _tchunk_is_fully_selected(int ndim,
		const H5Viewport *h5dset_vp,
		const H5Viewport *mem_vp)
{
	int along, h5along, not_fully;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		not_fully = h5dset_vp->h5dim[h5along] !=
			    (hsize_t) mem_vp->dim[along];
		if (not_fully)
			return 0;
	}
	return 1;
}

void _destroy_ChunkDataBuffer(ChunkDataBuffer *chunk_data_buf)
{
	if (chunk_data_buf->data_space_id != -1)
		H5Sclose(chunk_data_buf->data_space_id);
	if (chunk_data_buf->data != NULL)
		free(chunk_data_buf->data);
	if (chunk_data_buf->data_vp.h5off != NULL)
		free(chunk_data_buf->data_vp.h5off);
	if (chunk_data_buf->compressed_data != NULL)
		free(chunk_data_buf->compressed_data);
	return;
}

int _init_ChunkDataBuffer(ChunkDataBuffer *chunk_data_buf,
		const H5DSetDescriptor *h5dset, int use_Rtype)
{
	size_t data_length, data_type_size;
	hid_t data_type_id;
	int h5along;
	const H5TypeDescriptor *h5type;

	if (h5dset->h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("'h5dset->h5chunkdim' is NULL");
		return -1;
	}

	/* Initialize ChunkDataBuffer struct members that control
	   what _destroy_ChunkDataBuffer() needs to free or close. */
	chunk_data_buf->data_space_id = -1;
	chunk_data_buf->data = NULL;
	chunk_data_buf->data_vp.h5off = NULL;
	chunk_data_buf->compressed_data = NULL;

	/* Set struct member 'data_length'. */
	data_length = 1;
	for (h5along = 0; h5along < h5dset->ndim; h5along++)
		data_length *= h5dset->h5chunkdim[h5along];
	chunk_data_buf->data_length = data_length;

	/* Set struct members 'data_type_id' and 'data_type_size'. */
	h5type = h5dset->h5type;
	if (h5type->h5class == H5T_STRING) {
		data_type_id = h5type->h5type_id;
		data_type_size = h5type->h5type_size;
	} else if (use_Rtype ||
		   h5type->native_type_size >= h5type->Rtype_size)
	{
		/* Copying data from 'chunk_data_buf' to final R array will
		   require NO type casting. */
		data_type_id = h5type->native_type_id_for_Rtype;
		data_type_size = h5type->Rtype_size;
	} else {
		/* Copying data from 'chunk_data_buf' to final R array will
		   require type casting. */
		data_type_id = h5type->native_type_id;
		data_type_size = h5type->native_type_size;
	}
	chunk_data_buf->data_type_id = data_type_id;
	chunk_data_buf->data_type_size = data_type_size;

	/* Set struct member 'data_size'. */
	chunk_data_buf->data_size = data_length * data_type_size;
	return 0;
}

int _load_chunk(const ChunkIterator *chunk_iter,
		ChunkDataBuffer *chunk_data_buf,
		int use_H5Dread_chunk)
{
	const H5DSetDescriptor *h5dset;
	hid_t data_space_id;
	int ret;

	if (chunk_data_buf->data == NULL) {
		chunk_data_buf->data = malloc(chunk_data_buf->data_size);
		if (chunk_data_buf->data == NULL) {
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for 'chunk_data_buf->data'");
			return -1;
		}
	}
	h5dset = chunk_iter->h5dset;
	if (!use_H5Dread_chunk) {
		if (chunk_data_buf->data_space_id == -1) {
			data_space_id = H5Screate_simple(h5dset->ndim,
							 h5dset->h5chunkdim,
							 NULL);
			if (data_space_id < 0) {
				PRINT_TO_ERRMSG_BUF("H5Screate_simple() "
						    "returned an error");
				return -1;
			}
			chunk_data_buf->data_space_id = data_space_id;
		}
		if (chunk_data_buf->data_vp.h5off == NULL) {
			chunk_data_buf->data_vp.h5off =
				_alloc_hsize_t_buf(h5dset->ndim, 1,
					"'chunk_data_buf->data_vp.h5off'");
			if (chunk_data_buf->data_vp.h5off == NULL)
				return -1;
		}
		chunk_data_buf->data_vp.h5dim = chunk_iter->h5dset_vp.h5dim;
		ret = _read_H5Viewport(h5dset,
				&chunk_iter->h5dset_vp,
				chunk_data_buf->data_type_id,
				chunk_data_buf->data_space_id,
				chunk_data_buf->data,
				&chunk_data_buf->data_vp);
	} else {
		/* Experimental! */
		if (chunk_data_buf->compressed_data == NULL) {
			chunk_data_buf->compressed_data =
				malloc(chunk_data_buf->data_size +
				       CHUNK_COMPRESSION_OVERHEAD);
			if (chunk_data_buf->compressed_data == NULL) {
				PRINT_TO_ERRMSG_BUF(
					"failed to allocate memory for "
					"'chunk_data_buf->compressed_data'");
				return -1;
			}
		}
		ret = read_h5chunk(h5dset->dset_id,
				   &chunk_iter->h5dset_vp,
				   chunk_data_buf);
	}
	return ret;
}

int _reclaim_vlen_bufs(ChunkDataBuffer *chunk_data_buf)
{
	int ret;

	ret = H5Dvlen_reclaim(chunk_data_buf->data_type_id,
			      chunk_data_buf->data_space_id,
			      H5P_DEFAULT, chunk_data_buf->data);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dvlen_reclaim() returned an error");
	return ret;
}

