/****************************************************************************
 *            Exploring alternate rhdf5::h5read() implementations           *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************
 *
 * Some useful links:
 * - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5S.html
 * - Documentation of H5Dread():
 *     https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
 * - An H5Dread() example:
 *     https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
 */
#include "HDF5Array.h"

#include <zlib.h>

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcmp */

//#include <time.h>

typedef struct {
	hsize_t *h5off, *h5dim;
	int *off, *dim; // same as h5off and h5dim but stored as int
} Viewport;


/****************************************************************************
 * Low-level helpers
 */

/* 'mode' controls what fields we should allocate memory for:
 *     mode == 0: 'h5off' and 'h5dim' fields only
 *     mode == 1: 'off' and 'dim' fields only
 *     mode == 2: for all the fields
 */
static int alloc_Viewport(Viewport *vp, int ndim, int mode)
{
	vp->h5off = NULL;
	vp->off = NULL;
	if (mode != 1) {
		/* Allocate memory for the 'h5off' and 'h5dim' fields. */
		vp->h5off = _alloc_hsize_t_buf(2 * ndim, 0, "Viewport fields");
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
					    "for Viewport fields");
			return -1;
		}
		vp->dim = vp->off + ndim;
	}
	return 0;
}

static void free_Viewport(Viewport *vp)
{
	if (vp->h5off != NULL)
		free(vp->h5off);
	if (vp->off != NULL)
		free(vp->off);
	return;
}

static int set_hyperslab(hid_t space_id, const Viewport *vp)
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

static int add_hyperslab(hid_t space_id, const Viewport *vp)
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

static hid_t get_mem_space(int ndim, const int *ans_dim)
{
	hsize_t *h5dim;
	int along, h5along;
	hid_t mem_space_id;

	/* Allocate and set 'h5dim'. */
	h5dim = _alloc_hsize_t_buf(ndim, 0, "'h5dim'");
	if (h5dim == NULL)
		return -1;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		h5dim[h5along] = ans_dim[along];

	mem_space_id = H5Screate_simple(ndim, h5dim, NULL);
	if (mem_space_id < 0)
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
	free(h5dim);
	return mem_space_id;
}

static long long int check_selection_against_dset(
		const DSetHandle *dset_handle,
		SEXP starts, SEXP counts, int *selection_dim_buf)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = dset_handle->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) dset_handle->h5dim[h5along];
	return _check_selection(starts, counts, dim_buf->elts,
				selection_dim_buf);
}

static long long int check_ordered_selection_against_dset(
		const DSetHandle *dset_handle,
		SEXP starts, SEXP counts, int *selection_dim_buf,
		int *nstart_buf, int *nblock_buf,
		long long int *last_block_start_buf)
{
	int ndim, along, h5along;
	LLongAE *dim_buf;

	ndim = dset_handle->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--)
		dim_buf->elts[along] =
			(long long int) dset_handle->h5dim[h5along];
	return _check_ordered_selection(starts, counts, dim_buf->elts,
					selection_dim_buf,
					nstart_buf, nblock_buf,
					last_block_start_buf);
}

static int map_starts_to_chunks(const DSetHandle *dset_handle,
				SEXP starts,
				int *nstart_buf,
				IntAEAE *breakpoint_bufs,
				LLongAEAE *chunkidx_bufs)
{
	int ndim, along, h5along;
	LLongAE *dim_buf, *chunkdim_buf;

	ndim = dset_handle->ndim;
	dim_buf = new_LLongAE(ndim, ndim, 0);
	chunkdim_buf = new_LLongAE(ndim, ndim, 0);
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		dim_buf->elts[along] =
			(long long int) dset_handle->h5dim[h5along];
		chunkdim_buf->elts[along] =
			(long long int) dset_handle->h5chunk_spacings[h5along];
	}
	return _map_starts_to_chunks(starts, dim_buf->elts, chunkdim_buf->elts,
				     nstart_buf,
				     breakpoint_bufs, chunkidx_bufs);
}

static inline int next_midx(int ndim, const int *max_idx_plus_one,
			    int *midx_buf)
{
	int along, i;

	for (along = 0; along < ndim; along++) {
		i = midx_buf[along] + 1;
		if (i < max_idx_plus_one[along]) {
			midx_buf[along] = i;
			break;
		}
		midx_buf[along] = 0;
	}
	return along;
}


/****************************************************************************
 * read_data_1_2()
 *
 * A single call to H5Dread().
 *
 * More precisely:
 *   - First walk over the blocks described by 'starts' and 'counts' and
 *     add each block to the selection.
 *   - Then make a single call to H5Dread().
 */

static size_t set_nblock_buf(int ndim, SEXP starts,
			     int expand, const int *ans_dim, int *nblock_buf)
{
	size_t num_blocks;
	int along, nblock;
	SEXP start;

	num_blocks = 1;
	for (along = 0; along < ndim; along++) {
		start = VECTOR_ELT(starts, along);
		if (start != R_NilValue) {
			nblock = LENGTH(start);
		} else {
			nblock = expand ? ans_dim[along] : 1;
		}
		num_blocks *= nblock_buf[along] = nblock;
	}
	return num_blocks;
}

static void init_srcvp_buf(const DSetHandle *dset_handle,
			   SEXP starts, Viewport *srcvp_buf)
{
	int ndim, along, h5along;
	hsize_t d;

	ndim = dset_handle->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			srcvp_buf->h5off[h5along] = 0;
			d = dset_handle->h5dim[h5along];
		} else {
			d = 1;
		}
		srcvp_buf->h5dim[h5along] = d;
	}
	return;
}

static void update_srcvp_buf(int ndim,
			const int *midx, int moved_along,
			SEXP starts, SEXP counts,
			Viewport *srcvp_buf)
{
	int along, h5along, i;
	SEXP start, count;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		srcvp_buf->h5off[h5along] = _get_trusted_elt(start, i) - 1;
		if (counts == R_NilValue)
			continue;
		count = VECTOR_ELT(counts, along);
		if (count == R_NilValue)
			continue;
		srcvp_buf->h5dim[h5along] = _get_trusted_elt(count, i);
	}
	return;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_hyperslabs(const DSetHandle *dset_handle,
			SEXP starts, SEXP counts, const int *ans_dim,
			int *nblock_buf, int *midx_buf)
{
	int ret, ndim, moved_along;
	Viewport srcvp_buf;
	long long int num_hyperslabs;

	ret = H5Sselect_none(dset_handle->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = dset_handle->ndim;
	set_nblock_buf(ndim, starts, 0, ans_dim, nblock_buf);

	/* Allocate 'srcvp_buf'. */
	if (alloc_Viewport(&srcvp_buf, ndim, 0) < 0)
		return -1;

	init_srcvp_buf(dset_handle, starts, &srcvp_buf);

	/* Walk on the hyperslabs. */
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		update_srcvp_buf(ndim, midx_buf, moved_along,
				  starts, counts, &srcvp_buf);
		/* Add to current selection. */
		ret = add_hyperslab(dset_handle->space_id, &srcvp_buf);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nblock_buf, midx_buf);
	} while (moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);

	free_Viewport(&srcvp_buf);
	return ret < 0 ? -1 : num_hyperslabs;
}

static inline hsize_t *add_element(int ndim, const int *outer_midx,
				   SEXP starts, hsize_t *coord_p)
{
	int h5along, i;
	SEXP start;
	long long int coord;

	for (h5along = ndim - 1; h5along >= 0; h5along--) {
		i = outer_midx[h5along];
		start = VECTOR_ELT(starts, h5along);
		if (start == R_NilValue) {
			coord = i;
		} else {
			coord = _get_trusted_elt(start, i) - 1;
		}
		*(coord_p++) = (hsize_t) coord;
	}
	return coord_p;
}

/* Return nb of selected elements (or -1 on error). */
static long long int select_elements(const DSetHandle *dset_handle,
			SEXP starts, const int *ans_dim,
			int *nblock_buf, int *midx_buf)
{
	int ndim, outer_moved_along, ret;
	size_t num_elements;
	hsize_t *coord_buf, *coord_p;

	ndim = dset_handle->ndim;
	num_elements = set_nblock_buf(ndim, starts, 1, ans_dim, nblock_buf);

	/* Allocate 'coord_buf'. */
	coord_buf = _alloc_hsize_t_buf(num_elements * ndim, 0, "'coord_buf'");
	if (coord_buf == NULL)
		return -1;

	/* Walk on the selected elements. */
	coord_p = coord_buf;
	outer_moved_along = ndim;
	do {
		coord_p = add_element(ndim, midx_buf, starts, coord_p);
		outer_moved_along = next_midx(ndim, nblock_buf, midx_buf);
	} while (outer_moved_along < ndim);

	ret = H5Sselect_elements(dset_handle->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	free(coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

static int set_selection(const DSetHandle *dset_handle, int method,
			 SEXP starts, SEXP counts, const int *ans_dim)
{
	int ndim;
	IntAE *nblock_buf, *midx_buf;
	long long int num_blocks;

	ndim = dset_handle->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	//clock_t t0 = clock();
	if (method == 1) {
		num_blocks = select_hyperslabs(dset_handle,
					starts, counts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	} else {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "'method' is set to 2");
			return -1;
		}
		num_blocks = select_elements(dset_handle,
					starts, ans_dim,
					nblock_buf->elts, midx_buf->elts);
	}
	//double dt = (1.0 * clock() - t0) / CLOCKS_PER_SEC;
	//printf("time for setting selection: %e\n", dt);
	//printf("num_blocks: %lld, time per block: %e\n",
	//	num_blocks, dt / num_blocks);
	return num_blocks < 0 ? -1 : 0;
}

static int read_data_1_2(const DSetHandle *dset_handle, int method,
			 SEXP starts, SEXP counts, const int *ans_dim,
			 void *out, hid_t mem_space_id)
{
	int ret;

	ret = set_selection(dset_handle, method,
			    starts, counts, ans_dim);
	if (ret < 0)
		return -1;

	ret = H5Sselect_all(mem_space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_all() returned an error");
		return -1;
	}

	//clock_t t0 = clock();
	ret = H5Dread(dset_handle->dset_id,
		      dset_handle->mem_type_id, mem_space_id,
		      dset_handle->space_id, H5P_DEFAULT, out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	//printf("time for reading data from selection: %e\n",
	//	(1.0 * clock() - t0) / CLOCKS_PER_SEC);
	return ret;
}


/****************************************************************************
 * read_data_3()
 *
 * One call to H5Dread() per block in the selection.
 *
 * More precisely, walk over the blocks described by 'starts' and 'counts'
 * and make one call to H5Dread() per block.
 *
 */

static int read_hyperslab(const DSetHandle *dset_handle,
		SEXP starts, SEXP counts,
		const int *midx, int moved_along,
		Viewport *srcvp_buf, Viewport *destvp_buf,
		void *out, hid_t mem_space_id)
{
	int ndim, along, h5along, i, ret;
	SEXP start;

	ndim = dset_handle->ndim;

	/* Update 'destvp_buf' and 'srcvp_buf' IN THAT ORDER! */
	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		h5along = ndim - 1 - along;
		i = midx[along];
		if (i == 0) {
			destvp_buf->h5off[h5along] = 0;
		} else {
			destvp_buf->h5off[h5along] += srcvp_buf->h5dim[h5along];
		}
	}
	update_srcvp_buf(ndim, midx, moved_along, starts, counts, srcvp_buf);

	ret = set_hyperslab(dset_handle->space_id, srcvp_buf);
	if (ret < 0)
		return -1;
	ret = set_hyperslab(mem_space_id, destvp_buf);
	if (ret < 0)
		return -1;
	ret = H5Dread(dset_handle->dset_id,
		      dset_handle->mem_type_id, mem_space_id,
		      dset_handle->space_id, H5P_DEFAULT, out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int read_data_3(const DSetHandle *dset_handle,
		       SEXP starts, SEXP counts, const int *ans_dim,
		       void *out, hid_t mem_space_id)
{
	int ndim, moved_along, ret;
	Viewport srcvp_buf, destvp_buf;
	IntAE *nblock_buf, *midx_buf;
	long long int num_hyperslabs;

	ndim = dset_handle->ndim;
	nblock_buf = new_IntAE(ndim, ndim, 0);
	midx_buf = new_IntAE(ndim, ndim, 0);

	set_nblock_buf(ndim, starts, 0, ans_dim, nblock_buf->elts);

	/* Allocate 'srcvp_buf' and 'destvp_buf'. */
	if (alloc_Viewport(&srcvp_buf, ndim, 0) < 0)
		return -1;
	destvp_buf.h5off = _alloc_hsize_t_buf(ndim, 1, "'destvp_buf.h5off'");
	destvp_buf.h5dim = srcvp_buf.h5dim;
	if (destvp_buf.h5off == NULL) {
		free_Viewport(&srcvp_buf);
		return -1;
	}

	/* Initialize 'srcvp_buf' (this also initializes 'destvp_buf.h5dim'). */
	init_srcvp_buf(dset_handle, starts, &srcvp_buf);

	/* Walk on the hyperslabs. */
	num_hyperslabs = 0;
	moved_along = ndim;
	do {
		num_hyperslabs++;
		ret = read_hyperslab(dset_handle, starts, counts,
				     midx_buf->elts, moved_along,
				     &srcvp_buf, &destvp_buf,
				     out, mem_space_id);
		if (ret < 0)
			break;
		moved_along = next_midx(ndim, nblock_buf->elts, midx_buf->elts);
	} while (moved_along < ndim);

	//printf("nb of hyperslabs = %lld\n", num_hyperslabs);
	free(destvp_buf.h5off);
	free_Viewport(&srcvp_buf);
	return ret;
}


/****************************************************************************
 * read_data_4_5()
 *
 * One call to H5Dread() per chunk touched by the selection.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Make one call to H5Dread() to load the **entire** chunk data to an
 *     intermediate buffer.
 *   - Copy the selected data from the intermediate buffer to the output
 *     array.
 *
 * Assumes that 'dset_handle->h5chunk_spacings' and 'dset_handle->h5nchunk'
 * are NOT NULL. This is NOT checked!
 */

static void set_nchunk_buf(const DSetHandle *dset_handle,
			   const SEXP starts, const LLongAEAE *chunkidx_bufs,
			   IntAE *nchunk_buf)
{
	int ndim, along, h5along, nchunk;

	ndim = dset_handle->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) != R_NilValue) {
			nchunk = LLongAE_get_nelt(chunkidx_bufs->elts[along]);
		} else {
			nchunk = dset_handle->h5nchunk[h5along];
		}
		nchunk_buf->elts[along] = nchunk;
	}
	return;
}

static void update_chunkvp_buf(const DSetHandle *dset_handle,
			const int *midx, int moved_along,
			SEXP starts, const LLongAEAE *chunkidx_bufs,
			Viewport *chunkvp_buf)
{
	int ndim, along, h5along, i;
	long long int chunkidx;
	hsize_t spacing, off, d;

	ndim = dset_handle->ndim;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = midx[along];
		if (VECTOR_ELT(starts, along) != R_NilValue) {
			chunkidx = chunkidx_bufs->elts[along]->elts[i];
		} else {
			chunkidx = i;
		}
		spacing = dset_handle->h5chunk_spacings[h5along];
		off = chunkidx * spacing;
		d = dset_handle->h5dim[h5along] - off;
		if (d > spacing)
			d = spacing;
		chunkvp_buf->h5off[h5along] = off;
		chunkvp_buf->h5dim[h5along] = d;
	}
	//printf("# chunkvp_buf->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//	printf(" %llu", chunkvp_buf->h5off[h5along]);
	//printf("\n");
	//printf("# chunkvp_buf->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//	printf(" %llu", chunkvp_buf->h5dim[h5along]);
	//printf("\n");
	return;
}

static void update_destvp_buf(int ndim,
			const int *midx, int moved_along,
			SEXP starts, const IntAEAE *breakpoint_bufs,
			const Viewport *chunkvp, Viewport *destvp_buf)
{
	int along, h5along, i, off, d;
	const int *breakpoint;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		i = midx[along];
		if (VECTOR_ELT(starts, along) != R_NilValue ) {
			breakpoint = breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = chunkvp->h5off[h5along];
			d = chunkvp->h5dim[h5along];
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
	//	printf(" %d", destvp_buf->off[along]);
	//printf("\n");
	//printf("# destvp_buf (dims):");
	//for (along = 0; along < ndim; along++)
	//	printf(" %d", destvp_buf->dim[along]);
	//printf("\n");
	return;
}

/* It takes about 218s on my laptop to load all the chunks from the EH1040
 * dataset (big 10x Genomics brain dataset in dense format, chunks of 100x100,
 * wrapped in the TENxBrainData package). That's 60 microseconds per chunk! */
static int load_chunk(const DSetHandle *dset_handle,
		      const Viewport *chunkvp,
		      const Viewport *middlevp,
		      void *chunk_data_out, hid_t mem_space_id)
{
	int ret;

	ret = set_hyperslab(dset_handle->space_id, chunkvp);
	if (ret < 0)
		return -1;
	ret = set_hyperslab(mem_space_id, middlevp);
	if (ret < 0)
		return -1;

	ret = H5Dread(dset_handle->dset_id,
		      dset_handle->mem_type_id, mem_space_id,
		      dset_handle->space_id, H5P_DEFAULT, chunk_data_out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

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
#define COMPRESSION_OVERHEAD 8	// empirical (increase if necessary)

static int direct_load_chunk(const DSetHandle *dset_handle,
			     const Viewport *chunkvp,
			     void *chunk_data_out,
			     void *compressed_chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(dset_handle->dset_id,
					chunkvp->h5off,
					&chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > dset_handle->chunk_data_buf_size +
				 COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%lu + %d)",
				    chunk_storage_size,
				    dset_handle->chunk_data_buf_size,
				    COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(dset_handle->dset_id, H5P_DEFAULT,
			    chunkvp->h5off, &filters,
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
				     dset_handle->chunk_data_buf_size);
}

static void init_in_offset_and_out_offset(int ndim, SEXP starts,
			const int *outdim, const Viewport *destvp,
			const Viewport *chunkvp,
			const hsize_t *h5chunk_spacings,
			size_t *in_offset, size_t *out_offset)
{
	size_t in_off, out_off;
	int along, h5along, i;
	SEXP start;

	in_off = out_off = 0;
	for (along = ndim - 1, h5along = 0; along >= 0; along--, h5along++) {
		in_off *= h5chunk_spacings[h5along];
		out_off *= outdim[along];
		i = destvp->off[along];
		start = VECTOR_ELT(starts, along);
		if (start != R_NilValue)
			in_off += _get_trusted_elt(start, i) - 1 -
				  chunkvp->h5off[h5along];
		out_off += i;
	}
	*in_offset = in_off;
	*out_offset = out_off;
	//printf("# in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static inline void update_in_offset_and_out_offset(int ndim,
			const int *inner_midx, int inner_moved_along,
			SEXP starts,
			const int *outdim, const Viewport *destvp,
			const hsize_t *h5chunk_spacings,
			size_t *in_offset, size_t *out_offset)
{
	SEXP start;
	int i1, i0, along, h5along, di;
	long long int in_off_inc, out_off_inc;

	start = VECTOR_ELT(starts, inner_moved_along);
	if (start != R_NilValue) {
		i1 = destvp->off[inner_moved_along] +
		     inner_midx[inner_moved_along];
		i0 = i1 - 1;
		in_off_inc = _get_trusted_elt(start, i1) -
			     _get_trusted_elt(start, i0);
	} else {
		in_off_inc = 1;
	}
	out_off_inc = 1;
	if (inner_moved_along >= 1) {
		along = inner_moved_along - 1;
		h5along = ndim - inner_moved_along;
		do {
			in_off_inc *= h5chunk_spacings[h5along];
			out_off_inc *= outdim[along];
			di = 1 - destvp->dim[along];
			start = VECTOR_ELT(starts, along);
			if (start != R_NilValue) {
				i1 = destvp->off[along];
				i0 = i1 - di;
				in_off_inc += _get_trusted_elt(start, i1) -
					      _get_trusted_elt(start, i0);
			} else {
				in_off_inc += di;
			}
			out_off_inc += di;
			along--;
			h5along++;
		} while (along >= 0);
	}
	*in_offset += in_off_inc;
	*out_offset += out_off_inc;
	//printf("## in_offset = %lu out_offset = %lu\n",
	//       *in_offset, *out_offset);
	return;
}

static int gather_chunk_data(const DSetHandle *dset_handle,
			const int *outer_midx, int outer_moved_along,
			SEXP starts, const IntAEAE *breakpoint_bufs,
			const void *in,
			const Viewport *chunkvp,
			SEXP ans, const int *outdim,
			const Viewport *destvp,
			int *inner_midx_buf,
			hid_t mem_type_id)
{
	int ndim, inner_moved_along;
	size_t in_offset, out_offset, s_len;
	long long int num_elts;
	const char *s;
	SEXP ans_elt;

	ndim = dset_handle->ndim;
	init_in_offset_and_out_offset(ndim, starts,
			outdim, destvp,
			chunkvp, dset_handle->h5chunk_spacings,
			&in_offset, &out_offset);

	/* Walk on the selected elements in current chunk. */
	num_elts = 0;
	while (1) {
		num_elts++;
		switch (dset_handle->Rtype) {
		    case LGLSXP:
			LOGICAL(ans)[out_offset] = ((int *) in)[in_offset];
		    break;
		    case INTSXP:
			INTEGER(ans)[out_offset] = ((int *) in)[in_offset];
		    break;
		    case REALSXP:
			REAL(ans)[out_offset] = ((double *) in)[in_offset];
		    break;
		    case STRSXP:
			s = (char *) in + in_offset * dset_handle->size;
			for (s_len = 0; s_len < dset_handle->size; s_len++)
				if (s[s_len] == 0)
					break;
			ans_elt = PROTECT(mkCharLen(s, s_len));
			SET_STRING_ELT(ans, out_offset, ans_elt);
			UNPROTECT(1);
		    break;
		    default:
			PRINT_TO_ERRMSG_BUF("unsupported type: %s",
					    CHAR(type2str(dset_handle->Rtype)));
			return -1;
		}
		inner_moved_along = next_midx(ndim, destvp->dim,
					      inner_midx_buf);
		if (inner_moved_along == ndim)
			break;
		update_in_offset_and_out_offset(ndim,
				inner_midx_buf, inner_moved_along,
				starts,
				outdim, destvp,
				dset_handle->h5chunk_spacings,
				&in_offset, &out_offset);
	};
	//printf("# nb of selected elements in current chunk = %lld\n",
	//       num_elts);
	return 0;
}

/*
  WARNING: method 5 is not working properly on some datasets:
      library(HDF5Array)
      library(ExperimentHub)
      hub <- ExperimentHub()
      fname0 <- hub[["EH1039"]]
      h5mread(fname0, "mm10/barcodes", list(1), method=4L)
      # [1] "AAACCTGAGATAGGAG-1"
      h5mread(fname0, "mm10/barcodes", list(1), method=5L)
      # [1] "AAAAAAAAAAAAAAAAAAAA"
  Looks like the chunk data has been shuffled (transposed in that case)
  before being written to disk in order to improve compression.
  TODO: Investigate this further. I suspect we need to check whether a
  "Data shuffling filter" (H5Z_FILTER_SHUFFLE) was used at creation time.
  Check H5Pget_filter() for how to know whether this filter was used or not.
  There should be a way to retrieve information about how the data was
  shuffled.
 */
static int read_data_4_5(const DSetHandle *dset_handle, int method,
			 SEXP starts,
			 const IntAEAE *breakpoint_bufs,
			 const LLongAEAE *chunkidx_bufs,
			 const int *outdim, SEXP ans)
{
	int ndim, moved_along, ret;
	hid_t mem_space_id;
	void *chunk_data_buf, *compressed_chunk_data_buf = NULL;
	Viewport chunkvp_buf, middlevp_buf, destvp_buf;
	IntAE *nchunk_buf, *outer_midx_buf, *inner_midx_buf;
	long long int num_chunks, ndot;

	ndim = dset_handle->ndim;
	mem_space_id = H5Screate_simple(ndim, dset_handle->h5chunk_spacings,
					NULL);
	if (mem_space_id < 0) {
		PRINT_TO_ERRMSG_BUF("H5Screate_simple() returned an error");
		return -1;
	}

	/* Allocate 'chunkvp_buf', 'middlevp_buf', and 'destvp_buf'. */
	if (alloc_Viewport(&chunkvp_buf, ndim, 0) < 0) {
		H5Sclose(mem_space_id);
		return -1;
	}
	middlevp_buf.h5off = _alloc_hsize_t_buf(ndim, 1,
						"'middlevp_buf.h5off'");
	middlevp_buf.h5dim = chunkvp_buf.h5dim;
	if (middlevp_buf.h5off == NULL) {
		free_Viewport(&chunkvp_buf);
		H5Sclose(mem_space_id);
		return -1;
	}
	if (alloc_Viewport(&destvp_buf, ndim, 1) < 0) {
		free(middlevp_buf.h5off);
		free_Viewport(&chunkvp_buf);
		H5Sclose(mem_space_id);
		return -1;
	}

	if (method == 4) {
		chunk_data_buf = malloc(dset_handle->chunk_data_buf_size);
	} else {
		warning("method 5 is still experimental, use at your own risk");
		chunk_data_buf = malloc(2 * dset_handle->chunk_data_buf_size +
					COMPRESSION_OVERHEAD);
	}
	if (chunk_data_buf == NULL) {
		free_Viewport(&destvp_buf);
		free(middlevp_buf.h5off);
		free_Viewport(&chunkvp_buf);
		H5Sclose(mem_space_id);
		PRINT_TO_ERRMSG_BUF("failed to allocate memory "
				    "for 'chunk_data_buf'");
		return -1;
	}
	if (method != 4)
		compressed_chunk_data_buf = chunk_data_buf +
					    dset_handle->chunk_data_buf_size;

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(dset_handle, starts, chunkidx_bufs, nchunk_buf);
	outer_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);

	/* Walk on the chunks. */
	num_chunks = ndot = 0;
	moved_along = ndim;
	//clock_t t_load_chunk = 0, t_gather_chunk_data = 0, t0;
	do {
		num_chunks++;
		//if (num_chunks % 20000 == 0) {
		//	printf(".");
		//	if (++ndot % 50 == 0)
		//		printf(" [%lld chunks processed]\n",
		//		       num_chunks);
		//	fflush(stdout);
		//}
		update_chunkvp_buf(dset_handle,
			outer_midx_buf->elts, moved_along,
			starts, chunkidx_bufs, &chunkvp_buf);

		//t0 = clock();
		if (method == 4) {
			ret = load_chunk(dset_handle,
					 &chunkvp_buf, &middlevp_buf,
					 chunk_data_buf, mem_space_id);
		} else {
			ret = direct_load_chunk(dset_handle,
					&chunkvp_buf,
					chunk_data_buf,
					compressed_chunk_data_buf);
		}
		if (ret < 0)
			break;
		//t_load_chunk += clock() - t0;

		//t0 = clock();
		update_destvp_buf(ndim,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			&chunkvp_buf, &destvp_buf);

		ret = gather_chunk_data(dset_handle,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			chunk_data_buf,
			&chunkvp_buf,
			ans, outdim, &destvp_buf,
			inner_midx_buf->elts,
			dset_handle->mem_type_id);
		if (ret < 0)
			break;
		//t_gather_chunk_data += clock() - t0;

		moved_along = next_midx(ndim, nchunk_buf->elts,
					outer_midx_buf->elts);
	} while (moved_along < ndim);
	//printf("\n");
	//printf("total chunks processed = %lld\n", num_chunks);
	//double dt1 = 1.0 * t_load_chunk / CLOCKS_PER_SEC;
	//printf("total time spent loading chunks: %3.6f\n", dt1);
	//double dt2 = 1.0 * t_gather_chunk_data / CLOCKS_PER_SEC;
	//printf("total time spent gathering the chunk data: %3.6f\n", dt2);
	//printf("total time: %3.6f\n", dt1 + dt2);

	free(chunk_data_buf);
	free_Viewport(&destvp_buf);
	free(middlevp_buf.h5off);
	free_Viewport(&chunkvp_buf);
	H5Sclose(mem_space_id);
	return ret;
}


/****************************************************************************
 * read_data_6()
 *
 * One call to H5Dread() per chunk touched by the selection. No intermediate
 * buffer.
 *
 * More precisely, walk over the chunks touched by 'starts'. For each chunk:
 *   - Select the hyperslabs obtained by intersecting the selection with
 *     the current chunk.
 *   - Call H5Dread(). This loads the selected data **directly** to the
 *     output array.
 *
 * Assumes that 'dset_handle->h5chunk_spacings' and 'dset_handle->h5nchunk'
 * are NOT NULL. This is NOT checked!
 */

/*
static hsize_t *alloc_coord_buf(int ndim, const hsize_t *h5chunk_spacings)
{
	size_t max_inner_block_per_chunk;
	int along;

	max_inner_block_per_chunk = 1;
	for (along = 0; along < ndim; along++)
		max_inner_block_per_chunk *= (h5chunk_spacings[along] + 1) / 2;
	//printf("max_inner_block_per_chunk = %lu\n",
	//	 max_inner_block_per_chunk);
	return _alloc_hsize_t_buf(max_inner_block_per_chunk * ndim,
				  "'coord_buf'");
}
*/

static void update_inner_breakpoints(int ndim, int moved_along,
		SEXP starts,
		const Viewport *destvp,
		IntAEAE *inner_breakpoint_bufs, int *inner_nblock_buf)
{
	int along, d, off, nblock, i;
	IntAE *inner_breakpoint_buf;
	SEXP start;
	long long int s0, s1;

	for (along = 0; along < ndim; along++) {
		if (along > moved_along)
			break;
		inner_breakpoint_buf = inner_breakpoint_bufs->elts[along];
		IntAE_set_nelt(inner_breakpoint_buf, 0);
		d = destvp->dim[along];
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue) {
			IntAE_insert_at(inner_breakpoint_buf, 0, d);
			inner_nblock_buf[along] = 1;
			continue;
		}
		off = destvp->off[along];
		s1 = _get_trusted_elt(start, off);
		nblock = 0;
		for (i = 1; i < d; i++) {
			s0 = s1;
			s1 = _get_trusted_elt(start, off + i);
			if (s1 != s0 + 1)
				IntAE_insert_at(inner_breakpoint_buf,
						nblock++, i);
		}
		IntAE_insert_at(inner_breakpoint_buf, nblock++, d);
		inner_nblock_buf[along] = nblock;
	}
	return;
}

static void init_inner_srcvp_buf(int ndim,
			SEXP starts, const Viewport *chunkvp,
			Viewport *inner_srcvp_buf)
{
	int along, h5along;
	hsize_t d;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (VECTOR_ELT(starts, along) == R_NilValue) {
			inner_srcvp_buf->h5off[h5along] =
					chunkvp->h5off[h5along];
			d = chunkvp->h5dim[h5along];
		} else {
			d = 1;
		}
		inner_srcvp_buf->h5dim[h5along] = d;
	}
	return;
}

static void update_inner_srcvp_buf(int ndim,
			SEXP starts, const Viewport *destvp,
			const int *inner_midx, int inner_moved_along,
			const IntAEAE *inner_breakpoint_bufs,
			Viewport *inner_srcvp_buf)
{
	int along, h5along, idx, off, d, i;
	SEXP start;
	const int *inner_breakpoint;

	for (along = 0; along < ndim; along++) {
		if (along > inner_moved_along)
			break;
		start = VECTOR_ELT(starts, along);
		if (start == R_NilValue)
			continue;
		inner_breakpoint = inner_breakpoint_bufs->elts[along]->elts;
		idx = inner_midx[along];
		off = idx == 0 ? 0 : inner_breakpoint[idx - 1];
		d = inner_breakpoint[idx] - off;
		i = destvp->off[along] + off;
		h5along = ndim - 1 - along;
		inner_srcvp_buf->h5off[h5along] =
				_get_trusted_elt(start, i) - 1;
		inner_srcvp_buf->h5dim[h5along] = d;
	}
	return;
}

/* Return nb of selected elements (or -1 on error). */
static long long int select_elements_from_chunk(const DSetHandle *dset_handle,
			SEXP starts,
			const Viewport *destvp,
			int *inner_midx_buf, hsize_t *coord_buf)
{
	int ret, ndim, inner_moved_along, along, h5along, i;
	size_t num_elements;
	hsize_t *coord_p;
	long long int coord;
	SEXP start;

	ret = H5Sselect_none(dset_handle->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = dset_handle->ndim;

	/* Walk on the selected elements from the current chunk. */
	num_elements = 0;
	coord_p = coord_buf;
	inner_moved_along = ndim;
	do {
		num_elements++;
		for (along = ndim - 1, h5along = 0;
		     along >= 0;
		     along--, h5along++)
		{
			i = destvp->h5off[h5along] + inner_midx_buf[along];
			start = VECTOR_ELT(starts, along);
			if (start == R_NilValue) {
				coord = i;
			} else {
				coord = _get_trusted_elt(start, i) - 1;
			}
			*(coord_p++) = (hsize_t) coord;
		}
		inner_moved_along = next_midx(ndim, destvp->dim,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);

	ret = H5Sselect_elements(dset_handle->space_id, H5S_SELECT_APPEND,
				 num_elements, coord_buf);
	if (ret < 0)
		return -1;
	return (long long int) num_elements;
}

/* Return nb of hyperslabs (or -1 on error). */
static long long int select_intersection_of_blocks_with_chunk(
		const DSetHandle *dset_handle,
		SEXP starts,
		const Viewport *destvp, const Viewport *chunkvp,
		const IntAEAE *inner_breakpoint_bufs,
		const int *inner_nblock,
		int *inner_midx_buf,
		Viewport *inner_srcvp_buf)
{
	int ret, ndim, inner_moved_along;
	long long int num_hyperslabs;

	ret = H5Sselect_none(dset_handle->space_id);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Sselect_none() returned an error");
		return -1;
	}

	ndim = dset_handle->ndim;

	init_inner_srcvp_buf(ndim, starts, chunkvp, inner_srcvp_buf);

	/* Walk on the block intersections with the currrent chunk. */
	num_hyperslabs = 0;
	inner_moved_along = ndim;
	do {
		num_hyperslabs++;
		update_inner_srcvp_buf(ndim, starts, destvp,
				       inner_midx_buf, inner_moved_along,
				       inner_breakpoint_bufs,
				       inner_srcvp_buf);
		ret = add_hyperslab(dset_handle->space_id, inner_srcvp_buf);
		if (ret < 0)
			return -1;
		inner_moved_along = next_midx(ndim, inner_nblock,
					      inner_midx_buf);
	} while (inner_moved_along < ndim);
	//printf("nb of hyperslabs = %lld\n", num_hyperslabs); // = prod(inner_nblock)

	return num_hyperslabs;
}

static int read_selection(const DSetHandle *dset_handle,
			  const Viewport *destvp_buf,
			  void *out, hid_t mem_space_id)
{
	int ret;

	ret = set_hyperslab(mem_space_id, destvp_buf);
	if (ret < 0)
		return -1;

	ret = H5Dread(dset_handle->dset_id,
		      dset_handle->mem_type_id, mem_space_id,
		      dset_handle->space_id, H5P_DEFAULT, out);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dread() returned an error");
	return ret;
}

static int read_data_6(const DSetHandle *dset_handle,
		       SEXP starts,
		       const IntAEAE *breakpoint_bufs,
		       const LLongAEAE *chunkidx_bufs,
		       const int *outdim, SEXP ans)
{
	void *out;
	int ndim, moved_along, ret;
	hid_t mem_space_id;
	Viewport chunkvp_buf, inner_srcvp_buf, destvp_buf;
	//hsize_t *coord_buf;
	IntAE *nchunk_buf, *outer_midx_buf,
	      *inner_nblock_buf, *inner_midx_buf;
	IntAEAE *inner_breakpoint_bufs;
	long long int num_chunks, ndot;

	switch (dset_handle->Rtype) {
	    case LGLSXP:  out = LOGICAL(ans); break;
	    case INTSXP:  out = INTEGER(ans); break;
	    case REALSXP: out = REAL(ans);    break;
	    default:
		PRINT_TO_ERRMSG_BUF("unsupported type: %s",
				    CHAR(type2str(dset_handle->Rtype)));
		return -1;
	}
	ndim = dset_handle->ndim;
	mem_space_id = get_mem_space(ndim, outdim);
	if (mem_space_id < 0)
		return -1;

	/* Allocate 'chunkvp_buf', 'inner_srcvp_buf', and 'destvp_buf'. */
	if (alloc_Viewport(&chunkvp_buf, ndim, 0) < 0) {
		H5Sclose(mem_space_id);
		return -1;
	}
	if (alloc_Viewport(&inner_srcvp_buf, ndim, 0) < 0) {
		free_Viewport(&chunkvp_buf);
		H5Sclose(mem_space_id);
		return -1;
	}
	if (alloc_Viewport(&destvp_buf, ndim, 2) < 0) {
		free_Viewport(&inner_srcvp_buf);
		free_Viewport(&chunkvp_buf);
		H5Sclose(mem_space_id);
		return -1;
	}

	/* Allocate 'coord_buf'. */
	//coord_buf = alloc_coord_buf(ndim, dset_handle->h5chunk_spacings);
	//if (coord_buf == NULL) {
	//	free_Viewport(&destvp_buf);
	//	free_Viewport(&inner_srcvp_buf);
	//	free_Viewport(&chunkvp_buf);
	//	H5Sclose(mem_space_id);
	//	return -1;
	//}

	/* Prepare buffers. */
	nchunk_buf = new_IntAE(ndim, ndim, 0);
	set_nchunk_buf(dset_handle, starts, chunkidx_bufs, nchunk_buf);
	outer_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_nblock_buf = new_IntAE(ndim, ndim, 0);
	inner_midx_buf = new_IntAE(ndim, ndim, 0);
	inner_breakpoint_bufs = new_IntAEAE(ndim, ndim);

	/* Walk on the chunks. */
	num_chunks = ndot = 0;
	moved_along = ndim;
	//clock_t t_select_elements = 0, t_read_selection = 0, t0;
	do {
		num_chunks++;
		//if (num_chunks % 20000 == 0) {
		//	printf(".");
		//	if (++ndot % 50 == 0)
		//		printf(" [%lld chunks processed]\n",
		//		       num_chunks);
		//	fflush(stdout);
		//}
		update_chunkvp_buf(dset_handle,
			outer_midx_buf->elts, moved_along,
			starts, chunkidx_bufs,
			&chunkvp_buf);

		update_destvp_buf(ndim,
			outer_midx_buf->elts, moved_along,
			starts, breakpoint_bufs,
			&chunkvp_buf, &destvp_buf);

		update_inner_breakpoints(ndim, moved_along,
			starts, &destvp_buf,
			inner_breakpoint_bufs, inner_nblock_buf->elts);

		//t0 = clock();
		/* Having 'inner_nblock_buf->elts' identical to
		   'destvp_buf.dim' means that all the inner blocks
		   are single elements so we use select_elements_from_chunk()
		   which could be faster than
		   select_intersection_of_blocks_with_chunk() in that case.
		   NO IT'S NOT FASTER! */
		//ret = memcmp(inner_nblock_buf->elts, destvp_buf.dim,
		//	     ndim * sizeof(int));
		//printf("# chunk %lld: %s\n", num_chunks,
		//       ret == 0 ? "select_elements" : "select_hyperslabs");
		//if (ret == 0) {
		//	ret = select_elements_from_chunk(dset_handle,
		//		starts, &destvp_buf,
		//		inner_midx_buf->elts, coord_buf);
		//} else {
			ret = select_intersection_of_blocks_with_chunk(
				dset_handle, starts, &destvp_buf, &chunkvp_buf,
				inner_breakpoint_bufs, inner_nblock_buf->elts,
				inner_midx_buf->elts, &inner_srcvp_buf);
		//}
		if (ret < 0)
			break;
		//t_select_elements += clock() - t0;

		//t0 = clock();
		ret = read_selection(dset_handle, &destvp_buf,
				     out, mem_space_id);
		if (ret < 0)
			break;
		//t_read_selection += clock() - t0;

		moved_along = next_midx(ndim, nchunk_buf->elts,
					outer_midx_buf->elts);
	} while (moved_along < ndim);
	//printf("\n");
	//printf("total chunks processed = %lld\n", num_chunks);
	//double dt1 = 1.0 * t_select_elements / CLOCKS_PER_SEC;
	//printf("total time spent selecting elements: %3.6f\n", dt1);
	//double dt2 = 1.0 * t_read_selection / CLOCKS_PER_SEC;
	//printf("total time spent reading elements: %3.6f\n", dt2);
	//printf("total time: %3.6f\n", dt1 + dt2);

	//free(coord_buf);
	free_Viewport(&destvp_buf);
	free_Viewport(&inner_srcvp_buf);
	free_Viewport(&chunkvp_buf);
	H5Sclose(mem_space_id);
	return ret;
}


/****************************************************************************
 * C_h5mread()
 */

/* Return R_NilValue on error. */
static SEXP h5mread_1_2_3(const DSetHandle *dset_handle,
			  SEXP starts, SEXP counts, int noreduce,
			  int method, int *ans_dim)
{
	int ndim, ret;
	long long ans_len;
	IntAE *nstart_buf, *nblock_buf;
	LLongAE *last_block_start_buf;
	SEXP ans, reduced;
	void *out;
	hid_t mem_space_id;
	int nprotect = 0;

	ndim = dset_handle->ndim;
	if (noreduce || method == 2) {
		/* This call will populate 'ans_dim'. */
		ans_len = check_selection_against_dset(dset_handle,
					starts, counts, ans_dim);
		if (ans_len < 0)
			return R_NilValue;
	} else {
		nstart_buf = new_IntAE(ndim, ndim, 0);
		nblock_buf = new_IntAE(ndim, ndim, 0);
		last_block_start_buf = new_LLongAE(ndim, ndim, 0);
		/* This call will populate 'ans_dim', 'nblock_buf',
		   and 'last_block_start_buf'. */
		ans_len = check_ordered_selection_against_dset(dset_handle,
					starts, counts, ans_dim,
					nstart_buf->elts, nblock_buf->elts,
					last_block_start_buf->elts);
		if (ans_len < 0)
			return R_NilValue;
		if (_selection_can_be_reduced(ndim, nstart_buf->elts,
						    nblock_buf->elts)) {
			reduced = PROTECT(_reduce_selection(
						starts, counts, ans_dim,
						nblock_buf->elts,
						last_block_start_buf->elts));
			nprotect++;
			starts = VECTOR_ELT(reduced, 0);
			counts = VECTOR_ELT(reduced, 1);
		}
	}

	ans = PROTECT(allocVector(dset_handle->Rtype, (R_xlen_t) ans_len));
	nprotect++;

	if (ans_len != 0) {
		switch (dset_handle->Rtype) {
		    case LGLSXP:  out = LOGICAL(ans); break;
		    case INTSXP:  out = INTEGER(ans); break;
		    case REALSXP: out = REAL(ans);    break;
		    default:
			PRINT_TO_ERRMSG_BUF("unsupported type: %s",
					    CHAR(type2str(dset_handle->Rtype)));
			goto on_error;
		}
		mem_space_id = get_mem_space(ndim, ans_dim);
		if (mem_space_id < 0)
			goto on_error;
		if (method <= 2) {
			ret = read_data_1_2(dset_handle, method,
					starts, counts, ans_dim,
					out, mem_space_id);
		} else {
			ret = read_data_3(dset_handle,
					starts, counts, ans_dim,
					out, mem_space_id);
		}
		H5Sclose(mem_space_id);
		if (ret < 0)
			goto on_error;
	}

	UNPROTECT(nprotect);
	return ans;

    on_error:
	UNPROTECT(nprotect);
	return R_NilValue;
}

/* Return R_NilValue on error. */
static SEXP h5mread_4_5_6(const DSetHandle *dset_handle, SEXP starts,
			  int method, int *ans_dim)
{
	int ndim, ret, along;
	IntAEAE *breakpoint_bufs;
	LLongAEAE *chunkidx_bufs;
	R_xlen_t ans_len;
	SEXP ans;

	ndim = dset_handle->ndim;
	breakpoint_bufs = new_IntAEAE(ndim, ndim);
	chunkidx_bufs = new_LLongAEAE(ndim, ndim);

	/* This call will populate 'ans_dim', 'breakpoint_bufs',
	   and 'chunkidx_bufs'. */
	ret = map_starts_to_chunks(dset_handle, starts,
				   ans_dim, breakpoint_bufs, chunkidx_bufs);
	if (ret < 0)
		return R_NilValue;

	ans_len = 1;
	for (along = 0; along < ndim; along++)
		ans_len *= ans_dim[along];
	ans = PROTECT(allocVector(dset_handle->Rtype, ans_len));

	if (ans_len != 0) {
		if (method <= 5) {
			/* methods 4 and 5 */
			ret = read_data_4_5(dset_handle, method,
					starts, breakpoint_bufs, chunkidx_bufs,
					ans_dim, ans);
		} else {
			/* method 6 */
			ret = read_data_6(dset_handle,
					starts, breakpoint_bufs, chunkidx_bufs,
					ans_dim, ans);
		}
		if (ret < 0)
			goto on_error;
	}

	UNPROTECT(1);
	return ans;

    on_error:
	UNPROTECT(1);
	return R_NilValue;
}

/* Return R_NilValue on error. */
static SEXP h5mread(hid_t dset_id, SEXP starts, SEXP counts, int noreduce,
		    int as_int, int method)
{
	int ndim;
	DSetHandle dset_handle;
	SEXP ans, ans_dim;

	if (method < 0 || method > 6) {
		H5Dclose(dset_id);
		PRINT_TO_ERRMSG_BUF("'method' must be >= 0 and <= 6");
		return R_NilValue;
	}

	ndim = _shallow_check_selection(starts, counts);
	if (ndim < 0) {
		H5Dclose(dset_id);
		return R_NilValue;
	}

	/* _get_DSetHandle() will do H5Dclose(dset_id) in case of an error. */
	if (_get_DSetHandle(dset_id, as_int, 0, ndim, &dset_handle) < 0)
		return R_NilValue;

	ans = R_NilValue;

	if (dset_handle.Rtype == STRSXP) {
		/* Note that it should be easy to support contiguous string
		   data by treating it as if it was made of a single chunk. */
		if (dset_handle.h5chunk_spacings == NULL) {
			PRINT_TO_ERRMSG_BUF("contiguous (i.e. not chunked) "
					    "string data is not supported "
					    "at the moment");
			goto on_error;
		}
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF("'counts' must be NULL when "
					    "reading string data");
			goto on_error;
		}
		if (method == 0) {
			method = 4;
		} else if (method != 4 && method != 5) {
			PRINT_TO_ERRMSG_BUF("reading string data is only "
					    "supported by methods 4 and 5");
			goto on_error;
		}
	} else if (method == 0) {
		/* March 27, 2019: My early testing (from Nov 2018) seemed
		   to indicate that method 6 was a better choice over method 4
		   when the layout is chunked and 'counts' is NULL. Turns out
		   that doing more testing today seems to indicate the opposite
		   i.e. method 4 now seems to perform better than method 6 on
		   all the datasets I've tested so far, including those used by
		   Pete Hickey here:
		     https://github.com/Bioconductor/DelayedArray/issues/13
		   and those used in the examples in man/h5mread.Rd.
		   Note sure what happened between Nov 2018 and today. Did I
		   do something stupid in my early testing? Did something
		   change in Rhdf5lib?
		   Anyway thanks to Pete for providing such a useful report. */
		method = dset_handle.h5chunk_spacings != NULL &&
			 //counts == R_NilValue ? 6 : 1;
			 counts == R_NilValue ? 4 : 1;
	}

	ans_dim = PROTECT(NEW_INTEGER(ndim));

	if (method <= 3) {
		ans = h5mread_1_2_3(&dset_handle, starts, counts, noreduce,
				    method, INTEGER(ans_dim));
	} else if (counts != R_NilValue) {
		PRINT_TO_ERRMSG_BUF("'counts' must be NULL for "
				    "methods 4, 5, and 6");
	} else {
		ans = h5mread_4_5_6(&dset_handle, starts,
				    method, INTEGER(ans_dim));
	}
	if (ans != R_NilValue) {
		PROTECT(ans);
		SET_DIM(ans, ans_dim);
		UNPROTECT(1);
	}

	UNPROTECT(1);

    on_error:
	_close_DSetHandle(&dset_handle);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_h5mread(SEXP filepath, SEXP name,
	       SEXP starts, SEXP counts, SEXP noreduce,
	       SEXP as_integer, SEXP method)
{
	int noreduce0, as_int, method0;
	hid_t file_id, dset_id;
	SEXP ans;

	/* Check 'noreduce'. */
	if (!(IS_LOGICAL(noreduce) && LENGTH(noreduce) == 1))
		error("'noreduce' must be TRUE or FALSE");
	noreduce0 = LOGICAL(noreduce)[0];

	/* Check 'as_integer'. */
	if (!(IS_LOGICAL(as_integer) && LENGTH(as_integer) == 1))
		error("'as_integer' must be TRUE or FALSE");
	as_int = LOGICAL(as_integer)[0];

	/* Check 'method'. */
	if (!(IS_INTEGER(method) && LENGTH(method) == 1))
		error("'method' must be a single integer");
	method0 = INTEGER(method)[0];

	file_id = _get_file_id(filepath);
	dset_id = _get_dset_id(file_id, name, filepath);

	/* h5mread() will do H5Dclose(dset_id). */
	ans = PROTECT(h5mread(dset_id, starts, counts, noreduce0,
			      as_int, method0));
	H5Fclose(file_id);
	if (ans == R_NilValue) {
		UNPROTECT(1);
		error(_HDF5Array_errmsg_buf);
	}
	UNPROTECT(1);
	return ans;
}

