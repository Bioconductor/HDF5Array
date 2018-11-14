/****************************************************************************
 *                             C_simple_h5read()                            *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "HDF5Array.h"
#include "hdf5.h"

#include <stdlib.h>  /* for malloc, free */
#include <string.h>  /* for memcpy */

/*
h5readDataset <- function(h5dataset, index)
{
    h5spaceFile <- H5Dget_space(h5dataset)
    size <- H5Sselect_index(h5spaceFile, index)
    h5spaceMem <- H5Screate_simple(size)
    .Call("_H5Dread", h5dataset@ID, h5spaceFile@ID, h5spaceMem@ID,
                      PACKAGE="rhdf5")
}
*/

static char errmsg_buf[256];

static int check_start_and_count(int dset_rank, const hsize_t *dset_dims,
				 SEXP start, SEXP count,
				 int *total_count)
{
	int along, n, i, e, s, c;
	SEXP start_elt, count_elt;

	if (!(isVectorList(start) && LENGTH(start) == dset_rank)) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'start' must be a list with one "
			 "list element per dimension in the dataset");
		return -1;
	}
	if (!(isVectorList(count) && LENGTH(count) == dset_rank)) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'count' must be a list with one "
			 "list element per dimension in the dataset");
		return -1;
	}
	for (along = 0; along < dset_rank; along++) {
		start_elt = VECTOR_ELT(start, along);
		count_elt = VECTOR_ELT(count, along);
		if (!(IS_INTEGER(start_elt) &&
		      IS_INTEGER(count_elt) &&
		      (n = LENGTH(start_elt)) == LENGTH(count_elt)))
		{
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "'start[[%d]]' and 'count[[%d]]' must "
				 "be integer vectors of the same length",
				 along + 1, along + 1);
			return -1;
		}
		total_count[along] = e = 0;
		for (i = 0; i < n; i++) {
			s = INTEGER(start_elt)[i];
			c = INTEGER(count_elt)[i];
			if (s == NA_INTEGER || s <= e) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "start[[%d]][%d] is NA or <= 0 or < "
					 "start[[%d]][%d] + count[[%d]][%d]",
					 along + 1, i + 1,
					 along + 1, i,
					 along + 1, i);
				return -1;
			}
			if (c == NA_INTEGER || c <= 0) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "count[[%d]][%d] is NA or <= 0",
					 along + 1, i + 1);
				return -1;
			}
			e = s + c - 1;
			if (e > dset_dims[along]) {
				snprintf(errmsg_buf, sizeof(errmsg_buf),
					 "start[[%d]][%d] + count[[%d]][%d] "
					 "- 1 is greater than the "
					 "corresponding dimension in the "
					 "dataset ",
					 along + 1, i + 1, along + 1, i + 1);
				return -1;
			}
			total_count[along] += c;
		}
	}
	return 0;
}

/* Should we use H5Sselect_hyperslab() or H5Sselect_elements() for this?
   Useful links:
   - Documentation of H5Sselect_hyperslab() and H5Sselect_elements():
       https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
   - Documentation of H5Dread():
       https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Read
   - A useful example:
       https://support.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CheckAndReadExample
*/

static hid_t prepare_file_space(hid_t dset_id, SEXP start, SEXP count,
				int *total_count)
{
	hid_t file_space_id;
	int dset_rank, along, ret;
	hsize_t *dset_dims;

	file_space_id = H5Dget_space(dset_id);
	if (file_space_id < 0) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
		         "H5Dget_space() returned an error");
		return -1;
	}
	dset_rank = H5Sget_simple_extent_ndims(file_space_id);
	dset_dims = (hsize_t *) malloc(dset_rank * sizeof(hsize_t));
	if (dset_dims == NULL) {
		H5Sclose(file_space_id);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'dset_dims'");
		return -1;
	}
	if (H5Sget_simple_extent_dims(file_space_id, dset_dims, NULL) !=
	    dset_rank) {
		free(dset_dims);
		H5Sclose(file_space_id);
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "H5Sget_simple_extent_dims() returned "
			 "an unexpected value");
		return -1;
	}
	printf("dset_dims:");
	for (along = 0; along < dset_rank; along++)
		printf(" %llu", dset_dims[along]);
	printf("\n");
	ret = check_start_and_count(dset_rank, dset_dims, start, count,
				    total_count);
	if (ret < 0) {
		free(dset_dims);
		H5Sclose(file_space_id);
		return -1;
	}

	// prepare selection

	free(dset_dims);
	return file_space_id;
}

static hid_t prepare_mem_space(hid_t dset_id, SEXP start, SEXP count)
{
	hid_t mem_space_id;

	return mem_space_id;
}

/* Return NULL if error. */
static SEXP simple_h5read(hid_t dset_id, SEXP start, SEXP count)
{
	int count_len, *total_count, along, ret;
	hid_t file_space_id, mem_space_id;
	R_xlen_t ans_len;
	SEXP ans, ans_dim;

	/* Allocate 'total_count'. */
	if (!isVectorList(count)) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "'count' must be a list");
		return R_NilValue;
	}
	count_len = LENGTH(count);
	total_count = (int *) malloc(count_len * sizeof(int));
	if (total_count == NULL) {
		snprintf(errmsg_buf, sizeof(errmsg_buf),
			 "failed to allocate memory for 'total_count'");
		return R_NilValue;
	}

	/* Prepare 'file_space_id'. */
	file_space_id = prepare_file_space(dset_id, start, count, total_count);
	if (file_space_id < 0) {
		free(total_count);
		return R_NilValue;
	}

	/* Prepare 'mem_space_id'. */
	//mem_space_id = prepare_mem_space(dset_id, start, count);
	//if (mem_space_id < 0) {
	//	H5Sclose(file_space_id);
	//	free(total_count);
	//	return R_NilValue;
	//}

	/* Compute 'ans_len'. */
	ans_len = total_count[0];
	for (along = 1; along < count_len; along++)
		ans_len *= total_count[along];

	ans = PROTECT(NEW_INTEGER(ans_len));
	ans_dim = PROTECT(NEW_INTEGER(count_len));
	memcpy(INTEGER(ans_dim), total_count, count_len * sizeof(int));
	SET_DIM(ans, ans_dim);
	UNPROTECT(1);

	//ret = H5Dread(dset_id, H5T_NATIVE_INT,
	//	      mem_space_id, file_space_id,
	//	      H5P_DEFAULT, INTEGER(ans));

	//H5Sclose(mem_space_id);
	H5Sclose(file_space_id);
	free(total_count);
	UNPROTECT(1);
	//if (ret < 0) {
	//	snprintf(errmsg_buf, sizeof(errmsg_buf),
	//	         "H5Dread() returned an error");
	//	return R_NilValue;
	//}
	return ans;
}

SEXP C_simple_h5read(SEXP filepath, SEXP name, SEXP start, SEXP count)
{
	SEXP filepath0, name0, ans;
	hid_t file_id, dset_id;

	/* Check 'filepath'. */
	if (!(IS_CHARACTER(filepath) && LENGTH(filepath) == 1))
		error("'filepath' must be a single string");
	filepath0 = STRING_ELT(filepath, 0);
	if (filepath0 == NA_STRING)
		error("'filepath' cannot be NA");

	/* Check 'name'. */
	if (!(IS_CHARACTER(name) && LENGTH(name) == 1))
		error("'name' must be a single string");
	name0 = STRING_ELT(name, 0);
	if (name0 == NA_STRING)
		error("'name' cannot be NA");

	file_id = H5Fopen(CHAR(filepath0), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		error("failed to open file %s", CHAR(filepath0));
	dset_id = H5Dopen(file_id, CHAR(name0), H5P_DEFAULT);
	if (dset_id < 0) {
		H5Fclose(file_id);
		error("failed to open dataset %s from file %s",
		      CHAR(name0), CHAR(filepath0));
	}
	ans = PROTECT(simple_h5read(dset_id, start, count));
	H5Dclose(dset_id);
	H5Fclose(file_id);
	if (ans == R_NilValue) {
		UNPROTECT(1);
		error(errmsg_buf);
	}
	UNPROTECT(1);
	return ans;
}

