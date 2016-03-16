### =========================================================================
### Common operations on HDF5Matrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowSums(), colSums(), rowMeans(), colMeans()
###

.normarg_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop("\"", method, "\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
}

.HDF5Matrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowSums")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_colSums(t(x), na.rm=na.rm, dims=dims))

    APPLY <- function(submatrix) rowSums(submatrix, na.rm=na.rm)
    REDUCE <- `+`
    reduced <- integer(nrow(x))
    colblock_APPLY_REDUCE(x, APPLY, REDUCE, reduced)
}

.HDF5Matrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colSums")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))

    colsums_list <- colblock_APPLY(x, colSums, na.rm=na.rm)
    unlist(colsums_list, recursive=FALSE)
}

setMethod("rowSums", "HDF5Matrix", .HDF5Matrix_block_rowSums)
setMethod("colSums", "HDF5Matrix", .HDF5Matrix_block_colSums)

.HDF5Matrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowMeans")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))

    APPLY <- function(submatrix) {
        submatrix_sums <- rowSums(submatrix, na.rm=na.rm)
        submatrix_nvals <- ncol(submatrix)
        if (na.rm)
            submatrix_nvals <- submatrix_nvals - rowSums(is.na(submatrix))
        cbind(submatrix_sums, submatrix_nvals)
    }
    REDUCE <- `+`
    reduced <- cbind(
        numeric(nrow(x)),  # sums
        numeric(nrow(x))   # nvals
    )
    reduced <- colblock_APPLY_REDUCE(x, APPLY, REDUCE, reduced)
    reduced[ , 1L] / reduced[ , 2L]
}

.HDF5Matrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colMeans")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))

    colmeans_list <- colblock_APPLY(x, colMeans, na.rm=na.rm)
    unlist(colmeans_list, recursive=FALSE)
}

setMethod("rowMeans", "HDF5Matrix", .HDF5Matrix_block_rowMeans)
setMethod("colMeans", "HDF5Matrix", .HDF5Matrix_block_colMeans)

