### =========================================================================
### Common operations on HDF5Matrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowSums(), colSums(), rowMeans(), colMeans()
###

.HDF5Matrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    if (!identical(dims, 1))
        stop("\"rowSums\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
    if (!is.array(x) && x@is_transposed)
        return(.HDF5Matrix_block_colSums(t(x), na.rm=na.rm, dims=dims))
    blocks <- ArrayBlocks(dim(x), max(get_block_length(type(x)), nrow(x)))
    ans <- integer(nrow(x))  # init value
    for (i in seq_along(blocks)) {
        subarray <- extract_array_block(x, blocks, i)
        subarray_rowSums <- rowSums(as.matrix(subarray), na.rm=na.rm)
        ans <- ans + subarray_rowSums
    }
    ans
}

.HDF5Matrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    if (!identical(dims, 1))
        stop("\"colSums\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
    if (!is.array(x) && x@is_transposed)
        return(.HDF5Matrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))
    blocks <- ArrayBlocks(dim(x), max(get_block_length(type(x)), nrow(x)))
    ans <- lapply(seq_along(blocks),
        function(i) {
            subarray <- extract_array_block(x, blocks, i)
            colSums(as.matrix(subarray), na.rm=na.rm)
        })
    unlist(ans)
}

setMethod("rowSums", "HDF5Matrix", .HDF5Matrix_block_rowSums)
setMethod("colSums", "HDF5Matrix", .HDF5Matrix_block_colSums)

.HDF5Matrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    if (!identical(dims, 1))
        stop("\"rowMeans\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
    if (!is.array(x) && x@is_transposed)
        return(.HDF5Matrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))
    blocks <- ArrayBlocks(dim(x), max(get_block_length(type(x)), nrow(x)))
    sums <- nvals <- numeric(nrow(x))
    for (i in seq_along(blocks)) {
        subarray <- extract_array_block(x, blocks, i)
        tmp <- as.matrix(subarray)
        subarray_sums <- rowSums(tmp, na.rm=na.rm)
        subarray_nvals <- ncol(tmp)
        if (na.rm)
            subarray_nvals <- subarray_nvals - rowSums(is.na(tmp))
        sums <- sums + subarray_sums
        nvals <- nvals + subarray_nvals
    }
    sums / nvals
}

.HDF5Matrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    if (!identical(dims, 1))
        stop("\"colMeans\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
    if (!is.array(x) && x@is_transposed)
        return(.HDF5Matrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))
    blocks <- ArrayBlocks(dim(x), max(get_block_length(type(x)), nrow(x)))
    ans <- lapply(seq_along(blocks),
        function(i) {
            subarray <- extract_array_block(x, blocks, i)
            colMeans(as.matrix(subarray), na.rm=na.rm)
        })
    unlist(ans)
}

setMethod("rowMeans", "HDF5Matrix", .HDF5Matrix_block_rowMeans)
setMethod("colMeans", "HDF5Matrix", .HDF5Matrix_block_colMeans)

