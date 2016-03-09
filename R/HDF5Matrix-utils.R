### =========================================================================
### Common operations on HDF5Matrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowSums(), colSums(), rowMeans(), colMeans()
###

.HDF5Matrix_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    if (!identical(dims, 1))
        stop("\"rowSums\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
    if (x@transpose)
        return(.HDF5Matrix_colSums(t(x), na.rm=na.rm, dims=dims))
    blocks <- ArrayBlocks(dim(x), max(MAX_BLOCK_LENGTH, nrow(x)))
    ans <- integer(nrow(x))  # init value
    for (i in seq_along(blocks)) {
        block <- extract_array_block(x, blocks, i)
        ans <- ans + rowSums(as.matrix(block), na.rm=na.rm)
    }
    ans
}

.HDF5Matrix_colSums <- function(x, na.rm=FALSE, dims=1)
{
    if (!identical(dims, 1))
        stop("\"colSums\" method for HDF5Matrix objects ",
             "does not support the 'dims' argument yet")
    if (x@transpose)
        return(.HDF5Matrix_rowSums(t(x), na.rm=na.rm, dims=dims))
    blocks <- ArrayBlocks(dim(x), max(MAX_BLOCK_LENGTH, nrow(x)))
    ans <- lapply(seq_along(blocks),
        function(i) {
            block <- extract_array_block(x, blocks, i)
            colSums(as.matrix(block), na.rm=na.rm)
        })
    unlist(ans)
}

setMethod("rowSums", "HDF5Matrix", .HDF5Matrix_rowSums)
setMethod("colSums", "HDF5Matrix", .HDF5Matrix_colSums)

.HDF5Matrix_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    stop("rowMeans() not ready yet")
}

.HDF5Matrix_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    stop("colMeans() not ready yet")
}

setMethod("rowMeans", "HDF5Matrix", .HDF5Matrix_rowMeans)
setMethod("colMeans", "HDF5Matrix", .HDF5Matrix_colMeans)

