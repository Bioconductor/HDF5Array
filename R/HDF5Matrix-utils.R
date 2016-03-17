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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###
### We only support multiplication of an ordinary matrix (typically
### small) by an HDF5Matrix object (typically big). Multiplication of 2
### HDF5Matrix objects is not supported.
###

### Return an HDF5Matrix object. This object points to its own HDF5 dataset
### stored in a new file.
.HDF5Matrix_block_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "HDF5Matrix") || is.matrix(y),
              ncol(x) == nrow(y))

    out_file <- paste0(tempfile(), ".h5")
    out_name <- sprintf("%s %s %s", "x", "%*%", "y")
    ans_type <- typeof(match.fun(type(x))(1) * match.fun(type(y))(1))
    h5createFile(out_file)
    h5createDataset(out_file, out_name, c(nrow(x), ncol(y)),
                    storage.mode=ans_type)
    offset <- 0L  # offset of current block in nb of columns
    colblock_APPLY(y,
        function(submatrix) {
            z <- x %*% submatrix
            index <- list(NULL, seq_len(ncol(z)) + offset)
            h5write(z, out_file, out_name, index=index)
            offset <<- offset + ncol(z)
        }
    )
    HDF5Matrix(out_file, "/", out_name)
}

setMethod("%*%", c("HDF5Matrix", "matrix"),
    function(x, y) t(t(y) %*% t(x))
)

setMethod("%*%", c("matrix", "HDF5Matrix"),
    .HDF5Matrix_block_mult_by_left_matrix
)

