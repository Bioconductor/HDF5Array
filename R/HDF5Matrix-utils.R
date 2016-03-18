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

    REDUCE <- function(submatrix) rowSums(submatrix, na.rm=na.rm)
    COMBINE <- `+`
    init <- numeric(nrow(x))
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans, rownames(x))
}

.HDF5Matrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colSums")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))

    colsums_list <- colblock_APPLY(x, colSums, na.rm=na.rm,
                                   if_empty=numeric(0))
    unlist(colsums_list, recursive=FALSE)
}

setMethod("rowSums", "HDF5Matrix", .HDF5Matrix_block_rowSums)
setMethod("colSums", "HDF5Matrix", .HDF5Matrix_block_colSums)

.HDF5Matrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowMeans")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))

    REDUCE <- function(submatrix) {
        submatrix_sums <- rowSums(submatrix, na.rm=na.rm)
        submatrix_nvals <- ncol(submatrix)
        if (na.rm)
            submatrix_nvals <- submatrix_nvals - rowSums(is.na(submatrix))
        cbind(submatrix_sums, submatrix_nvals)
    }
    COMBINE <- `+`
    init <- cbind(
        numeric(nrow(x)),  # sums
        numeric(nrow(x))   # nvals
    )
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans[ , 1L] / ans[ , 2L], rownames(x))
}

.HDF5Matrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colMeans")
    if (is(x, "HDF5Array") && x@is_transposed)
        return(.HDF5Matrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))

    colmeans_list <- colblock_APPLY(x, colMeans, na.rm=na.rm,
                                    if_empty=numeric(0))
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

### Write a new HDF5 dataset to disk. Return an HDF5Matrix object that points
### to this new dataset.
.HDF5Matrix_block_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "HDF5Matrix") || is.matrix(y),
              ncol(x) == nrow(y))

    out_file <- getHDF5ArrayOutputFile()
    out_name <- getHDF5ArrayOutputName()
    on.exit(setHDF5ArrayOutputName())

    ans_type <- typeof(match.fun(type(x))(1) * match.fun(type(y))(1))
    h5createDataset(out_file, out_name, c(nrow(x), ncol(y)),
                    storage.mode=ans_type)

    colblock_APPLY(y,
        function(submatrix) x %*% submatrix,
        out_file=out_file,
        out_name=out_name
    )

    ## TODO: Investigate the possiblity to store the dimnames in the HDF5 file
    ## so the HDF5Matrix() constructor can bring them back. Then we wouldn't
    ## need to explicitely set them on 'ans' like we do below.
    ans <- HDF5Matrix(out_file, "/", out_name, type=ans_type)
    ans_rownames <- rownames(x)
    ans_colnames <- colnames(y)
    if (!(is.null(ans_rownames) && is.null(ans_colnames)))
        dimnames(ans) <- list(ans_rownames, ans_colnames)
    ans
}

setMethod("%*%", c("HDF5Matrix", "matrix"),
    function(x, y) t(t(y) %*% t(x))
)

setMethod("%*%", c("matrix", "HDF5Matrix"),
    .HDF5Matrix_block_mult_by_left_matrix
)

setMethod("%*%", c("HDF5Matrix", "HDF5Matrix"),
    function(x, y)
        stop(wmsg("multiplication of 2 HDF5Matrix objects is not supported, ",
                  "only multiplication of an ordinary matrix by an ",
                  "HDF5Matrix object at the moment"))
)

