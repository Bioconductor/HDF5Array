### =========================================================================
### Common operations on DelayedMatrix objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RowBinder and ColBinder objects
###
### These classes are for internal use only and are not exported.
###
 
setClass("MatrixBinder",
    representation(
        "VIRTUAL",
        seeds="list"  # List of matrix-like objects to bind. Each object
                      # is expected to satisfy the "seed contract" i.e. to
                      # support dim(), dimnames(), and subset_seed_as_array().
    ),
    prototype(
        seeds=list(new("matrix"))
    )
)

setClass("RowBinder", contains="MatrixBinder")
setClass("ColBinder", contains="MatrixBinder")

.objects_are_rbindable <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    if (!all(ndims == 2L))
        return(FALSE)
    ncols <- matrix(unlist(dims, use.names=FALSE), nrow=2L)[2L, ]
    all(ncols == ncols[[1L]])
}

.objects_are_cbindable <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    if (!all(ndims == 2L))
        return(FALSE)
    nrows <- matrix(unlist(dims, use.names=FALSE), nrow=2L)[1L, ]
    all(nrows == nrows[[1L]])
}

.validate_RowBinder <- function(x)
{
    if (length(x@seeds) == 0L)
        return(wmsg("'x@seeds' cannot be empty"))
    if (!.objects_are_rbindable(x@seeds))
        return(wmsg("'x@seeds' must be a list of matrix-like objects ",
                    "with the same number of columns"))
    TRUE
}

.validate_ColBinder <- function(x)
{
    if (length(x@seeds) == 0L)
        return(wmsg("'x@seeds' cannot be empty"))
    if (!.objects_are_cbindable(x@seeds))
        return(wmsg("'x@seeds' must be a list of matrix-like objects ",
                    "with the same number of rows"))
    TRUE
}

setValidity2("RowBinder", .validate_RowBinder)
setValidity2("ColBinder", .validate_ColBinder)

.new_RowBinder <- function(seeds) new2("RowBinder", seeds=seeds)

.new_ColBinder <- function(seeds) new2("ColBinder", seeds=seeds)

### Implement the "seed contract" i.e. dim(), dimnames(), and
### subset_seed_as_array().

.get_matrices_dims <- function(matrices)
{
    matrix(unlist(lapply(matrices, dim), use.names=FALSE), nrow=2L)
}

.get_RowBinder_dim <- function(x)
{
    dims <- .get_matrices_dims(x@seeds)
    nrow <- sum(dims[1L, ])
    ncol <- dims[2L, 1L]
    c(nrow, ncol)
}

.get_ColBinder_dim <- function(x)
{
    dims <- .get_matrices_dims(x@seeds)
    nrow <- dims[1L, 1L]
    ncol <- sum(dims[2L, ])
    c(nrow, ncol)
}

setMethod("dim", "RowBinder", .get_RowBinder_dim)
setMethod("dim", "ColBinder", .get_ColBinder_dim)

.get_RowBinder_dimnames <- function(x)
{
    combine_dimnames_along(x@seeds, .get_matrices_dims(x@seeds), 1L)
}

.get_ColBinder_dimnames <- function(x)
{
    combine_dimnames_along(x@seeds, .get_matrices_dims(x@seeds), 2L)
}

setMethod("dimnames", "RowBinder", .get_RowBinder_dimnames)
setMethod("dimnames", "ColBinder", .get_ColBinder_dimnames)

.subset_MatrixBinder_as_array <- function(x, index, bind.along)
{
    breakpoints <- cumsum(.get_matrices_dims(x@seeds)[bind.along, ])
    part_idx <- get_part_index(index[[bind.along]], breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    if (bind.along == 1L) {
        FUN <- function(i)
                 subset_seed_as_array(
                     x@seeds[[i]],
                     list(split_part_idx[[i]], index[[2L]]))
    } else {
        FUN <- function(i)
                 subset_seed_as_array(
                     x@seeds[[i]],
                     list(index[[1L]], split_part_idx[[i]]))
    }
    tmp <- lapply(seq_along(x@seeds), FUN)

    ## Recombine the matrices in 'tmp'.
    rev_idx <- get_rev_index(part_idx)
    if (bind.along == 1L) {
        do.call("rbind", tmp)[rev_idx , , drop=FALSE]
    } else {
        do.call("cbind", tmp)[ , rev_idx, drop=FALSE]
    }
}

setMethod("subset_seed_as_array", "RowBinder",
    function(seed, index) .subset_MatrixBinder_as_array(seed, index, 1L)
)
setMethod("subset_seed_as_array", "ColBinder",
    function(seed, index) .subset_MatrixBinder_as_array(seed, index, 2L)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rbind() and cbind()
###

.DelayedMatrix_rbind <- function(...)
{
    objects <- unname(list(...))
    if (!.objects_are_rbindable(objects))
        stop(wmsg("can only rbind() matrix-like objects ",
                  "with the same number of columns"))
    DelayedArray(.new_RowBinder(objects))
}

.DelayedMatrix_cbind <- function(...)
{
    objects <- unname(list(...))
    if (!.objects_are_cbindable(objects))
        stop(wmsg("can only cbind() matrix-like objects ",
                  "with the same number of rows"))
    DelayedArray(.new_ColBinder(objects))
}

setMethod("rbind", "DelayedMatrix", .DelayedMatrix_rbind)
setMethod("cbind", "DelayedMatrix", .DelayedMatrix_cbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowSums(), colSums(), rowMeans(), colMeans()
###

.normarg_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop("\"", method, "\" method for DelayedMatrix objects ",
             "does not support the 'dims' argument yet")
}

.DelayedMatrix_block_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_colSums(t(x), na.rm=na.rm, dims=dims))

    REDUCE <- function(submatrix) rowSums(submatrix, na.rm=na.rm)
    COMBINE <- `+`
    init <- numeric(nrow(x))
    ans <- colblock_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init)
    setNames(ans, rownames(x))
}

.DelayedMatrix_block_colSums <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colSums")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowSums(t(x), na.rm=na.rm, dims=dims))

    colsums_list <- colblock_APPLY(x, colSums, na.rm=na.rm,
                                   if_empty=numeric(0))
    unlist(colsums_list, recursive=FALSE)
}

setMethod("rowSums", "DelayedMatrix", .DelayedMatrix_block_rowSums)
setMethod("colSums", "DelayedMatrix", .DelayedMatrix_block_colSums)

.DelayedMatrix_block_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "rowMeans")
    if (is(x, "DelayedMatrix") && x@is_transposed)
        return(.DelayedMatrix_block_colMeans(t(x), na.rm=na.rm, dims=dims))

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

.DelayedMatrix_block_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    .normarg_dims(dims, "colMeans")
    if (is(x, "DelayedArray") && x@is_transposed)
        return(.DelayedMatrix_block_rowMeans(t(x), na.rm=na.rm, dims=dims))

    colmeans_list <- colblock_APPLY(x, colMeans, na.rm=na.rm,
                                    if_empty=numeric(0))
    unlist(colmeans_list, recursive=FALSE)
}

setMethod("rowMeans", "DelayedMatrix", .DelayedMatrix_block_rowMeans)
setMethod("colMeans", "DelayedMatrix", .DelayedMatrix_block_colMeans)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###
### We only support multiplication of an ordinary matrix (typically
### small) by a DelayedMatrix object (typically big). Multiplication of 2
### DelayedMatrix objects is not supported.
###

### Write a new HDF5 dataset to disk. Return an HDF5Matrix object that points
### to this new dataset.
.DelayedMatrix_block_mult_by_left_matrix <- function(x, y)
{
    stopifnot(is.matrix(x),
              is(y, "DelayedMatrix") || is.matrix(y),
              ncol(x) == nrow(y))

    out_file <- getHDF5DumpFile()
    out_name <- getHDF5DumpName()

    ans_type <- typeof(match.fun(type(x))(1) * match.fun(type(y))(1))
    h5createDataset2(out_file, out_name, c(nrow(x), ncol(y)),
                     storage.mode=ans_type)
    on.exit(setHDF5DumpName())

    colblock_APPLY(y,
        function(submatrix) x %*% submatrix,
        out_file=out_file,
        out_name=out_name
    )

    ## TODO: Investigate the possiblity to store the dimnames in the HDF5 file
    ## so the HDF5Array() constructor can bring them back. Then we wouldn't
    ## need to explicitely set them on 'ans' like we do below.
    ans <- HDF5Array(out_file, out_name, type=ans_type)
    ans_rownames <- rownames(x)
    ans_colnames <- colnames(y)
    if (!(is.null(ans_rownames) && is.null(ans_colnames)))
        dimnames(ans) <- list(ans_rownames, ans_colnames)
    ans
}

setMethod("%*%", c("DelayedMatrix", "matrix"),
    function(x, y) t(t(y) %*% t(x))
)

setMethod("%*%", c("matrix", "DelayedMatrix"),
    .DelayedMatrix_block_mult_by_left_matrix
)

setMethod("%*%", c("DelayedMatrix", "DelayedMatrix"),
    function(x, y)
        stop(wmsg("multiplication of 2 DelayedMatrix objects is not ",
                  "supported, only multiplication of an ordinary matrix by ",
                  "a DelayedMatrix object at the moment"))
)

