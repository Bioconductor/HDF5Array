### =========================================================================
### Bind DelayedArray objects along their rows or columns
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
    IRanges:::combine_dimnames_along(x@seeds, .get_matrices_dims(x@seeds), 1L)
}

.get_ColBinder_dimnames <- function(x)
{
    IRanges:::combine_dimnames_along(x@seeds, .get_matrices_dims(x@seeds), 2L)
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
### "rbind" and "cbind" methods for DelayedMatrix objects
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
### "rbind" and "cbind" methods for DelayedArray objects
###

.as_DelayedMatrix_objects <- function(objects)
{
    lapply(objects,
        function(object) {
            if (length(dim(object)) != 2L)
                stop(wmsg("cbind() and rbind() are not supported on ",
                          "DelayedArray objects that don't have exactly ",
                          "2 dimensions"))
            as(object, "DelayedMatrix")
        })
}

.DelayedArray_rbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("rbind", objects)
}

.DelayedArray_cbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("cbind", objects)
}

setMethod("rbind", "DelayedArray", .DelayedArray_rbind)
setMethod("cbind", "DelayedArray", .DelayedArray_cbind)

