### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Math" group generic
###
### Math members: abs, sign, sqrt, ceiling, floor, and many more...
###

setMethod("Math", "HDF5Array",
    function(x) register_delayed_op(x, .Generic)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Ops" group generics
###
### Arith members: "+", "-", "*", "/", "^", "%%", "%/%"
### Compare members: ==, !=, <=, >=, <, >
### Logic members: &, |
###

.HDF5Array_delayed_Ops_with_right_vector <- function(.Generic, e1, e2)
{
    e1_class <- class(e1)
    e2_class <- class(e2)
    if (!is.vector(e2))
        e2 <- as.vector(e2)
    if (!is.atomic(e2))
        stop(wmsg("`", .Generic, "` between ", e1_class, " and ",
                  e2_class, " objects is not supported"))
    e2_len <- length(e2)
    if (e2_len == 1L)
        return(register_delayed_op(e1, .Generic, Rargs=list(e2)))
    e1_len <- length(e1)
    if (e2_len > e1_len)
        stop(wmsg("right object is longer than left object"))
    e1_nrow <- nrow(e1)
    if (e1_nrow != 0L) {
        if (e2_len == 0L || e1_nrow %% e2_len != 0L)
            stop(wmsg("length of right object is not a divisor ",
                      "of number of rows of left object"))
        e2 <- rep(e2, length.out=nrow(e1))
    }
    register_delayed_op(e1, .Generic, Rargs=list(e2),
                                      recycle_along_last_dim=e1@is_transposed)
}

.HDF5Array_delayed_Ops_with_left_vector <- function(.Generic, e1, e2)
{
    e1_class <- class(e1)
    e2_class <- class(e2)
    if (!is.vector(e1))
        e1 <- as.vector(e1)
    if (!is.atomic(e1))
        stop(wmsg("`", .Generic, "` between ", e1_class, " and ",
                  e2_class, " objects is not supported"))
    e1_len <- length(e1)
    if (e1_len == 1L)
        return(register_delayed_op(e2, .Generic, Largs=list(e1)))
    e2_len <- length(e2)
    if (e1_len > e2_len)
        stop(wmsg("left object is longer than right object"))
    e2_nrow <- nrow(e2)
    if (e2_nrow != 0L) {
        if (e1_len == 0L || e2_nrow %% e1_len != 0L)
            stop(wmsg("length of left object is not a divisor ",
                      "of number of rows of right object"))
        e1 <- rep(e1, length.out=nrow(e2))
    }
    register_delayed_op(e2, .Generic, Largs=list(e1),
                                      recycle_along_last_dim=e2@is_transposed)
}

.HDF5Array_block_Ops <- function(.Generic, e1, e2)
{
    if (!identical(dim(e1), dim(e2)))
        stop("non-conformable arrays")
    GENERIC <- match.fun(.Generic)
    res_list <- block_MAPPLY(.Generic, e1, e2)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    dim(ans) <- dim(e1)
    ans
}

setMethod("Ops", c("HDF5Array", "vector"),
    function(e1, e2)
        .HDF5Array_delayed_Ops_with_right_vector(.Generic, e1, e2)
)

setMethod("Ops", c("vector", "HDF5Array"),
    function(e1, e2)
        .HDF5Array_delayed_Ops_with_left_vector(.Generic, e1, e2)
)

setMethod("Ops", c("HDF5Array", "HDF5Array"),
    function(e1, e2)
    {
        e1_dim <- dim(e1)
        e2_dim <- dim(e2)
        if (identical(e1_dim, e2_dim))
            return(.HDF5Array_block_Ops(.Generic, e1, e2))
        ## Effective dimensions.
        effdim_idx1 <- which(e1_dim != 1L)
        effdim_idx2 <- which(e2_dim != 1L)
        if ((length(effdim_idx1) == 1L) == (length(effdim_idx2) == 1L))
            stop("non-conformable arrays")
        if (length(effdim_idx1) == 1L) {
            .HDF5Array_delayed_Ops_with_left_vector(.Generic, e1, e2)
        } else {
            .HDF5Array_delayed_Ops_with_right_vector(.Generic, e1, e2)
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A low-level utility for putting HDF5Array object in a "straight" form
###
### Untranspose the HDF5Array object and put its rows and columns in their
### "native" order. The goal is to put the matrix elements in their "native"
### order (i.e. in the same order as on disk) so as.vector() is faster on the
### resulting object is faster than on the original object.
###

.straighten_index <- function(i)
{
    i_len <- length(i)
    if (i_len == 0L)
        return(i)
    i_max <- max(i)
    ## Threshold is a rough estimate obtained empirically.
    ## TODO: Refine this.
    if (i_max <= 2L * i_len * log(i_len))
        which(as.logical(tabulate(i, nbins=i_max)))
    else
        sort(unique(i))
}

.straighten <- function(x, untranspose=FALSE, straighten.index=FALSE)
{
    if (is.array(x))
        return(x)
    if (untranspose)
        x@is_transposed <- FALSE
    if (!straighten.index)
        return(x)
    x_index <- index(x)
    index_was_touched <- FALSE
    for (n in seq_along(x_index)) {
        if (isStrictlySorted(x_index[[n]]))
            next
        x_index[[n]] <- .straighten_index(x_index[[n]])
        index_was_touched <- TRUE
    }
    if (index_was_touched)
        index(x) <- x_index
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

.HDF5Array_block_anyNA <- function(x, recursive=FALSE)
{
    APPLY <- anyNA
    REDUCE <- `||`
    reduced <- FALSE
    BREAKIF <- identity

    x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
    block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, BREAKIF)
}

setMethod("anyNA", "HDF5Array", .HDF5Array_block_anyNA)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Summary" group generic
###
### Members: max, min, range, sum, prod, any, all
###

.HDF5Array_block_Summary <- function(.Generic, x, ..., na.rm=FALSE)
{
    GENERIC <- match.fun(.Generic)
    APPLY <- function(subarray) {
        ## We get a warning if 'subarray' is empty (which can't happen, blocks
        ## can't be empty) or if 'na.rm' is TRUE and 'subarray' contains only
        ## NA's or NaN's.
        val <- tryCatch(GENERIC(subarray, na.rm=na.rm), warning=identity)
        if (is(val, "warning"))
            return(NULL)
        val
    }
    REDUCE <- function(reduced, val) {
        if (is.null(reduced) && is.null(val))
            return(NULL)
        GENERIC(reduced, val)
    }
    reduced <- NULL
    BREAKIF <- function(reduced) {
        if (is.null(reduced))
            return(FALSE)
        switch(.Generic,
            max=         is.na(reduced) || reduced == Inf,
            min=         is.na(reduced) || reduced == -Inf,
            range=       is.na(reduced[[1L]]) || all(reduced == c(-Inf, Inf)),
            sum=, prod=  is.na(reduced),
            any=         identical(reduced, TRUE),
            all=         identical(reduced, FALSE),
            FALSE)  # fallback (actually not needed)
    }

    if (missing(x)) {
        objects <- list(...)
    } else {
        objects <- list(x, ...)
    }
    for (x in objects) {
        if (.Generic %in% c("sum", "prod")) {
            x <- .straighten(x, untranspose=TRUE)
        } else {
            x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
        }
        reduced <- block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, BREAKIF)
    }
    if (is.null(reduced))
        reduced <- GENERIC()
    reduced
}

setMethod("Summary", "HDF5Array",
    function(x, ..., na.rm=FALSE)
        .HDF5Array_block_Summary(.Generic, x, ..., na.rm=na.rm)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### Same arguments as base::mean.default().
.HDF5Array_block_mean <- function(x, trim=0, na.rm=FALSE)
{
    if (!identical(trim, 0))
        stop("\"mean\" method for HDF5Array objects ",
             "does not support the 'trim' argument yet")

    APPLY <- function(subarray) {
        tmp <- as.vector(subarray, mode="numeric")
        subarray_sum <- sum(tmp, na.rm=na.rm)
        subarray_nval <- length(tmp)
        if (na.rm)
            subarray_nval <- subarray_nval - sum(is.na(tmp))
        c(subarray_sum, subarray_nval)
    }
    REDUCE <- `+`
    reduced <- numeric(2)  # sum and nval
    BREAKIF <- function(reduced) is.na(reduced[[1L]])

    x <- .straighten(x, untranspose=TRUE)
    reduced <- block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, BREAKIF)
    reduced[[1L]] / reduced[[2L]]
}

### S3/S4 combo for mean.HDF5Array
mean.HDF5Array <- function(x, trim=0, na.rm=FALSE, ...)
    .HDF5Array_block_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "HDF5Array", .HDF5Array_block_mean)

