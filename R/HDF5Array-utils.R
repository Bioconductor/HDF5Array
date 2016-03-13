### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------


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
    init <- FALSE
    APPLY <- anyNA
    REDUCE <- `||`
    BREAKIF <- identity

    x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
    block_APPLY_REDUCE(x, init, APPLY, REDUCE, BREAKIF)
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
    init <- NULL
    APPLY <- function(subarray) {
        ## We get a warning if 'subarray' is empty (which can't happen, blocks
        ## can't be empty) or if 'na.rm' is TRUE and 'subarray' contains only
        ## NA's or NaN's.
        val <- tryCatch(GENERIC(subarray, na.rm=na.rm), warning=identity)
        if (is(val, "warning"))
            return(NULL)
        val
    }
    REDUCE <- function(init, val) {
        if (is.null(init) && is.null(val))
            return(NULL)
        GENERIC(init, val)
    }
    BREAKIF <- function(init) {
        if (is.null(init))
            return(FALSE)
        switch(.Generic,
            max=         is.na(init) || init == Inf,
            min=         is.na(init) || init == -Inf,
            range=       is.na(init[[1L]]) || all(init == c(-Inf, Inf)),
            sum=, prod=  is.na(init),
            any=         identical(init, TRUE),
            all=         identical(init, FALSE),
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
        init <- block_APPLY_REDUCE(x, init, APPLY, REDUCE, BREAKIF)
    }
    if (is.null(init))
        init <- GENERIC()
    init
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

    init <- numeric(2)  # sum and nval
    APPLY <- function(subarray) {
        tmp <- as.vector(subarray, mode="numeric")
        subarray_sum <- sum(tmp, na.rm=na.rm)
        subarray_nval <- length(tmp)
        if (na.rm)
            subarray_nval <- subarray_nval - sum(is.na(tmp))
        c(subarray_sum, subarray_nval)
    }
    REDUCE <- `+`
    BREAKIF <- function(init) is.na(init[[1L]])

    x <- .straighten(x, untranspose=TRUE)
    init <- block_APPLY_REDUCE(x, init, APPLY, REDUCE, BREAKIF)
    init[[1L]] / init[[2L]]
}

### S3/S4 combo for mean.HDF5Array
mean.HDF5Array <- function(x, trim=0, na.rm=FALSE, ...)
    .HDF5Array_block_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "HDF5Array", .HDF5Array_block_mean)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Compare" group generic
###
### Members: ==, !=, <=, >=, <, >
###

#setMethod("Compare", "HDF5Array",
#    function(e1, e2)
#    {
#        e1_len <- length(e1)
#        e2_len <- length(e2)
#    }
#)

