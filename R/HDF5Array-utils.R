### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Arith" and "Math" group generics
###
### Arith members: "+", "-", "*", "/", "^", "%%", "%/%"
### Math members: abs, sign, sqrt, ceiling, floor, and many more...
###

.HDF5Array_delayed_Arith <- function(.Generic, e1, e2)
{
    if (isSingleNumberOrNA(e2)) {
        ans <- register_delayed_op(e1, .Generic, Rargs=list(e2))
    } else if (isSingleNumberOrNA(e1)) {
        ans <- register_delayed_op(e2, .Generic, Largs=list(e1))
    } else {
        stop(wmsg("arithmetic operations on HDF5Array objects are supported ",
                  "only when one operand is a single number for now"))
    }
    ans
}

setMethod("Arith", c("HDF5Array", "numeric"),
    function(e1, e2) .HDF5Array_delayed_Arith(.Generic, e1, e2)
)

setMethod("Arith", c("numeric", "HDF5Array"),
    function(e1, e2) .HDF5Array_delayed_Arith(.Generic, e1, e2)
)

setMethod("Arith", c("HDF5Array", "HDF5Array"),
    function(e1, e2) .HDF5Array_delayed_Arith(.Generic, e1, e2)
)

setMethod("Math", "HDF5Array",
    function(x) register_delayed_op(x, .Generic)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Compare" group generic
###
### Members: ==, !=, <=, >=, <, >
###

.HDF5Array_block_Compare <- function(.Generic, e1, e2)
{
    if (!identical(dim(e1), dim(e2)))
        stop("comparing 2 ", class(e1), " objects with different ",
             "dimensions is not supported yet")
    GENERIC <- match.fun(.Generic)
    res_list <- block_MAPPLY(.Generic, e1, e2)
    ans <- unlist(res_list, recursive=FALSE, use.names=FALSE)
    dim(ans) <- dim(e1)
    ans
}

setMethod("Compare", "HDF5Array",
    function(e1, e2)
        .HDF5Array_block_Compare(.Generic, e1, e2)
)

