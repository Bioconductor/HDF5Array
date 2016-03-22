### =========================================================================
### DelayedMatrix objects
### -------------------------------------------------------------------------


### Extending DataTable gives us a few things for free (head(), tail(),
### etc...)
setClass("DelayedMatrix",
    contains=c("DelayedArray", "DataTable"),
    prototype=prototype(
        seeds=list(new("matrix")),
        index=list(integer(0), integer(0)),
        subindex=1:2
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.from_DelayedArray_to_DelayedMatrix <- function(from)
{
    from_dim <- dim(from)
    from_ndim <- length(from_dim)
    if (from_ndim < 2L)
        stop(wmsg(class(from), " object with less than 2 dimensions cannot ",
                  "be coerced to a DelayedMatrix object at the moment"))
    if (from_ndim == 2L)
        return(new2("DelayedMatrix", from))
    idx <- which(from_dim != 1L)
    if (length(idx) > 2L)
        stop(wmsg("Array-like object with more than 2 effective ",
                  "dimensions cannot be coerced to a DelayedMatrix ",
                  "object. ", slicing_tip))
    if (length(idx) == 2L) {
        n1n2 <- idx
    } else if (length(idx) == 1L && idx[[1L]] != 1L) {
        n1n2 <- c(1L, idx[[1L]])
    } else {
        n1n2 <- 1:2
    }
    new2("DelayedMatrix", from, subindex=from@subindex[n1n2])
}

setAs("DelayedArray", "DelayedMatrix", .from_DelayedArray_to_DelayedMatrix)

### array -> DelayedMatrix

.from_array_to_DelayedMatrix <- function(from)
{
    as(as(from, "DelayedArray"), "DelayedMatrix")
}

setAs("array", "DelayedMatrix", .from_array_to_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Specialized constructor
###

### 'x' must be an array-like object with 3 dimensions.
### 'MARGIN' is the dimension to drop.
make_DelayedMatrix_from_3D_DelayedArray <- function(x, MARGIN)
{
    if (!is(x, "DelayedArray"))
        x <- as(x, "DelayedArray")

    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (x_ndim != 3L)
        stop("'x' must have 3 dimensions")
    if (!isSingleNumber(MARGIN))
        stop("'MARGIN' must be a single integer")
    if (!is.integer(MARGIN))
        MARGIN <- as.integer(MARGIN)
    if (MARGIN < 1L || MARGIN > x_ndim)
        stop("'MARGIN' must be >= 1 and <= length(dim(x))")
    if (x_dim[[MARGIN]] != 1L)
        stop("'dim(x)[[MARGIN]]' must be 1")

    if (x@is_transposed)
        MARGIN <- x_ndim + 1L - MARGIN
    n1n2 <- seq_len(x_ndim)[-MARGIN]
    new2("DelayedMatrix", x, subindex=x@subindex[n1n2])
}

