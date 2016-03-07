### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------
###

SLICE_LENGTH <- 8000000L

setClass("ArraySlicing",
    representation(
        dim="integer",
        N="integer",
        by="integer"
    )
)

### Return an ArraySlicing where each slice is guaranteed to have a
### length <= 'slice_len'.
ArraySlicing <- function(dim, slice_len)
{
    p <- cumprod(dim)
    stopifnot(p[[length(p)]] > slice_len)
    N <- which(p >= slice_len)[[1L]]
    if (p[[N]] == slice_len) {
        N <- N + 1L
        by <- 1L
    } else if (N == 1L) {
        by <- slice_len
    } else {
        by <- slice_len %/% as.integer(p[[N - 1L]])
    }
    new("ArraySlicing", dim=dim, N=N, by=by)
}

.get_ArraySlicing_inner_length <- function(x)
{
    inner_length <- x@dim[[x@N]] %/% x@by
    bottom_slice_size <- x@dim[[x@N]] %% x@by
    if (bottom_slice_size != 0L)
        inner_length <- inner_length + 1L
    inner_length
}

.get_ArraySlicing_outer_length <- function(x)
{
    if (x@N < length(x@dim)) {
        outer_dim <- x@dim[(x@N + 1L):length(x@dim)]
        outer_length <- prod(outer_dim)
    } else {
        outer_length <- 1L
    }
    outer_length
}

### Return the number of slices in 'x'.
setMethod("length", "ArraySlicing",
    function(x)
        .get_ArraySlicing_inner_length(x) * .get_ArraySlicing_outer_length(x)
)

.make_array_slice_index <- function(slicing, i)
{
    #index <- vector("list", length(slicing@dim))
    index <- rep(alist(foo=), length(slicing@dim))
    names(index) <- NULL

    i <- i - 1L
    if (slicing@N < length(slicing@dim)) {
        inner_length <- .get_ArraySlicing_inner_length(slicing)
        i1 <- i %% inner_length
        i2 <- i %/% inner_length
    } else {
        i1 <- i
    }

    k1 <- i1 * slicing@by
    k2 <- k1 + slicing@by
    k1 <- k1 + 1L
    if (k2 > slicing@dim[[slicing@N]])
        k2 <- slicing@dim[[slicing@N]]
    index[[slicing@N]] <- k1:k2

    if (slicing@N < length(slicing@dim)) {
        outer_dim <- slicing@dim[(slicing@N + 1L):length(slicing@dim)]
        subindex <- arrayInd(i2 + 1L, outer_dim)
        index[(slicing@N + 1L):length(slicing@dim)] <- as.list(subindex)
    }
    index
}

get_array_slice <- function(x, slicing, i)
{
    index <- .make_array_slice_index(slicing, i)
    do.call(`[`, c(list(x), index, drop=FALSE))
}

### Should be a no-op. This is a sanity check for the slicing mechanism.
.HDF5Array_slice_and_glue <- function(x, slice_len)
{
    slicing <- ArraySlicing(dim(x), slice_len)
    all_slices <- lapply(seq_along(slicing),
                         function(i) get_array_slice(x, slicing, i))
    ans <- unlist(all_slices, recursive=FALSE)
    dim(ans) <- dim(x)
    ans
}

.HDF5Array_anyNA <- function(x, recursive=FALSE)
{
    if (length(x) <= SLICE_LENGTH)
        return(anyNA(as.array(x)))
    slicing <- ArraySlicing(dim(x), SLICE_LENGTH)
    for (i in seq_along(slicing)) {
        slice <- get_array_slice(x, slicing, i)
        if (anyNA(slice))
            return(TRUE)
    }
    FALSE
}

setMethod("anyNA", "HDF5Array", .HDF5Array_anyNA)

