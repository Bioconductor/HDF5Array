### =========================================================================
### Internal utilities for block processing of an array
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


DEFAULT_BLOCK_SIZE <- 40000000L

get_block_length <- function(type)
{
    type_size <- switch(type,
        logical=, integer=4L,
        numeric=, double=8L,
        complex=16L,
        character=8L,  # just the overhead of a CHARSXP, doesn't count for
                       # the string data itself
        raw=1L,
        stop("type ", type, " is not supported")
    )
    block_size <- getOption("HDF5Array.block.size",
                            default=.DEFAULT_BLOCK_SIZE)
    as.integer(block_size / type_size)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() getter
###
### For internal use only.
###

setGeneric("type", function(x) standardGeneric("type"))

setMethod("type", "array", function(x) typeof(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayBlocks objects
###

setClass("ArrayBlocks",
    representation(
        dim="integer",
        N="integer",
        by="integer"
    )
)

### Return an ArrayBlocks object i.e. a collection of subarrays of the
### original array with the following properties:
###   (a) The collection of blocks is a partitioning of the original array
###       i.e. the blocks fully cover it and don't overlap each other.
###   (b) Each block is made of adjacent elements in the original array.
###   (c) Each block has a length (i.e. nb of elements) <= 'max_block_len'.
.ArrayBlocks <- function(dim, max_block_len)
{
    ndim <- length(dim)
    p <- cumprod(dim)
    x_len <- p[[ndim]]
    if (max_block_len > x_len) {
        N <- ndim + 1L
        by <- 1L
    } else {
        N <- which(p >= max_block_len)[[1L]]
        if (p[[N]] == max_block_len) {
            N <- N + 1L
            by <- 1L
        } else if (N == 1L) {
            by <- max_block_len
        } else {
            by <- max_block_len %/% as.integer(p[[N - 1L]])
        }
    }
    new("ArrayBlocks", dim=dim, N=N, by=by)
}

.get_ArrayBlocks_inner_length <- function(x)
{
    ndim <- length(x@dim)
    if (x@N > ndim)
        return(1L)
    inner_len <- x@dim[[x@N]] %/% x@by
    last_inner_block_len <- x@dim[[x@N]] %% x@by
    if (last_inner_block_len != 0L)
        inner_len <- inner_len + 1L
    inner_len
}

.get_ArrayBlocks_outer_length <- function(x)
{
    ndim <- length(x@dim)
    if (x@N < ndim) {
        outer_dim <- x@dim[(x@N + 1L):ndim]
        outer_len <- prod(outer_dim)
    } else {
        outer_len <- 1L
    }
    outer_len
}

### Return the number of blocks in 'x'.
setMethod("length", "ArrayBlocks",
    function(x)
        .get_ArrayBlocks_inner_length(x) * .get_ArrayBlocks_outer_length(x)
)

### Return a multidimensional subscript as a list with 1 subscript per
### dimension in the original array.
.get_array_block_subscript <- function(blocks, i, expand.RangeNSBS=FALSE)
{
    nblock <- length(blocks)
    stopifnot(isSingleInteger(i), i >= 1, i <= nblock)

    ndim <- length(blocks@dim)
    subscript <- rep(alist(foo=), ndim)
    names(subscript) <- NULL

    if (nblock == 1L)
        return(subscript)

    i <- i - 1L
    if (blocks@N < ndim) {
        inner_len <- .get_ArrayBlocks_inner_length(blocks)
        i1 <- i %% inner_len
        i2 <- i %/% inner_len
    } else {
        i1 <- i
    }

    k1 <- i1 * blocks@by
    k2 <- k1 + blocks@by
    k1 <- k1 + 1L
    upper_bound <- blocks@dim[[blocks@N]]
    if (k2 > upper_bound)
        k2 <- upper_bound
    if (expand.RangeNSBS) {
        subscript_N <- k1:k2  # same as doing as.integer() on the RangeNSBS
                              # object below
    } else {
        subscript_N <- new2("RangeNSBS", subscript=c(k1, k2),
                                         upper_bound=upper_bound,
                                         check=FALSE)
    }
    subscript[[blocks@N]] <- subscript_N

    if (blocks@N < ndim) {
        outer_dim <- blocks@dim[(blocks@N + 1L):ndim]
        subindex <- arrayInd(i2 + 1L, outer_dim)
        subscript[(blocks@N + 1L):ndim] <- as.list(subindex)
    }
    subscript
}

extract_array_block <- function(x, blocks, i)
{
    subscript <- .get_array_block_subscript(blocks, i,
                                            expand.RangeNSBS=is.array(x))
    do.call(`[`, c(list(x), subscript, drop=FALSE))
}

### NOT exported. Used in unit tests.
split_array_in_blocks <- function(x, max_block_len)
{
    blocks <- .ArrayBlocks(dim(x), max_block_len)
    lapply(seq_along(blocks),
           function(i) extract_array_block(x, blocks, i))
}

### NOT exported. Used in unit tests.
### Rebuild the original array from the subarrays obtained by
### split_array_in_blocks() as an *ordinary* array.
### So if 'x' is an ordinary array, then:
###
###   unsplit_array_from_blocks(split_array_in_blocks(x, max_block_len), x)
###
### should be a no-op for any 'max_block_len' < 'length(x)'.
unsplit_array_from_blocks <- function(subarrays, x)
{
    ans <- combine_array_objects(subarrays)
    dim(ans) <- dim(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Block processing
###

### The core block-processing engine.
block_APPLY_REDUCE <- function(x, APPLY, REDUCE, reduced,
                               BREAKIF=NULL, block_len=NULL)
{
    APPLY <- match.fun(APPLY)
    REDUCE <- match.fun(REDUCE)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    if (is.null(block_len))
        block_len <- get_block_length(type(x))
    blocks <- .ArrayBlocks(dim(x), block_len)
    for (i in seq_along(blocks)) {
        subarray <- extract_array_block(x, blocks, i)
        if (!is.array(subarray))
            subarray <- as.array(subarray)
        val <- APPLY(subarray)
        reduced <- REDUCE(reduced, val)
        if (!is.null(BREAKIF) && BREAKIF(reduced))
            break
    }
    reduced
}

### Processing a matrix-like object by blocks of columns.
colblock_APPLY_REDUCE <- function(x, APPLY, REDUCE, reduced)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    block_len <- max(get_block_length(type(x)), x_dim[[1L]])
    block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, block_len=block_len)
}


colblock_APPLY <- function(x, APPLY, ...)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    APPLY <- match.fun(APPLY)
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    block_len <- max(get_block_length(type(x)), x_dim[[1L]])
    blocks <- .ArrayBlocks(x_dim, block_len)
    lapply(seq_along(blocks),
        function(i) {
            submatrix <- extract_array_block(x, blocks, i)
            if (!is.matrix(submatrix))
                submatrix <- as.matrix(submatrix)
            APPLY(submatrix, ...)
        })
}

