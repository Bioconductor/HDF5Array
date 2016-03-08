### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------
###


MAX_BLOCK_LENGTH <- 10000000L


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
ArrayBlocks <- function(dim, max_block_len)
{
    p <- cumprod(dim)
    x_len <- p[[length(p)]]
    stopifnot(x_len > max_block_len)
    N <- which(p >= max_block_len)[[1L]]
    if (p[[N]] == max_block_len) {
        N <- N + 1L
        by <- 1L
    } else if (N == 1L) {
        by <- max_block_len
    } else {
        by <- max_block_len %/% as.integer(p[[N - 1L]])
    }
    new("ArrayBlocks", dim=dim, N=N, by=by)
}

.get_ArrayBlocks_inner_length <- function(x)
{
    inner_length <- x@dim[[x@N]] %/% x@by
    last_inner_block_len <- x@dim[[x@N]] %% x@by
    if (last_inner_block_len != 0L)
        inner_length <- inner_length + 1L
    inner_length
}

.get_ArrayBlocks_outer_length <- function(x)
{
    if (x@N < length(x@dim)) {
        outer_dim <- x@dim[(x@N + 1L):length(x@dim)]
        outer_length <- prod(outer_dim)
    } else {
        outer_length <- 1L
    }
    outer_length
}

### Return the number of blocks in 'x'.
setMethod("length", "ArrayBlocks",
    function(x)
        .get_ArrayBlocks_inner_length(x) * .get_ArrayBlocks_outer_length(x)
)

### Return a multidimensional subscript as a list with 1 subscript per
### dimension in the original array.
.get_array_block_subscript <- function(blocks, i)
{
    subscript <- rep(alist(foo=), length(blocks@dim))
    names(subscript) <- NULL

    i <- i - 1L
    if (blocks@N < length(blocks@dim)) {
        inner_length <- .get_ArrayBlocks_inner_length(blocks)
        i1 <- i %% inner_length
        i2 <- i %/% inner_length
    } else {
        i1 <- i
    }

    k1 <- i1 * blocks@by
    k2 <- k1 + blocks@by
    k1 <- k1 + 1L
    if (k2 > blocks@dim[[blocks@N]])
        k2 <- blocks@dim[[blocks@N]]
    subscript[[blocks@N]] <- k1:k2

    if (blocks@N < length(blocks@dim)) {
        outer_dim <- blocks@dim[(blocks@N + 1L):length(blocks@dim)]
        subindex <- arrayInd(i2 + 1L, outer_dim)
        subscript[(blocks@N + 1L):length(blocks@dim)] <- as.list(subindex)
    }
    subscript
}

extract_array_block <- function(x, blocks, i)
{
    subscript <- .get_array_block_subscript(blocks, i)
    do.call(`[`, c(list(x), subscript, drop=FALSE))
}

### NOT exported. Used in unit tests.
split_array_in_blocks <- function(x, max_block_len)
{
    blocks <- ArrayBlocks(dim(x), max_block_len)
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
### anyNA()
###

.HDF5Array_anyNA <- function(x, recursive=FALSE)
{
    x <- straight(x)
    if (length(x) <= MAX_BLOCK_LENGTH)
        return(anyNA(as.vector(x)))
    blocks <- ArrayBlocks(dim(x), MAX_BLOCK_LENGTH)
    for (i in seq_along(blocks)) {
        block <- extract_array_block(x, blocks, i)
        if (anyNA(as.vector(block)))
            return(TRUE)
    }
    FALSE
}

setMethod("anyNA", "HDF5Array", .HDF5Array_anyNA)

