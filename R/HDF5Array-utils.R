### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------
###


### TODO: Let the user control this via a global option.
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
    inner_length <- x@dim[[x@N]] %/% x@by
    last_inner_block_len <- x@dim[[x@N]] %% x@by
    if (last_inner_block_len != 0L)
        inner_length <- inner_length + 1L
    inner_length
}

.get_ArrayBlocks_outer_length <- function(x)
{
    ndim <- length(x@dim)
    if (x@N < ndim) {
        outer_dim <- x@dim[(x@N + 1L):ndim]
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
    nblock <- length(blocks)
    stopifnot(isSingleInteger(i), i >= 1, i <= nblock)

    ndim <- length(blocks@dim)
    subscript <- rep(alist(foo=), ndim)
    names(subscript) <- NULL

    if (nblock == 1L)
        return(subscript)

    i <- i - 1L
    if (blocks@N < ndim) {
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

    if (blocks@N < ndim) {
        outer_dim <- blocks@dim[(blocks@N + 1L):ndim]
        subindex <- arrayInd(i2 + 1L, outer_dim)
        subscript[(blocks@N + 1L):ndim] <- as.list(subindex)
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
    blocks <- ArrayBlocks(dim(x), MAX_BLOCK_LENGTH)
    for (i in seq_along(blocks)) {
        block <- extract_array_block(x, blocks, i)
        if (anyNA(as.vector(block)))
            return(TRUE)
    }
    FALSE
}

setMethod("anyNA", "HDF5Array", .HDF5Array_anyNA)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summary group generic
###
### max(), min(), range(), sum(), prod(), any(), all().
###

setMethod("Summary", "HDF5Array",
    function(x, ..., na.rm=FALSE)
    {
        if (missing(x)) {
            objects <- list(...)
        } else {
            objects <- list(x, ...)
        }
        ans <- suppressWarnings(callGeneric(NULL))  # init value
        for (x in objects) {
            x <- straight(x)
            blocks <- ArrayBlocks(dim(x), MAX_BLOCK_LENGTH)
            for (i in seq_along(blocks)) {
                block <- extract_array_block(x, blocks, i)
                block_ans <- callGeneric(as.vector(block), na.rm=na.rm)
                ## Early bailout for any() and all().
                if (.Generic == "any") {
                    if (identical(block_ans, TRUE))
                        return(TRUE)
                } else if (.Generic == "all") {
                    if (identical(block_ans, FALSE))
                        return(FALSE)
                }
                ans <- callGeneric(ans, block_ans)
            }
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### Same arguments as base::mean.default().
.mean.HDF5Array <- function(x, trim=0, na.rm=FALSE)
{
    if (!identical(trim, 0))
        stop("\"mean\" method for HDF5Array objects ",
             "does not support the 'trim' argument yet")
    x <- straight(x)
    blocks <- ArrayBlocks(dim(x), MAX_BLOCK_LENGTH)
    sum <- nval <- 0
    for (i in seq_along(blocks)) {
        block <- extract_array_block(x, blocks, i)
        block <- as.vector(block, mode="numeric")
        block_sum <- sum(block, na.rm=na.rm)
        if (is.na(block_sum))
            return(block_sum)
        if (na.rm) {
            block_nval <- sum(!is.na(block))
        } else {
            block_nval <- length(block)
        }
        sum <- sum + block_sum
        nval <- nval + block_nval
    }
    sum / nval
}

### S3/S4 combo for mean.HDF5Array
mean.HDF5Array <- function(x, trim=0, na.rm=FALSE, ...)
    .mean.HDF5Array(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "HDF5Array", .mean.HDF5Array)

