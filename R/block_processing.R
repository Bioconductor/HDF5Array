### =========================================================================
### Internal utilities for block processing of an array
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### Default block size in bytes.
DEFAULT_BLOCK_SIZE <- 4500000L  # 4.5 Mb

### Atomic type sizes in bytes.
.TYPE_SIZES <- c(
    logical=4L,
    integer=4L,
    numeric=8L,
    double=8L,
    complex=16L,
    character=8L,  # just the overhead of a CHARSXP; doesn't account for the
                   # string data itself
    raw=1L
)

get_block_length <- function(type)
{
    type_size <- .TYPE_SIZES[type]
    idx <- which(is.na(type_size))
    if (length(idx) != 0L) {
        unsupported_types <- unique(type[idx])
        in1string <- paste0(unsupported_types, collapse=", ")
        stop("unsupported type(s): ",  in1string)
    }
    block_size <- getOption("HDF5Array.block.size",
                            default=DEFAULT_BLOCK_SIZE)
    as.integer(block_size / type_size)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Walking on the blocks
###
### 3 utility functions to process array-like objects by block.
###

.as_array_or_matrix <- function(x)
{
    if (length(dim(x)) == 2L)
        return(as.matrix(x))
    as.array(x)
}

### An lapply-like function.
block_APPLY <- function(x, APPLY, ..., if_empty=NULL,
                        out_file=NULL, out_name=NULL, block_len=NULL)
{
    APPLY <- match.fun(APPLY)
    if (is.null(block_len))
        block_len <- get_block_length(type(x))
    blocks <- ArrayBlocks(dim(x), block_len)
    nblock <- length(blocks)
    if (nblock == 0L)
        return(if_empty)
    expand_RangeNSBS <- is.array(x) || !is.null(out_file)
    lapply(seq_len(nblock),
        function(i) {
            subscripts <- DelayedArray:::get_array_block_subscripts(blocks, i,
                                                   expand_RangeNSBS)
            subarray <- DelayedArray:::subset_by_subscripts(x, subscripts)
            if (!is.array(subarray))
                subarray <- .as_array_or_matrix(subarray)
            block_ans <- APPLY(subarray, ...)
            if (is.null(out_file))
                return(block_ans)
            h5write2(block_ans, out_file, out_name, subscripts)
        })
}

### A mapply-like function for conformable arrays.
block_MAPPLY <- function(MAPPLY, ..., if_empty=NULL,
                         out_file=NULL, out_name=NULL, block_len=NULL)
{
    MAPPLY <- match.fun(MAPPLY)
    dots <- unname(list(...))
    dims <- sapply(dots, dim)
    x_dim <- dims[ , 1L]
    if (!all(dims == x_dim))
        stop("non-conformable arrays")
    if (is.null(block_len)) {
        types <- unlist(lapply(dots, type))
        block_len <- min(get_block_length(types))
    }
    blocks <- ArrayBlocks(x_dim, block_len)
    nblock <- length(blocks)
    if (nblock == 0L)
        return(if_empty)
    lapply(seq_len(nblock),
        function(i) {
            subscripts <- DelayedArray:::get_array_block_subscripts(blocks, i,
                                                   TRUE)
            subarrays <- lapply(dots,
                function(x) {
                    subarray <- DelayedArray:::subset_by_subscripts(x,
                                                                    subscripts)
                    if (!is.array(subarray))
                        subarray <- .as_array_or_matrix(subarray)
                    subarray
                })
            block_ans <- do.call(MAPPLY, subarrays)
            if (is.null(out_file))
                return(block_ans)
            h5write2(block_ans, out_file, out_name, subscripts)
        })
}

### A Reduce-like function.
block_REDUCE_and_COMBINE <- function(x, REDUCE, COMBINE, init,
                                        BREAKIF=NULL, block_len=NULL)
{
    REDUCE <- match.fun(REDUCE)
    COMBINE <- match.fun(COMBINE)
    if (!is.null(BREAKIF))
        BREAKIF <- match.fun(BREAKIF)
    if (is.null(block_len))
        block_len <- get_block_length(type(x))
    blocks <- ArrayBlocks(dim(x), block_len)
    for (i in seq_along(blocks)) {
        subarray <- DelayedArray:::extract_array_block(x, blocks, i)
        if (!is.array(subarray))
            subarray <- .as_array_or_matrix(subarray)
        reduced <- REDUCE(subarray)
        init <- COMBINE(i, subarray, init, reduced)
        if (!is.null(BREAKIF) && BREAKIF(init))
            break
    }
    init
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Walking on the blocks of columns
###
### 2 convenience wrappers around block_APPLY() and block_REDUCE_and_COMBINE()
### to process a matrix-like object by block of columns.
###

colblock_APPLY <- function(x, APPLY, ..., if_empty=NULL,
                           out_file=NULL, out_name=NULL)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    APPLY <- match.fun(APPLY)
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    block_len <- max(get_block_length(type(x)), x_dim[[1L]])
    block_APPLY(x, APPLY, ..., if_empty=if_empty,
                out_file=out_file, out_name=out_name, block_len=block_len)
}

colblock_REDUCE_and_COMBINE <- function(x, REDUCE, COMBINE, init)
{
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    ## We're going to walk along the columns so need to increase the block
    ## length so each block is made of at least one column.
    block_len <- max(get_block_length(type(x)), x_dim[[1L]])
    block_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init, block_len=block_len)
}

