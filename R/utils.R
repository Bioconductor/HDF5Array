### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Simple wrappers around rhdf5::h5read() and rhdf5::h5write()
###

.set_missing_subscripts_to_NULL <- function(subscripts)
{
    stopifnot(is.list(subscripts))
    ans <- vector("list", length(subscripts))
    not_missing_idx <- which(!vapply(subscripts, is.name, logical(1)))
    ans[not_missing_idx] <- subscripts[not_missing_idx]
    ans
}

.make_index_from_subscripts <- function(subscripts)
{
    stopifnot(is.list(subscripts))
    subscripts <- DelayedArray:::expand_RangeNSBS_subscripts(subscripts)
    .set_missing_subscripts_to_NULL(subscripts)
}

h5read2 <- function(file, name, index=NULL)
{
    if (!is.null(index))
        index <- .make_index_from_subscripts(index)
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    suppressWarnings(h5read(file, name, index=index))
}

h5write2 <- function(obj, file, name, index=NULL)
{
    if (!is.null(index))
        index <- .make_index_from_subscripts(index)
    h5write(obj, file, name, index=index)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A simple wrapper around rhdf5::h5createDataset() that automatically
### chooses the chunk geometry
###

### Here is the trade-off: The shorter the chunks, the snappier the "show"
### method feels (on my laptop, it starts to feel sloppy with a chunk
### length > 10 millions). OTOH small chunks tend to slow down methods that
### do block processing (e.g. sum(), range(), etc...). Setting the default
### to 1 million seems a good compromise.
.chunk_as_hypercube <- function(dims, chunk_len=1000000L)
{
    if (prod(dims) <= chunk_len)
        return(dims)
    ndim <- length(dims)

    ## The perfect chunk is the hypercube.
    chunk <- as.integer(round(rep.int(chunk_len ^ (1 / ndim), ndim)))

    ## But it could have dimensions that are greater than 'dims'. In that case
    ## we need to reshape it.
    while (any(chunk > dims)) {
        chunk <- pmin(chunk, dims)
        r <- chunk_len / prod(chunk)  # > 1
        extend_along <- which(chunk < dims)
        extend_factor <- r ^ (1 / length(extend_along))
        chunk[extend_along] <- as.integer(chunk[extend_along] * extend_factor)
    }
    chunk
}

.chunk_as_subblock <- function(dims, storage.mode="double", ratio=75L)
{
    block_len <- DelayedArray:::get_block_length(storage.mode)
    chunk_len <- as.integer(ceiling(block_len / ratio))
    ## 'block_len' must be a multiple of 'chunk_len'.
    stopifnot(block_len %% chunk_len == 0L)
    chunks <- ArrayBlocks(dims, chunk_len)
    chunk <- chunks@dim
    ndim <- length(chunk)
    if (chunks@N > ndim)
        return(chunk)
    chunk[[chunks@N]] <- chunks@by
    if (chunks@N == ndim)
        return(chunk)
    chunk[(chunks@N+1L):ndim] <- 1L
    chunk
}

### A simple wrapper around h5createDataset() that automatically chooses the
### chunk geometry.
h5createDataset2 <- function(file, dataset, dims, storage.mode="double")
{
    if (storage.mode == "character") {
        size <- max(nchar(dataset, type="width"))
    } else {
        size <- NULL
    }
    #chunk <- .chunk_as_hypercube(dims)
    chunk <- .chunk_as_subblock(dims, storage.mode)
    ok <- h5createDataset(file, dataset, dims, storage.mode=storage.mode,
                                               size=size,
                                               chunk=chunk)
    if (!ok)
        stop(wmsg("failed to create dataset '", dataset, "' ",
                  "in file '", file, "'"), call.=FALSE)
}

