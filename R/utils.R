### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5dim()
###

h5dim <- function(file, name)
{
    if (substr(name, 1L, 1L) != "/")
        name <- paste0("/", name)
    group <- gsub("(.*/)[^/]*$", "\\1", name)
    name <- gsub(".*/([^/]*)$", "\\1", name)
    f <- H5Fopen(file, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(f))
    g <- H5Gopen(f, group)
    on.exit(H5Gclose(g), add=TRUE)
    d <- H5Dopen(g, name)
    on.exit(H5Dclose(d), add=TRUE)
    H5Sget_simple_extent_dims(H5Dget_space(d))$size
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Simple wrappers around rhdf5::h5read() and rhdf5::h5write()
###

h5read2 <- function(file, name, index=NULL)
{
    if (!is.null(index))
        index <- DelayedArray:::expand_Nindex_RangeNSBS(index)
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    suppressWarnings(h5read(file, name, index=index))
}

h5write2 <- function(obj, file, name, index=NULL)
{
    if (!is.null(index))
        index <- DelayedArray:::expand_Nindex_RangeNSBS(index)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Detect and trim trailing slahes in a character vector
###

has_trailing_slash <- function(x)
{
    stopifnot(is.character(x))
    #nc <- nchar(x)
    #substr(x, start=nc, stop=nc) == "/"
    grepl("/$", x)  # seems slightly faster than the above
}

trim_trailing_slashes <- function(x)
{
    sub("/*$", "", x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A simple mechanism to lock/unlock a file so processes can get temporary
### exclusive access to it
###

.locked_path <- function(filepath)
{
    if (!isSingleString(filepath) || filepath == "")
        stop("'filepath' must be a single non-empty string")
    paste0(filepath, "-locked")
}

.safe_file_rename <- function(from, to)
{
    !file.exists(to) && suppressWarnings(file.rename(from, to))
}

lock_file <- function(filepath)
{
    locked_path <- .locked_path(filepath)
    ## Must wait if the file is already locked.
    while (TRUE) {
        if (.safe_file_rename(filepath, locked_path))
            break
        Sys.sleep(0.01)
    }
    locked_path
}

unlock_file <- function(filepath)
{
    locked_path <- .locked_path(filepath)
    if (!.safe_file_rename(locked_path, filepath))
        stop("failed to unlock '", filepath, "' file")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A global counter that is safe to use in the context of parallelized
### execution
###

.read_counter <- function(filepath)
{
    counter <- readLines(filepath)
    stopifnot(length(counter) == 1L)
    counter <- suppressWarnings(as.integer(counter))
    if (is.na(counter))
        stop("file '", filepath, "' does not contain a counter")
    counter
}

### Will overwrite an existing file.
.write_counter <- function(counter, filepath)
{
    writeLines(as.character(counter), filepath)
    counter
}

### NOT safe to use in the context of parallel execution!
init_global_counter <- function(filepath, counter=1L)
{
    if (!isSingleString(filepath) || filepath == "")
        stop("'filepath' must be a single non-empty string")
    if (file.exists(filepath))
        stop("file '", filepath, "' already exists")
    if (!isSingleNumber(counter))
        stop("'counter' must be a single number")
    if (!is.integer(counter))
        counter <- as.integer(counter)
    .write_counter(counter, filepath)
}

### Use a lock mechanism to prevent several processes from trying to increment
### the counter simultaneously. So is safe to use in the context of parallel
### execution e.g.
###
###   library(BiocParallel)
###   filepath <- tempfile()
###   init_global_counter(filepath)
###   bplapply(1:10, function(i) get_global_counter(filepath, increment=TRUE))
###
get_global_counter <- function(filepath, increment=FALSE)
{
    if (!isTRUEorFALSE(increment))
        stop("'increment' must be TRUE or FALSE")
    locked_path <- lock_file(filepath)
    on.exit(unlock_file(filepath))
    counter <- .read_counter(locked_path)
    if (increment)
        .write_counter(counter + 1L, locked_path)
    counter
}

