### =========================================================================
### Manage settings for writing output to HDF5
### -------------------------------------------------------------------------
###


.HDF5_dump_settings_envir <- new.env(hash=TRUE, parent=emptyenv())

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpFile <- function(file=paste0(tempfile(), ".h5"))
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path ",
                  "to the current output HDF5 file, that is, to the HDF5 ",
                  "file where all newly created datasets shall be written"))
    if (file.exists(file)) {
        show_content <- TRUE
    } else {
        h5createFile(file)
        show_content <- FALSE
    }
    ## If h5ls() fails (e.g. file exists but is not HDF5), "file" setting must
    ## remain untouched.
    content <- h5ls(file)
    assign("file", file, envir=.HDF5_dump_settings_envir)
    if (show_content)
        return(content)
    return(invisible(content))
}

getHDF5DumpFile <- function()
    get("file", envir=.HDF5_dump_settings_envir)

### A convenience wrapper.
lsHDF5DumpFile <- function() h5ls(getHDF5DumpFile())

assign("auto_inc_ID", 0L, envir=.HDF5_dump_settings_envir)

setHDF5DumpName <- function(name)
{
    if (missing(name)) {
        suppressWarnings(rm(list="name", envir=.HDF5_dump_settings_envir))
        auto_inc_ID <- get("auto_inc_ID", envir=.HDF5_dump_settings_envir)
        return(assign("auto_inc_ID", auto_inc_ID + 1L,
                      envir=.HDF5_dump_settings_envir))
    }
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name of ",
                  "the next dataset to be written to the current output ",
                  "HDF5 file"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
    assign("name", name, envir=.HDF5_dump_settings_envir)
}

getHDF5DumpName <- function()
{
    name <- try(get("name", envir=.HDF5_dump_settings_envir), silent=TRUE)
    if (is(name, "try-error")) {
        auto_inc_ID <- get("auto_inc_ID", envir=.HDF5_dump_settings_envir)
        name <- sprintf("/HDF5ArrayDataset%05d", auto_inc_ID)
    }
    name
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level non-exported stuff
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
    block_len <- get_block_length(storage.mode)
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

