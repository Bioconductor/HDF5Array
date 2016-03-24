### =========================================================================
### Manage settings for writing output to HDF5
### -------------------------------------------------------------------------
###


.HDF5_output_settings_envir <- new.env(hash=TRUE, parent=emptyenv())

### Called by .onLoad() hook (see zzz.R file).
setHDF5OutputFile <- function(file=paste0(tempfile(), ".h5"))
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
    assign("file", file, envir=.HDF5_output_settings_envir)
    if (show_content)
        return(content)
    return(invisible(content))
}

getHDF5OutputFile <- function()
    get("file", envir=.HDF5_output_settings_envir)

### A convenience wrapper.
lsHDF5OutputFile <- function() h5ls(getHDF5OutputFile())

assign("auto_inc_ID", 0L, envir=.HDF5_output_settings_envir)

setHDF5OutputName <- function(name)
{
    if (missing(name)) {
        suppressWarnings(rm(list="name", envir=.HDF5_output_settings_envir))
        auto_inc_ID <- get("auto_inc_ID", envir=.HDF5_output_settings_envir)
        return(assign("auto_inc_ID", auto_inc_ID + 1L,
                      envir=.HDF5_output_settings_envir))
    }
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name of ",
                  "the next dataset to be written to the current output ",
                  "HDF5 file"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
    assign("name", name, envir=.HDF5_output_settings_envir)
}

getHDF5OutputName <- function()
{
    name <- try(get("name", envir=.HDF5_output_settings_envir), silent=TRUE)
    if (is(name, "try-error")) {
        auto_inc_ID <- get("auto_inc_ID", envir=.HDF5_output_settings_envir)
        name <- sprintf("/HDF5ArrayDataset%05d", auto_inc_ID)
    }
    name
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level non-exported stuff
###

### Here is the trade-off: The shorter the chunks, the snappier the "show"
### method feels (it starts to feel sloppy with a chunk length >= 10 millions).
### OTOH small chunks make methods that walk on the entire dataset (e.g.
### range()) slower. Setting the default to 1 million seems a good compromise.
.compute_chunk <- function(dims, chunk_len=1000000L)
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

### A simple wrapper around h5createDataset() that tries to automatically
### choose a reasonable chunk size.
h5createDataset2 <- function(file, dataset, dims, storage.mode="double")
{
    chunk <- .compute_chunk(dims)
    ok <- h5createDataset(file, dataset, dims, storage.mode=storage.mode,
                                               chunk=chunk)
    if (!ok)
        stop(wmsg("failed to create dataset '", dataset, "' ",
                  "in file '", file, "'"), call.=FALSE)
}

