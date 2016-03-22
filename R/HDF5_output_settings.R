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
                  "to the HDF5 file where output will be written"))
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
                  "the dataset in the HDF5 file to which output will be ",
                  "written"))
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

