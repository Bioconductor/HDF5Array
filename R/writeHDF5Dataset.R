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
### HDF5DatasetDump objects
###

setClass("HDF5DatasetDump",
    contains="OnDiskArrayDump",
    representation(
        file="character",  # Single string.
        name="character",  # Dataset name.
        dim="integer",
        dimnames="list",
        type="character"   # Single string.
    )
)

setMethod("dimnames", "HDF5DatasetDump",
    function(x)
    {
        ans <- x@dimnames
        if (all(S4Vectors:::sapply_isNULL(ans)))
            return(NULL)
        ans
    }
)

### HDF5DatasetDump object created with an earlier call to HDF5DatasetDump()
### should be closed before calling HDF5DatasetDump() again.
### FIXME: Investigate the possiblity to write the dimnames to the HDF5 file.
HDF5DatasetDump <- function(dim, dimnames=NULL, type="double")
{
    file <- getHDF5DumpFile()
    name <- getHDF5DumpName()
    if (is(try(h5createDataset2(file, name, dim, type)), "try-error"))
        stop(wmsg("Failed to create a new HDF5DatasetDump object. Make sure ",
                  "to close() the previously created HDF5DatasetDump object ",
                  "first. Alternatively call setHDF5DumpName() before trying ",
                  "to call HDF5DatasetDump() again."))
    if (is.null(dimnames)) {
        dimnames <- vector("list", length(dim))
    } else {
        ## TODO: Write the dimnames to the HDF5 file.
    }
    new2("HDF5DatasetDump", file=file, name=name,
                            dim=dim, dimnames=dimnames, type=type)
}

setMethod("write_to_dump", c("array", "HDF5DatasetDump"),
    function(x, dump, subscripts=NULL)
        h5write2(x, dump@file, dump@name, index=subscripts)
)

setMethod("close", "HDF5DatasetDump",
    function(con, ...) setHDF5DumpName()
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### FIXME: This needs to propagate the dimnames. Unfortunately this is not
### possible at the moment. See FIXME right before definition of
### HDF5DatasetDump() above in this file and right before definition of
### HDF5Dataset() in HDF5Array-class.R about this.
.from_HDF5DatasetDump_to_HDF5Dataset <- function(from)
{
    HDF5Dataset(from@file, from@name, type=from@type)
}

setAs("HDF5DatasetDump", "HDF5Dataset", .from_HDF5DatasetDump_to_HDF5Dataset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeHDF5Dataset()
###

### Return an invisible HDF5Dataset object pointing to the newly written HDF5
### dataset on disk.
### FIXME: This needs to propagate the dimnames. Unfortunately this is not
### possible at the moment. See various FIXMEs above in this file about this.
writeHDF5Dataset <- function(x, file, name)
{
    old_dump_file <- getHDF5DumpFile()
    old_dump_name <- getHDF5DumpName()
    setHDF5DumpFile(file)
    on.exit({setHDF5DumpFile(old_dump_file); setHDF5DumpName(old_dump_name)})
    setHDF5DumpName(name)
    dump <- HDF5DatasetDump(dim(x), dimnames(x), type(x))
    write_to_dump(x, dump)
    invisible(as(dump, "HDF5Dataset"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercing an array-like object to HDF5Dataset, HDF5Array, or HDF5Matrix
### dumps it to disk.
###

.dump_as_HDF5Dataset <- function(from)
{
    dump <- HDF5DatasetDump(dim(from), dimnames(from), type(from))
    on.exit(close(dump))
    write_to_dump(from, dump)
    as(dump, "HDF5Dataset")
}

setAs("ANY", "HDF5Dataset", .dump_as_HDF5Dataset)

.dump_as_HDF5Array <- function(from)
{
    ans <- as(as(from, "HDF5Dataset"), "HDF5Array")
    ## Temporarily needed because coercion from HDF5DatasetDump to HDF5Dataset
    ## doesn't propagate the dimnames at the moment. See FIXME above.
    ## TODO: Remove line below when FIXME above is addressed.
    dimnames(ans) <- dimnames(from)
    ans
}

setAs("HDF5DatasetDump", "DelayedArray", .dump_as_HDF5Array)
setAs("ANY", "HDF5Array", .dump_as_HDF5Array)

