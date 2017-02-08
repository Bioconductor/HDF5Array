### =========================================================================
### Manage settings for writing output to HDF5
### -------------------------------------------------------------------------
###


.check_HDF5_dump_file <- function(file)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path ",
                  "to an existing HDF5 file or to a new file"))
    if (file.exists(file))
        return(h5ls(file))
    h5createFile(file)
    return(NULL)
}

.check_HDF5_dump_name <- function(name)
{
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the HDF5 dataset to write"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
}

.HDF5_dump_settings_envir <- new.env(parent=emptyenv())

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpFile <- function(file=paste0(tempfile(), ".h5"))
{
    file_content <- .check_HDF5_dump_file(file)
    assign("file", file, envir=.HDF5_dump_settings_envir)
    if (is.null(file_content))
        return(invisible(file_content))
    return(file_content)
}

getHDF5DumpFile <- function()
    get("file", envir=.HDF5_dump_settings_envir)

### A convenience wrapper.
lsHDF5DumpFile <- function() h5ls(getHDF5DumpFile())

assign("auto_inc_ID", 0L, envir=.HDF5_dump_settings_envir)

.get_auto_inc_ID <- function()
{
    get("auto_inc_ID", envir=.HDF5_dump_settings_envir)
}

.set_HDF5_dump_name_to_next_auto_inc_ID <- function()
{
    suppressWarnings(rm(list="name", envir=.HDF5_dump_settings_envir))
    auto_inc_ID <- .get_auto_inc_ID() + 1L
    assign("auto_inc_ID", auto_inc_ID, envir=.HDF5_dump_settings_envir)
}

setHDF5DumpName <- function(name)
{
    if (missing(name))
        return(.set_HDF5_dump_name_to_next_auto_inc_ID())
    .check_HDF5_dump_name(name)
    assign("name", name, envir=.HDF5_dump_settings_envir)
}

getHDF5DumpName <- function()
{
    name <- try(get("name", envir=.HDF5_dump_settings_envir), silent=TRUE)
    if (is(name, "try-error")) {
        auto_inc_ID <- .get_auto_inc_ID()
        name <- sprintf("/HDF5ArrayAUTO%05d", auto_inc_ID)
    }
    name
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5RealizationSink objects
###

setClass("HDF5RealizationSink",
    contains="RealizationSink",
    representation(
        dim="integer",
        dimnames="list",
        type="character",  # Single string.
        file="character",  # Single string.
        name="character"   # Dataset name.
    )
)

setMethod("dimnames", "HDF5RealizationSink",
    function(x)
    {
        ans <- x@dimnames
        if (all(S4Vectors:::sapply_isNULL(ans)))
            return(NULL)
        ans
    }
)

### FIXME: Investigate the possiblity to write the dimnames to the HDF5 file.
HDF5RealizationSink <- function(dim, dimnames=NULL, type="double",
                                file=NULL, name=NULL)
{
    if (is.null(file)) {
        file <- getHDF5DumpFile()
    } else {
        .check_HDF5_dump_file(file)
    }
    if (is.null(name)) {
        use_HDF5_dump_name <- TRUE
        name <- getHDF5DumpName()
    } else {
        use_HDF5_dump_name <- FALSE
        .check_HDF5_dump_name(name)
    }
    h5createDataset2(file, name, dim, type)
    if (use_HDF5_dump_name)
        .set_HDF5_dump_name_to_next_auto_inc_ID()
    if (is.null(dimnames)) {
        dimnames <- vector("list", length(dim))
    } else {
        ## TODO: Write the dimnames to the HDF5 file.
    }
    new2("HDF5RealizationSink", dim=dim, dimnames=dimnames, type=type,
                                file=file, name=name)
}

setMethod("write_to_sink", c("array", "HDF5RealizationSink"),
    function(x, sink, offsets=NULL)
    {
        if (is.null(offsets)) {
            stopifnot(all(dim(x) == sink@dim))
            index <- NULL
        } else {
            block_ranges <- IRanges(offsets, width=dim(x))
            index <- DelayedArray:::make_subscripts_from_ranges(
                                        block_ranges,
                                        sink@dim,
                                        expand.RangeNSBS=TRUE)
        }
        h5write2(x, sink@file, sink@name, index=index)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### FIXME: This needs to propagate the dimnames. Unfortunately this is not
### possible at the moment. See FIXME right before definition of
### HDF5RealizationSink() above in this file and right before definition of
### HDF5ArraySeed() in HDF5Array-class.R about this.
.from_HDF5RealizationSink_to_HDF5ArraySeed <- function(from)
{
    HDF5ArraySeed(from@file, from@name, type=from@type)
}

setAs("HDF5RealizationSink", "HDF5ArraySeed",
      .from_HDF5RealizationSink_to_HDF5ArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeHDF5Array()
###

### Return an invisible HDF5ArraySeed object pointing to the newly written
### HDF5 dataset on disk.
### FIXME: This needs to propagate the dimnames. Unfortunately this is not
### possible at the moment. See various FIXMEs above in this file about this.
writeHDF5Array <- function(x, file, name)
{
    sink <- HDF5RealizationSink(dim(x), dimnames(x), type(x),
                                file=file, name=name)
    write_to_sink(x, sink)
    invisible(HDF5Array(as(sink, "HDF5ArraySeed")))
}

writeHDF5Dataset <- function(...)
{
    .Deprecated("writeHDF5Array")
    writeHDF5Array(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercing an array-like object to HDF5ArraySeed dumps it to disk.
###

.dump_as_HDF5ArraySeed <- function(from)
{
    sink <- HDF5RealizationSink(dim(from), dimnames(from), type(from))
    write_to_sink(from, sink)
    as(sink, "HDF5ArraySeed")
}

setAs("ANY", "HDF5ArraySeed", .dump_as_HDF5ArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to HDF5Array.
###
### Note that unless the object to coerce is a HDF5ArraySeed object, coercing
### it to HDF5Array dumps it to disk.
###

.as_HDF5Array <- function(from)
{
    ans <- HDF5Array(as(from, "HDF5ArraySeed"))
    ## Temporarily needed because coercion from HDF5RealizationSink to
    ## HDF5ArraySeed doesn't propagate the dimnames at the moment. See FIXME
    ## above.
    ## TODO: Remove line below when FIXME above is addressed.
    dimnames(ans) <- dimnames(from)
    ans
}

setAs("HDF5RealizationSink", "DelayedArray", .as_HDF5Array)
setAs("ANY", "HDF5Array", .as_HDF5Array)

### Automatic coercion method from DelayedArray to HDF5Array silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "HDF5Array", .as_HDF5Array)

