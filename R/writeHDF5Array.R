### =========================================================================
### writeHDF5Array() and HDF5 dump management
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
### Coercing a HDF5RealizationSink object.
###

### FIXME: This coercion needs to propagate the dimnames *thru* the HDF5 file.
### For more details about this, see FIXME right before definition of
### HDF5RealizationSink() above in this file and right before definition of
### HDF5ArraySeed() in HDF5Array-class.R.
setAs("HDF5RealizationSink", "HDF5ArraySeed",
    function(from) HDF5ArraySeed(from@file, from@name, type=from@type)
)

### Note that this coercion currently drops the dimnames but will naturally
### propagate them when coercion from HDF5RealizationSink to HDF5ArraySeed
### propagates them. See FIXME above.
setAs("HDF5RealizationSink", "HDF5Array",
    function(from) HDF5Array(as(from, "HDF5ArraySeed"))
)

setAs("HDF5RealizationSink", "DelayedArray",
    function(from)
    {
        ans <- HDF5Array(as(from, "HDF5ArraySeed"))
        ## Temporarily needed because coercion from HDF5RealizationSink to
        ## HDF5ArraySeed does not propagate the dimnames at the moment. See
        ## FIXME above.
        ## TODO: Remove line below when FIXME above is addressed.
        dimnames(ans) <- dimnames(from)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeHDF5Array()
###

### Write the dataset to the current dump if 'file' and 'name' are not
### specified.
### Return a HDF5Array object pointing to the newly written HDF5 dataset on
### disk.
### FIXME: This needs to write the dimnames to the file. See various FIXMEs
### above in this file about this.
writeHDF5Array <- function(x, file=NULL, name=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    sink <- HDF5RealizationSink(dim(x), dimnames(x), type(x),
                                file=file, name=name)
    if (verbose) {
        old_verbose <- DelayedArray:::set_verbose_block_processing(verbose)
        on.exit(DelayedArray:::set_verbose_block_processing(old_verbose))
    }
    write_to_sink(x, sink)
    as(sink, "HDF5Array")
}

writeHDF5Dataset <- function(...)
{
    .Deprecated("writeHDF5Array")
    writeHDF5Array(...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to HDF5Array.
###
### The methods below write the object to disk. Note that coercion from
### HDF5RealizationSink to HDF5Array is already taken care of by the specific
### method above and doesn't write anything to disk. So coercing to HDF5Array
### in general writes the object to disk *except* when the object to coerce is
### a HDF5RealizationSink object.
###

.as_HDF5Array <- function(from) writeHDF5Array(from)  # write to current dump

setAs("ANY", "HDF5Array", .as_HDF5Array)

### Automatic coercion method from DelayedArray to HDF5Array silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "HDF5Array", .as_HDF5Array)

