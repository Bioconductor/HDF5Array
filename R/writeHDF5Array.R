### =========================================================================
### writeHDF5Array() and HDF5 dump management
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 global internal counters: one for the dump files, one for the dump
### names
###
### The 2 counters are safe to use in the context of parallel execution e.g.
###
###   library(BiocParallel)
###   bplapply(1:5, function(i) .get_dump_files_global_counter(increment=TRUE))
###   bplapply(1:5, function(i) .get_dump_names_global_counter(increment=TRUE))
###

.get_dump_files_global_counter_filepath <- function()
{
    file.path(tempdir(), "HDF5Array_dump_files_global_counter")
}
 
.get_dump_names_global_counter_filepath <- function()
{
    file.path(tempdir(), "HDF5Array_dump_names_global_counter")
}

### Called by .onLoad() hook (see zzz.R file). 
init_HDF5_dump_files_global_counter <- function()
{
    filepath <- .get_dump_files_global_counter_filepath()
    init_global_counter(filepath)
}

### Called by .onLoad() hook (see zzz.R file).
init_HDF5_dump_names_global_counter <- function()
{
    filepath <- .get_dump_names_global_counter_filepath()
    init_global_counter(filepath)
}

.get_dump_files_global_counter <- function(increment=FALSE)
{
    filepath <- .get_dump_files_global_counter_filepath()
    get_global_counter(filepath, increment=increment)
}
.get_dump_names_global_counter <- function(increment=FALSE)
{
    filepath <- .get_dump_names_global_counter_filepath()
    get_global_counter(filepath, increment=increment)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getHDF5DumpFile(), lsHDF5DumpFile(), set/getHDF5DumpName()
###

.check_dump_file <- function(file)
{
    if (!isSingleString(file) || file == "")
        stop(wmsg("'file' must be a single string specifying the path ",
                  "to a new or existing HDF5 file"))
    if (file.exists(file))
        return(h5ls(file))
    h5createFile(file)
    return(NULL)
}

.check_dump_name <- function(name)
{
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the HDF5 dataset to write"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
}

.dump_settings_envir <- new.env(parent=emptyenv())

.set_dump_dir <- function(dir)
{
    if (!dir.exists(dir)) {
        if (file.exists(dir))
            stop(wmsg("\"", dir, "\" already exists and is a file, ",
                      "not a directory"))
        if (!suppressWarnings(dir.create(dir)))
            stop("cannot create directory \"", dir, "\"")
    }
    dir <- file_path_as_absolute(dir)
    assign("dir", dir, envir=.dump_settings_envir)
}

.get_dump_dir <- function()
{
    dir <- try(get("dir", envir=.dump_settings_envir), silent=TRUE)
    if (is(dir, "try-error")) {
        dir <- file.path(tempdir(), "HDF5Array_dump")
        .set_dump_dir(dir)
    }
    dir
}

### Return the user-specified file of the dump or an error if the user didn't
### specify a file.
.get_dump_specfile <- function()
{
    get("specfile", envir=.dump_settings_envir)
}

.get_dump_autofile <- function(increment=FALSE)
{
    counter <- .get_dump_files_global_counter(increment=increment)
    file <- file.path(.get_dump_dir(), sprintf("auto%05d.h5", counter))
    if (!file.exists(file))
        h5createFile(file)
    file
}

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpFile <- function(file)
{
    if (missing(file)) {
        suppressWarnings(rm(list="specfile", envir=.dump_settings_envir))
        file <- .get_dump_autofile()
        file_content <- .check_dump_file(file)
    } else {
        if (!isSingleString(file) || file == "")
            stop("'file' must be a single non-empty string")
        nc <- nchar(file)
        if (substr(file, start=nc, stop=nc) == "/") {
            if (nc >= 2L)
                file <- substr(file, start=1L, stop=nc-1L)
            .set_dump_dir(file)
            file <- .get_dump_autofile()
            file_content <- .check_dump_file(file)
        } else {
            file_content <- .check_dump_file(file)
            file <- file_path_as_absolute(file)
            assign("specfile", file, envir=.dump_settings_envir)
        }
    }
    if (is.null(file_content))
        return(invisible(file_content))
    file_content
}

### Return the *absolute path* to the dump file.
getHDF5DumpFile <- function()
{
    file <- try(.get_dump_specfile(), silent=TRUE)
    if (is(file, "try-error"))
        file <- .get_dump_autofile()
    file
}

.get_dump_file_for_use <- function()
{
    file <- try(.get_dump_specfile(), silent=TRUE)
    if (is(file, "try-error"))
        file <- .get_dump_autofile(increment=TRUE)
    file
}

### A convenience wrapper.
lsHDF5DumpFile <- function() h5ls(getHDF5DumpFile())

### Return the user-specified name of the dump or an error if the user didn't
### specify a name.
.get_dump_specname <- function()
{
    get("specname", envir=.dump_settings_envir)
}

.get_dump_autoname <- function(increment=FALSE)
{
    counter <- .get_dump_names_global_counter(increment=increment)
    sprintf("/HDF5ArrayAUTO%05d", counter)
}

setHDF5DumpName <- function(name)
{
    if (missing(name)) {
        suppressWarnings(rm(list="specname", envir=.dump_settings_envir))
        name <- .get_dump_autoname()
        return(invisible(name))
    }
    .check_dump_name(name)
    assign("specname", name, envir=.dump_settings_envir)
}

getHDF5DumpName <- function()
{
    name <- try(.get_dump_specname(), silent=TRUE)
    if (is(name, "try-error"))
        name <- .get_dump_autoname()
    name
}

.get_dump_name_for_use <- function()
{
    name <- try(.get_dump_specname(), silent=TRUE)
    if (is(name, "try-error")) {
        name <- .get_dump_autoname(increment=TRUE)
    } else {
        setHDF5DumpName()
    }
    name
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5RealizationSink objects
###
### The HDF5RealizationSink class is a concrete RealizationSink subclass that
### implements an HDF5 realization sink.
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
        file <- .get_dump_file_for_use()
    } else {
        .check_dump_file(file)
    }
    if (is.null(name)) {
        name <- .get_dump_name_for_use()
    } else {
        .check_dump_name(name)
    }
    h5createDataset2(file, name, dim, type)
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
        x_dim <- dim(x)
        sink_dim <- sink@dim
        if (is.null(offsets)) {
            stopifnot(identical(x_dim, sink_dim))
            index <- NULL
        } else {
            stopifnot(length(x_dim) == length(sink_dim))
            block_ranges <- IRanges(offsets, width=x_dim)
            index <- DelayedArray:::make_Nindex_from_block_ranges(
                                         block_ranges, sink_dim,
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

