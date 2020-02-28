### =========================================================================
### HDF5 dump management
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
### Normalization (with basic checking) of an HDF5 file path or dataset name
###

### Return the *absolute path* to the dump file.
### Has the side effect of creating the file as an empty HDF5 file if it does
### not exist yet.
normalize_dump_filepath <- function(filepath)
{
    if (!isSingleString(filepath) || filepath == "")
        stop(wmsg("'filepath' must be a non-empty string specifying the ",
                  "path to a new or existing HDF5 file"))
    if (!file.exists(filepath))
        h5createFile(filepath)
    file_path_as_absolute(filepath)
}

normalize_dump_name <- function(name)
{
    if (!isSingleString(name) || name == "")
        stop(wmsg("'name' must be a non-empty string specifying the name ",
                  "of the HDF5 dataset to write"))
    trim_trailing_slashes(name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Very low-level stuff used in this file only
###

.dump_settings_envir <- new.env(parent=emptyenv())

### Create directory 'dir' if it doesn't exist yet.
.set_dump_dir <- function(dir)
{
    ## Even though file_path_as_absolute() will trim the trailing slashes,
    ## we need to do this early. Otherwise, checking for the existence of a
    ## file of the same name as the to-be-created directory will fail.
    if (nchar(dir) > 1L)
        dir <- trim_trailing_slashes(dir)
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

.set_dump_autofiles_mode <- function()
{
    suppressWarnings(rm(list="specfile", envir=.dump_settings_envir))
}

### Create file as an empty HDF5 file if it doesn't exist yet.
.set_dump_specfile <- function(filepath)
{
    filepath <- normalize_dump_filepath(filepath)
    assign("specfile", filepath, envir=.dump_settings_envir)
}

.set_dump_autonames_mode <- function()
{
    suppressWarnings(rm(list="specname", envir=.dump_settings_envir))
}

.set_dump_specname <- function(name)
{
    assign("specname", name, envir=.dump_settings_envir)
}

### Return the user-specified file of the dump or an error if the user didn't
### specify a file.
.get_dump_specfile <- function()
{
    get("specfile", envir=.dump_settings_envir)
}

.get_dump_autoname <- function(increment=FALSE)
{
    counter <- .get_dump_names_global_counter(increment=increment)
    sprintf("/HDF5ArrayAUTO%05d", counter)
}

### Return the user-specified name of the dump or an error if the user didn't
### specify a name.
.get_dump_specname <- function()
{
    get("specname", envir=.dump_settings_envir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get/setHDF5DumpDir()
###

getHDF5DumpDir <- function()
{
    get("dir", envir=.dump_settings_envir)
}

### Create auto file as an empty HDF5 file if it doesn't exist yet.
.get_dump_autofile <- function(increment=FALSE)
{
    counter <- .get_dump_files_global_counter(increment=increment)
    filepath <- file.path(getHDF5DumpDir(), sprintf("auto%05d.h5", counter))
    if (!file.exists(filepath))
        h5createFile(filepath)
    filepath
}

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpDir <- function(dir)
{
    if (missing(dir)) {
        dir <- file.path(tempdir(), "HDF5Array_dump")
    } else if (!isSingleString(dir) || dir == "") {
        stop(wmsg("'dir' must be a non-empty string specifying the path ",
                  "to a new or existing directory"))
    }
    dir <- .set_dump_dir(dir)
    .set_dump_autofiles_mode()
    .get_dump_autofile()
    invisible(dir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getHDF5DumpFile()
###

### Set the current HDF5 dump file. Create it as an empty HDF5 file if it
### doesn't exist yet.
setHDF5DumpFile <- function(filepath)
{
    if (missing(filepath)) {
        .set_dump_autofiles_mode()
        filepath <- .get_dump_autofile()
    } else {
        if (!isSingleString(filepath) || filepath == "")
            stop("'filepath' must be a non-empty string")
        if (has_trailing_slash(filepath)) {
            setHDF5DumpDir(filepath)
            filepath <- .get_dump_autofile()
        } else {
            filepath <- .set_dump_specfile(filepath)
        }
    }
    file_content <- h5ls(filepath)
    if (nrow(file_content) == 0L)
        return(invisible(file_content))
    file_content
}

### Return the *absolute path* to the dump file.
getHDF5DumpFile <- function(for.use=FALSE)
{
    if (!isTRUEorFALSE(for.use))
        stop("'for.use' must be TRUE or FALSE")
    filepath <- try(.get_dump_specfile(), silent=TRUE)
    if (is(filepath, "try-error"))
        filepath <- .get_dump_autofile(increment=for.use)
    filepath
}

### A convenience wrapper.
lsHDF5DumpFile <- function() h5ls(getHDF5DumpFile())


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getHDF5DumpName()
###

setHDF5DumpName <- function(name)
{
    if (missing(name)) {
        .set_dump_autonames_mode()
        name <- .get_dump_autoname()
        return(invisible(name))
    }
    name <- normalize_dump_name(name)
    .set_dump_specname(name)
}

getHDF5DumpName <- function(for.use=FALSE)
{
    if (!isTRUEorFALSE(for.use))
        stop("'for.use' must be TRUE or FALSE")
    name <- try(.get_dump_specname(), silent=TRUE)
    if (is(name, "try-error")) {
        name <- .get_dump_autoname(increment=for.use)
    } else if (for.use) {
        ## If the dump file is a user-specified file, we switch back to
        ## automatic dump names.
        filepath <- try(.get_dump_specfile(), silent=TRUE)
        if (!is(filepath, "try-error"))
            .set_dump_autonames_mode()
    }
    name
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getHDF5DumpChunkLength() and set/getHDF5DumpChunkShape()
###

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpChunkLength <- function(length=1000000L)
{
    assign("chunk_len", length, envir=.dump_settings_envir)
}

getHDF5DumpChunkLength <- function()
{
    get("chunk_len", envir=.dump_settings_envir)
}

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpChunkShape <- function(shape="scale")
{
    assign("chunk_shape", shape, envir=.dump_settings_envir)
}

getHDF5DumpChunkShape <- function()
{
    get("chunk_shape", envir=.dump_settings_envir)
}

### TODO: Replace 'dim' argument with 'x'. This will allow the default
### **write** chunk dim to be set to the **read** chunk dim (i.e. to
### 'chunkdim(x)') if this is known (i.e. if 'chunkdim(x)' is not NULL).
getHDF5DumpChunkDim <- function(dim)
{
    chunk_len <- getHDF5DumpChunkLength()
    chunk_shape <- getHDF5DumpChunkShape()
    makeCappedVolumeBox(chunk_len, dim, chunk_shape)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### set/getHDF5DumpCompressionLevel
###

normalize_compression_level <- function(level)
{
    if (!isSingleNumber(level))
        stop("'level' must be a single number")
    if (!is.integer(level))
        level <- as.integer(level)
    if (level < 0L || level > 9L)
        stop(wmsg("'level' must be between 0 (no compression) ",
                  "and 9 (highest and slowest compression)"))
    level
}

### Called by .onLoad() hook (see zzz.R file).
setHDF5DumpCompressionLevel <- function(level=6L)
{
    level <- normalize_compression_level(level)
    assign("level", level, envir=.dump_settings_envir)
}

getHDF5DumpCompressionLevel <- function()
{
    get("level", envir=.dump_settings_envir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Dump log
###

### Called by .onLoad() hook (see zzz.R file).
get_HDF5_dump_logfile <- function()
{
    file.path(tempdir(), "HDF5Array_dump_log")
}

.get_dataset_creation_global_counter_filepath <- function()
{
    file.path(tempdir(), "HDF5Array_dataset_creation_global_counter")
}

### Called by .onLoad() hook (see zzz.R file).
init_HDF5_dataset_creation_global_counter <- function()
{
    filepath <- .get_dataset_creation_global_counter_filepath()
    init_global_counter(filepath)
}

.get_dataset_creation_global_counter <- function(increment=FALSE)
{
    filepath <- .get_dataset_creation_global_counter_filepath()
    get_global_counter(filepath, increment=increment)
}

### Use a lock mechanism so is safe to use in the context of parallel
### execution.
appendDatasetCreationToHDF5DumpLog <- function(filepath, name, dim, type,
                                               chunkdim, level)
{
    logfile <- get_HDF5_dump_logfile()
    locked_path <- lock_file(logfile)
    on.exit(unlock_file(logfile))
    counter <- .get_dataset_creation_global_counter(increment=TRUE)
    dims <- paste0(dim, collapse="x")
    chunkdims <- paste0(chunkdim, collapse="x")
    cat(as.character(Sys.time()), counter,
        filepath, name, dims, type, chunkdims, level,
        sep="\t", file=locked_path, append=TRUE)
    cat("\n", file=locked_path, append=TRUE)
}

showHDF5DumpLog <- function()
{
    COLNAMES <- c("time", "counter",
                  "filepath", "name", "dims", "type", "chunkdims", "level")
    ## The nb of lines in the log file is the current value of the dataset
    ## creation counter minus one.
    counter <- .get_dataset_creation_global_counter()
    if (counter == 1L) {
        dump_log <- data.frame(time=character(0),
                               counter=integer(0),
                               filepath=character(0),
                               name=character(0),
                               dims=character(0),
                               type=character(0),
                               chunkdims=character(0),
                               level=integer(0),
                               stringsAsFactors=FALSE)
    } else {
        logfile <- get_HDF5_dump_logfile()
        locked_path <- lock_file(logfile)
        on.exit(unlock_file(logfile))
        dump_log <- read.table(locked_path,
                               sep="\t", stringsAsFactors=FALSE)
        colnames(dump_log) <- COLNAMES
        fmt <- "[%s] #%d In file '%s': creation of dataset '%s' (%s:%s, chunkdims=%s, level=%d)"
        message(paste0(sprintf(fmt,
                               dump_log$time, dump_log$counter,
                               dump_log$filepath, dump_log$name,
                               dump_log$dims, dump_log$type,
                               dump_log$chunkdims, dump_log$level),
                       "\n"),
                appendLF=FALSE)
    }
    invisible(dump_log)
}

create_and_log_HDF5_dataset <- function(filepath, name, dim, maxdim=dim,
                                        type="double", H5type=NULL, size=NULL,
                                        chunkdim=dim, level=6L)
{
    h5createDataset2(filepath, name, dim, maxdim=maxdim,
                     type=type, H5type=H5type, size=size,
                     chunkdim=chunkdim, level=level)
    appendDatasetCreationToHDF5DumpLog(filepath, name, dim,
                                       type, chunkdim, level)
}

