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
normalize_dump_file <- function(file)
{
    if (!isSingleString(file) || file == "")
        stop(wmsg("'file' must be a non-empty string specifying the path ",
                  "to a new or existing HDF5 file"))
    if (!file.exists(file))
        h5createFile(file)
    file_path_as_absolute(file)
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
.set_dump_specfile <- function(file)
{
    file <- normalize_dump_file(file)
    assign("specfile", file, envir=.dump_settings_envir)
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
    file <- file.path(getHDF5DumpDir(), sprintf("auto%05d.h5", counter))
    if (!file.exists(file))
        h5createFile(file)
    file
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
setHDF5DumpFile <- function(file)
{
    if (missing(file)) {
        .set_dump_autofiles_mode()
        file <- .get_dump_autofile()
    } else {
        if (!isSingleString(file) || file == "")
            stop("'file' must be a non-empty string")
        if (has_trailing_slash(file)) {
            setHDF5DumpDir(file)
            file <- .get_dump_autofile()
        } else {
            file <- .set_dump_specfile(file)
        }
    }
    file_content <- h5ls(file)
    if (nrow(file_content) == 0L)
        return(invisible(file_content))
    file_content
}

### Return the *absolute path* to the dump file.
getHDF5DumpFile <- function(for.use=FALSE)
{
    if (!isTRUEorFALSE(for.use))
        stop("'for.use' must be TRUE or FALSE")
    file <- try(.get_dump_specfile(), silent=TRUE)
    if (is(file, "try-error"))
        file <- .get_dump_autofile(increment=for.use)
    file
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
        file <- try(.get_dump_specfile(), silent=TRUE)
        if (!is(file, "try-error"))
            .set_dump_autonames_mode()
    }
    name
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get the chunk dimensions
###

getHDF5DumpChunkDim <- function(dim, type, ratio=75)
{
    max_block_len <- DelayedArray:::get_max_block_length(type)
    chunk_len <- as.integer(ceiling(max_block_len / ratio))
    ## 'max_block_len' must be a multiple of 'chunk_len'.
    stopifnot(max_block_len %% chunk_len == 0L)
    chunks <- ArrayBlocks(dim, chunk_len)
    chunk_dim <- chunks@dim
    ndim <- length(chunk_dim)
    if (chunks@N > ndim)
        return(chunk_dim)
    chunk_dim[[chunks@N]] <- chunks@by
    if (chunks@N == ndim)
        return(chunk_dim)
    chunk_dim[(chunks@N+1L):ndim] <- 1L
    chunk_dim
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
appendDatasetCreationToHDF5DumpLog <- function(file, name, dim, type,
                                               chunk_dim, level)
{
    logfile <- get_HDF5_dump_logfile()
    locked_path <- lock_file(logfile)
    on.exit(unlock_file(logfile))
    counter <- .get_dataset_creation_global_counter(increment=TRUE)
    dims <- paste0(dim, collapse="x")
    chunk_dims <- paste0(chunk_dim, collapse="x")
    cat(as.character(Sys.time()), counter,
        file, name, dims, type, chunk_dims, level,
        sep="\t", file=locked_path, append=TRUE)
    cat("\n", file=locked_path, append=TRUE)
}

showHDF5DumpLog <- function()
{
    COLNAMES <- c("time", "counter",
                  "file", "name", "dims", "type", "chunk_dims", "level")
    ## The nb of lines in the log file is the current value of the dataset
    ## creation counter minus one.
    counter <- .get_dataset_creation_global_counter()
    if (counter == 1L) {
        dump_log <- data.frame(time=character(0),
                               counter=integer(0),
                               file=character(0),
                               name=character(0),
                               dims=character(0),
                               type=character(0),
                               chunk_dims=character(0),
                               level=integer(0),
                               stringsAsFactors=FALSE)
    } else {
        logfile <- get_HDF5_dump_logfile()
        locked_path <- lock_file(logfile)
        on.exit(unlock_file(logfile))
        dump_log <- read.table(locked_path,
                               sep="\t", stringsAsFactors=FALSE)
        colnames(dump_log) <- COLNAMES
        fmt <- "[%s] #%d In file '%s': creation of dataset '%s' (%s:%s, chunk_dims=%s, level=%d)"
        message(paste0(sprintf(fmt,
                               dump_log$time, dump_log$counter,
                               dump_log$file, dump_log$name,
                               dump_log$dims, dump_log$type,
                               dump_log$chunk_dims, dump_log$level),
                       "\n"),
                appendLF=FALSE)
    }
    invisible(dump_log)
}

