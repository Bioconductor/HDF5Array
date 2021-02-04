### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


normarg_path <- function(path, what1, what2)
{
    if (!isSingleString(path))
        stop(wmsg(what1, " must be a single string specifying the path ",
                  "to the file where the ", what2, " is located"))
    file_path_as_absolute(path)  # return absolute path in canonical form
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

add_prefix_to_basename <- function(name, prefix=".")
{
    stopifnot(isSingleString(name), isSingleString(prefix))
    slash_idx <- which(safeExplode(name) == "/")
    if (length(slash_idx) == 0L) {
        dname <- ""
        bname <- name
    } else {
        last_slash_idx <- max(slash_idx)
        dname <- substr(name, start=1L, stop=last_slash_idx)
        bname <- substr(name, start=last_slash_idx+1L, stop=nchar(name))
    }
    paste0(dname, prefix, bname)
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

