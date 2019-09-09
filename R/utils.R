### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5exists()
###

h5exists <- function(filepath, name)
{
    fid <- H5Fopen(filepath)
    on.exit(H5Fclose(fid))
    H5Lexists(fid, name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5dim() and h5chunkdim()
###

### Return an object of class H5IdComponent representing an H5 dataset ID.
.get_h5dataset <- function(filepath, name)
{
    if (substr(name, 1L, 1L) != "/")
        name <- paste0("/", name)
    group <- gsub("(.*/)[^/]*$", "\\1", name)
    name <- gsub(".*/([^/]*)$", "\\1", name)
    fid <- H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    gid <- H5Gopen(fid, group)
    on.exit(H5Gclose(gid), add=TRUE)
    H5Dopen(gid, name)
}

.dim_as_integer <- function(dim, filepath, name, what="dimensions")
{
    if (is.integer(dim))
        return(dim)
    if (any(dim > .Machine$integer.max)) {
        dim_in1string <- paste0(dim, collapse=" x ")
        stop(wmsg("The ", what, " (", dim_in1string, ") ",
                  "of HDF5 dataset '", name, "' ",
                  "from file '", filepath, "' are too big.\n\n",
                  "The HDF5Array package only supports datasets with ",
                  "all ", what, " <= 2^31-1 (= ", .Machine$integer.max, ") ",
                  "at the moment."))
    }
    as.integer(dim)
}

### The TENxMatrixSeed() constructor calls h5dim() with 'as.integer=FALSE'
### in order to get the dimension of a monodimensional array of length >= 2^31.
h5dim <- function(filepath, name, as.integer=TRUE)
{
    did <- .get_h5dataset(filepath, name)
    on.exit(H5Dclose(did), add=TRUE)
    sid <- H5Dget_space(did)
    on.exit(H5Sclose(sid), add=TRUE)
    dim <- H5Sget_simple_extent_dims(sid)$size
    if (as.integer)
        dim <- .dim_as_integer(dim, filepath, name)
    dim
}

### Return NULL or an integer vector parallel to 'h5dim(filepath, name)'.
h5chunkdim <- function(filepath, name, adjust=FALSE)
{
    did <- .get_h5dataset(filepath, name)
    on.exit(H5Dclose(did), add=TRUE)
    pid <- H5Dget_create_plist(did)
    on.exit(H5Pclose(pid), add=TRUE)
    if (H5Pget_layout(pid) != "H5D_CHUNKED")
        return(NULL)
    ## We use rev() to invert the order of the dimensions returned by
    ## H5Pget_chunk(). It seems that H5Pget_chunk() should take care of
    ## this though, for consistency with how rhdf5 handles the order of the
    ## dimensions everywhere else (e.g. see ?H5Sget_simple_extent_dims).
    chunkdim <- rev(H5Pget_chunk(pid))
    chunkdim <- .dim_as_integer(chunkdim, filepath, name,
                                what="chunk dimensions")
    if (adjust) {
        dim <- h5dim(filepath, name, as.integer=FALSE)
        ## A sanity check that should never fail.
        if (length(chunkdim) != length(dim) || any(chunkdim > pmax(dim, 1L)))
            stop(wmsg("the chunk dimensions of HDF5 dataset '", name, "' ",
                      "from file '", filepath, "' are not compatible with ",
                      "its dimensions"))
        chunkdim <- as.integer(pmin(dim, chunkdim))
    }
    chunkdim
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A thin wrapper around h5mread()
###

h5read2 <- function(filepath, name, index=NULL)
{
    if (!is.null(index))
        index <- DelayedArray:::expand_Nindex_RangeNSBS(index)
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    #suppressWarnings(h5read(filepath, name, index=index))
    h5mread(filepath, name, starts=index)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Manipulate one-dimensional HDF5 datasets
###

### Length of a one-dimensional HDF5 dataset.
### Return the length as a single integer (if < 2^31) or numeric (if >= 2^31).
h5length <- function(filepath, name)
{
    len <- h5dim(filepath, name, as.integer=FALSE)
    stopifnot(length(len) == 1L)
    len
}

### Append data to a one-dimensional HDF5 dataset.
### Return the length of the extended dataset.
h5append <- function(data, filepath, name)
{
    old_len <- as.numeric(h5length(filepath, name))
    data_len <- length(data)
    new_len <- old_len + data_len
    h5set_extent(filepath, name, new_len)
    h5write(data, filepath, name, start=old_len+1, count=data_len)
    new_len
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A simple wrapper around rhdf5::h5createDataset()
###

### A simple wrapper around rhdf5::h5createDataset().
h5createDataset2 <- function(filepath, name, dim, maxdim=dim,
                             type="double", H5type=NULL,
                             chunkdim=dim, level=6L)
{
    if (type == "character") {
        size <- max(nchar(name, type="width"))
    } else {
        size <- NULL
    }
    ## If h5createDataset() fails, it will leave an HDF5 file handle opened.
    ## Calling H5close() will close all opened HDF5 object handles.
    #on.exit(H5close())
    ok <- h5createDataset(filepath, name, dim, maxdims=maxdim,
                          storage.mode=type, H5type=H5type, size=size,
                          chunk=chunkdim, level=level)
    if (!ok)
        stop(wmsg("failed to create dataset '", name, "' ",
                  "in file '", filepath, "'"), call.=FALSE)
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

