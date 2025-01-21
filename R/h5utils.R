### =========================================================================
### Some low-level HDF5 utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used in validity methods
###

### 'path' is expected to be the **absolute** path to a local HDF5 file.
validate_h5_absolute_path <- function(path, what="'path'")
{
    if (!(isSingleString(path) && nzchar(path)))
        return(paste0(what, " must be a single non-empty string"))

    ## Check that 'path' points to an HDF5 file that is accessible.
    if (!file.exists(path))
        return(paste0(what, " (\"", path, "\") must be the path to ",
                      "an existing HDF5 file"))
    if (dir.exists(path))
        return(paste0(what, " (\"", path, "\") must be the path to ",
                      "an HDF5 file, not a directory"))
    h5_content <- try(h5ls(path), silent=TRUE)
    if (inherits(h5_content, "try-error"))
        return(paste0(what, " (\"", path, "\") doesn't seem to be ",
                      "the path to a valid HDF5 file"))
    if (path != file_path_as_absolute(path))
        return(paste0(what, " (\"", path, "\") must be the absolute ",
                      "canonical path the HDF5 file"))
    TRUE
}

validate_h5_dataset_name <- function(path, name, what="'name'")
{
    if (!(isSingleString(name) && nzchar(name)))
        return(paste0(what, " must be a single non-empty string"))

    if (!h5exists(path, name))
        return(paste0(what, " (\"", name, "\") doesn't exist ",
                      "in HDF5 file \"", path, "\""))
    if (!h5isdataset(path, name))
        return(paste0(what, " (\"", name, "\") is not a dataset ",
                      "in HDF5 file \"", path, "\""))
    h5_dim <- try(h5dim(path, name), silent=TRUE)
    if (inherits(h5_dim, "try-error"))
        return(paste0(what, " (\"", name, "\") is a dataset with ",
                      "no dimensions in HDF5 file \"", path, "\""))
    TRUE
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

### Compute the value to pass to the 'size' argument of HDF5RealizationSink(),
### which will be passed all the way down to h5createDataset2() below, and
### finally to rhdf5::h5createDataset().
compute_max_string_size <- function(x)
{
    ## We want this to work on any array-like object, not just ordinary
    ## arrays, so we must use type() instead of is.character().
    if (type(x) != "character")
        return(NULL)
    if (length(x) == 0L)
        return(0L)
    ## Calling nchar() on 'x' will trigger block processing if 'x' is a
    ## DelayedArray object, so it could take a while.
    max(nchar(x, type="bytes", keepNA=FALSE))
}

h5createDataset2 <- function(filepath, name, dim, maxdim=dim,
                             type="double", H5type=NULL, size=NULL,
                             chunkdim=dim, level=6L)
{
    stopifnot(is.numeric(dim),
              is.numeric(maxdim), length(maxdim) == length(dim))
    if (!is.null(chunkdim)) {
        stopifnot(is.numeric(chunkdim), length(chunkdim) == length(dim))
        chunkdim <- pmin(chunkdim, maxdim)
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

