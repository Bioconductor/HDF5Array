### =========================================================================
### Some low-level HDF5 utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5exists()
###

h5exists <- function(filepath, name)
{
    fid <- H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    H5Lexists(fid, name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5isdataset() and h5isgroup()
###

h5isdataset <- function(filepath, name)
{
    fid <- H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    did <- try(H5Dopen(fid, name), silent=TRUE)
    ans <- !inherits(did, "try-error")
    if (ans)
        H5Dclose(did)
    ans
}

h5isgroup <- function(filepath, name)
{
    fid <- H5Fopen(filepath, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(fid))
    gid <- try(H5Gopen(fid, name), silent=TRUE)
    ans <- !inherits(gid, "try-error")
    if (ans)
        H5Gclose(gid)
    ans
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
    if (is(filepath, "H5File")) {
        fid <- as(filepath, "H5IdComponent")
    } else {
        fid <- H5Fopen(filepath, flags="H5F_ACC_RDONLY")
        on.exit(H5Fclose(fid))
    }
    gid <- H5Gopen(fid, group)
    on.exit(H5Gclose(gid), add=TRUE)
    H5Dopen(gid, name)
}

dim_as_integer <- function(dim, filepath, name, what="HDF5 dataset")
{
    if (is.integer(dim))
        return(dim)
    if (any(dim > .Machine$integer.max)) {
        dim_in1string <- paste0(dim, collapse=" x ")
        if (is(filepath, "H5File"))
            filepath <- path(filepath)
        stop(wmsg("Dimensions of ", what, " are too big: ", dim_in1string),
             "\n\n  ",
             wmsg("(This error is about HDF5 dataset '", name, "' ",
                  "from file '", filepath, "'.)"),
             "\n\n  ",
             wmsg("Please note that the HDF5Array package only ",
                  "supports datasets where each dimension is ",
                  "<= '.Machine$integer.max' (= 2**31 - 1)."))
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
        dim <- dim_as_integer(dim, filepath, name)
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
    chunkdim <- dim_as_integer(chunkdim, filepath, name,
                               what="HDF5 dataset chunks")
    if (adjust) {
        dim <- h5dim(filepath, name, as.integer=FALSE)
        ## A sanity check that should never fail.
        stopifnot(length(chunkdim) == length(dim))
        chunkdim <- as.integer(pmin(dim, chunkdim))
    }
    chunkdim
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normarg_h5_filepath() and normarg_h5_name()
###

normarg_h5_filepath <- function(path, what1="'filepath'", what2="the dataset")
{
    if (!isSingleString(path))
        stop(wmsg(what1, " must be a single string specifying the path ",
                  "to the HDF5 file where ", what2, " is located"))
    file_path_as_absolute(path)  # return absolute path in canonical form
}

normarg_h5_name <- function(name, what1="'name'",
                                  what2="the name of a dataset",
                                  what3="")
{
    if (!isSingleString(name))
        stop(wmsg(what1, " must be a single string specifying ",
                  what2, " in the HDF5 file", what3))
    if (name == "")
        stop(wmsg(what1, " cannot be the empty string"))
    if (substr(name, start=1L, stop=1L) == "/") {
        name <- sub("^/*", "/", name)  # only keep first leading slash
    } else {
        name <- paste0("/", name)
    }
    name
}


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

