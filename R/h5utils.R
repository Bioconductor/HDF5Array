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
        stopifnot(length(chunkdim) == length(dim))
        chunkdim <- as.integer(pmin(dim, chunkdim))
    }
    chunkdim
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A thin wrapper around h5mread()
###

h5read2 <- function(filepath, name, index=NULL, as.integer=FALSE)
{
    if (!is.null(index))
        index <- DelayedArray:::expand_Nindex_RangeNSBS(index)
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    #suppressWarnings(h5read(filepath, name, index=index))
    h5mread(filepath, name, starts=index, as.integer=as.integer)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5writeDimnames() / h5readDimnames()
###

.check_dimnames <- function(dimnames, dim, name)
{
    ndim <- length(dim)
    stopifnot(is.list(dimnames), length(dimnames) <= ndim)
    not_is_NULL <- !S4Vectors:::sapply_isNULL(dimnames)
    for (along in which(not_is_NULL)) {
        dn <- dimnames[[along]]
        if (!(is.vector(dn) && is.atomic(dn) || is.factor(dn)))
            stop(wmsg("each list element in the supplied 'dimnames' ",
                      "must be NULL, or an atomic vector, or a factor"))
        if (length(dn) != dim[[along]])
            stop(wmsg("length of 'dimnames[[", along, "]]' ",
                      "(", length(dn), ") must equal the extent ",
                      "of the corresponding dimension in HDF5 ",
                      "dataset '", name, "' (", dim[[along]], ")"))
    }
    labels <- names(dimnames)
    if (!is.null(labels) && any(is.na(labels)))
        stop(wmsg("'names(dimnames)' cannot contain NAs"))
    not_is_NULL
}

.normarg_dsnames <- function(dsnames, not_is_NULL, filepath, name)
{
    stopifnot(is.character(dsnames), length(dsnames) == length(not_is_NULL))
    if (any(not_is_NULL & is.na(dsnames)))
        stop(wmsg("'dsnames' cannot have NAs associated with dimensions ",
                  "in HDF5 dataset '", name, "' for which the corresponding ",
                  "list elements in 'dimnames' are NULL"))
    for (along in which(not_is_NULL)) {
        dsname <- dsnames[[along]]
        if (h5exists(filepath, dsname))
            stop(wmsg("dataset '", dsname, "' already exists"))
    }
    dsnames[!not_is_NULL] <- NA_character_
    dsnames
}

### name:     The name of the HDF5 dataset on which to set the dimnames.
### dimnames: A list (possibly named) with 1 list element per dimension in
###           dataset 'name'.
### dsnames:  A character vector containing the names of the HDF5 datasets
###           (1 per dimension in dataset 'name') where to write the dimnames.
###           Names associated with dimensions for which the corresponding
###           list elements in 'dimnames' are NULL are ignored (hence can be
###           NAs).
h5writeDimnames <- function(filepath, name, dimnames, dsnames)
{
    ## Fail if dataset 'name' is not pristine.
    scales <- h5getdimscales(filepath, name, "dimnames")
    if (!all(is.na(scales)))
        stop(wmsg("the dimnames for HDF5 dataset '", name, "' are ",
                  "already stored in the following datasets: ",
                  paste(scales, collapse=", ")))
    labels <- h5getdimlabels(filepath, name)
    if (!is.null(labels))
        stop(wmsg("HDF5 dataset '", name, "' already has dimension labels"))

    ## Check 'dimnames'.
    dim <- h5dim(filepath, name)
    not_is_NULL <- .check_dimnames(dimnames, dim, name)

    ## Check 'dsnames'.
    dsnames <- .normarg_dsnames(dsnames, not_is_NULL, filepath, name)

    ## Write dimnames.
    for (along in which(not_is_NULL)) {
        dn <- dimnames[[along]]
        dsname <- dsnames[[along]]
        h5write(dn, filepath, dsname)
    }

    ## Attach new datasets to dimensions of dataset 'name'.
    h5setdimscales(filepath, name, dsnames, "dimnames")

    ## Set the dimension labels.
    labels <- names(dimnames)
    if (!is.null(labels) && any(nzchar(labels)))
        h5setdimlabels(filepath, name, labels)
}

h5readDimnames <- function(filepath, name)
{
    scales <- h5getdimscales(filepath, name, "dimnames")
    labels <- h5getdimlabels(filepath, name)
    if (all(is.na(scales)) && is.null(labels))
        return(NULL)
    lapply(setNames(scales, labels),
           function(scale) {
               if (is.na(scale))
                   return(NULL)
               as.character(h5mread(filepath, scale))
           })
}

