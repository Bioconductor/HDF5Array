### =========================================================================
### Some low-level HDF5 utilities
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
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
### write_h5dimnames() / read_h5dimnames()
###

.check_filepath_and_name <- function(filepath, name)
{
    ## Fail if 'name' is a Dimension Scale dataset or has Dimension Scales
    ## on it.
    if (h5isdimscale(filepath, name))
        stop(wmsg("cannot write dimnames for an HDF5 dataset '", name, "' ",
                  "that contains the dimnames of another dataset in ",
                  "the HDF5 file"))
    dimscales <- h5getdimscales(filepath, name, scalename="dimnames")
    if (!all(is.na(dimscales))) {
        dimscales <- dimscales[!is.na(dimscales)]
        stop(wmsg("the dimnames for HDF5 dataset '", name, "' are ",
                  "already stored in the following dataset(s): ",
                  paste(paste0("'", dimscales, "'"), collapse=", ")))
    }
    dimlabels <- h5getdimlabels(filepath, name)
    if (!is.null(dimlabels))
        stop(wmsg("HDF5 dataset '", name, "' already has dimension labels"))
}

.check_dimnames <- function(dimnames, filepath, name)
{
    dim <- h5dim(filepath, name)
    ndim <- length(dim)
    stopifnot(is.list(dimnames), length(dimnames) <= ndim)
    not_is_NULL <- !S4Vectors:::sapply_isNULL(dimnames)
    for (along in which(not_is_NULL)) {
        dn <- dimnames[[along]]
        if (!(is.vector(dn) && is.atomic(dn)))
            stop(wmsg("each list element in the supplied 'dimnames' ",
                      "must an atomic vector or a NULL"))
        if (length(dn) != dim[[along]])
            stop(wmsg("length of 'dimnames[[", along, "]]' ",
                      "(", length(dn), ") must equal the ",
                      "extent of the corresponding dimension in ",
                      "HDF5 dataset '", name, "' (", dim[[along]], ")"))
    }
    dimlabels <- names(dimnames)
    if (!is.null(dimlabels) && any(is.na(dimlabels)))
        stop(wmsg("'names(dimnames)' cannot contain NAs"))
    not_is_NULL
}

.normarg_group <- function(group, name)
{
    if (!isSingleStringOrNA(group))
        stop(wmsg("'group' must be a single string or NA"))
    if (is.na(group))
        group <- sprintf("%s/.%s_dimnames", dirname(name), basename(name))
    group
}

.normarg_dimscales <- function(dimscales, group, not_is_NULL, filepath, name)
{
    ndim <- length(not_is_NULL)
    if (is.null(dimscales)) {
        ## Generate automatic dataset names.
        digits <- as.integer(log10(ndim + 0.5)) + 1L
        fmt <- paste0("%0", digits, "d")
        dimscales <- sprintf(fmt, seq_len(ndim))
    } else {
        if (!is.character(dimscales) || length(dimscales) != ndim)
            stop(wmsg("'dimscales' must be a character vector containing ",
                      "the names of the HDF5 datasets (1 per list element ",
                      "in 'dimnames') where to write the dimnames"))
        if (any(not_is_NULL & is.na(dimscales)))
            stop(wmsg("'dimscales' cannot have NAs associated with ",
                      "list elements in 'dimnames' that are not NULL"))
    }
    if (nzchar(group))
        dimscales <- paste0(group, "/", dimscales)
    dimscales[!not_is_NULL] <- NA_character_
    for (along in which(not_is_NULL)) {
        dimscale <- dimscales[[along]]
        if (h5exists(filepath, dimscale))
            stop(wmsg("dataset '", dimscale, "' already exists"))
    }
    dimscales
}

### Exported!
### dimnames:  A list (possibly named) with 1 list element per dimension in
###            dataset 'name'.
### name:      The name of the HDF5 dataset on which to set the dimnames.
### group:     The name of the HDF5 group where to write the dimnames.
###            If NA, the group name is automatically generated from 'name'.
###            An empty string ("") means that no group should be used.
###            Otherwise, the names in 'dimscales' must be relative to the
###            specified group name.
### dimscales: A character vector containing the names of the HDF5 datasets
###            (1 per list element in 'dimnames') where to write the dimnames.
###            Names associated with NULL list elements in 'dimnames' are
###            ignored.
write_h5dimnames <- function(dimnames, filepath, name, group=NA, dimscales=NULL)
{
    ## 1. Lots of checks.

    ## Before we start writing to the file we want some guarantees that
    ## the full operation will succeed. The checks we make access the file
    ## in read-only mode.
    .check_filepath_and_name(filepath, name)

    not_is_NULL <- .check_dimnames(dimnames, filepath, name)

    group <- .normarg_group(group, name)

    dimscales <- .normarg_dimscales(dimscales, group, not_is_NULL,
                                    filepath, name)

    ## 2. Write to the HDF5 file.

    ## Create group if needed.
    if (!is.na(group) && !h5exists(filepath, group))
        h5createGroup(filepath, group)

    ## Write dimnames.
    for (along in which(not_is_NULL)) {
        dn <- dimnames[[along]]
        dimscale <- dimscales[[along]]
        h5write(dn, filepath, dimscale)
    }

    ## Attach new datasets to dimensions of dataset 'name'.
    h5setdimscales(filepath, name, dimscales, scalename="dimnames")

    ## Set the dimension labels.
    dimlabels <- names(dimnames)
    if (!is.null(dimlabels) && any(nzchar(dimlabels)))
        h5setdimlabels(filepath, name, dimlabels)
}

### Exported!
read_h5dimnames <- function(filepath, name)
{
    dimscales <- h5getdimscales(filepath, name, scalename="dimnames")
    dimlabels <- h5getdimlabels(filepath, name)
    if (all(is.na(dimscales)) && is.null(dimlabels))
        return(NULL)
    lapply(setNames(dimscales, dimlabels),
           function(dimscale) {
               if (is.na(dimscale))
                   return(NULL)
               DelayedArray:::set_dim(h5mread(filepath, dimscale), NULL)
           })
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### validate_lengths_of_h5dimnames()
###

validate_lengths_of_h5dimnames <- function(filepath, name)
{
    dimscales <- h5getdimscales(filepath, name, scalename="dimnames")
    dim <- h5dim(filepath, name)
    for (along in which(!is.na(dimscales))) {
        dimscale <- dimscales[[along]]
        dimscale_len <- prod(h5dim(filepath, dimscale))
        if (dimscale_len != dim[[along]])
            return(paste0("HDF5 dataset '", name, "' has invalid ",
                          "dimnames: length of dataset '", dimscale, "' ",
                          "(", dimscale_len, ") is not equal to ",
                          "the extent of dimension ", along, " in ",
                          "dataset '", name, "' (", dim[[along]], ")"))
    }
    TRUE
}

