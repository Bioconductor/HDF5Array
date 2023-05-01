### =========================================================================
### h5writeDimnames() and other related functions
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_h5dimnames() / set_h5dimnames()
###

get_h5dimnames <- function(filepath, name)
{
    h5getdimscales(filepath, name, scalename="dimnames")
}

### Fail if 'name' is a Dimension Scale dataset or has Dimension Scales on it.
.check_filepath_and_name <- function(filepath, name)
{
    if (h5isdimscale(filepath, name))
        stop(wmsg("HDF5 dataset \"", name, "\" contains the dimnames for ",
                  "another dataset in the HDF5 file so dimnames cannot ",
                  "be set on it"))
    current_h5dimnames <- get_h5dimnames(filepath, name)
    if (!all(is.na(current_h5dimnames))) {
        ds <- current_h5dimnames[!is.na(current_h5dimnames)]
        stop(wmsg("the dimnames for HDF5 dataset \"", name, "\" are already ",
                  "stored in HDF5 file \"", filepath, "\" (in dataset(s): ",
                  paste(paste0("\"", ds, "\""), collapse=", "), ")"))
    }
    dimlabels <- h5getdimlabels(filepath, name)
    if (!is.null(dimlabels))
        stop(wmsg("HDF5 dataset \"", name, "\" already has dimension labels"))
}

.validate_h5dimnames_lengths <- function(filepath, name, h5dimnames)
{
    dim <- h5dim(filepath, name)
    for (along in which(!is.na(h5dimnames))) {
        h5dn <- h5dimnames[[along]]
        h5dn_len <- prod(h5dim(filepath, h5dn))
        if (h5dn_len != dim[[along]])
            return(paste0("length of HDF5 dataset \"", h5dn, "\" ",
                          "(", h5dn_len, ") is not equal to the ",
                          "extent of dimension ", along, " in HDF5 ",
                          "dataset \"", name, "\" (", dim[[along]], ")"))
    }
    TRUE
}

.check_h5dimnames <- function(filepath, name, h5dimnames)
{
    dim <- h5dim(filepath, name)
    ndim <- length(dim)
    if (!is.character(h5dimnames))
        stop(wmsg("'h5dimnames' must be a character vector containing ",
                  "the names of the HDF5 datasets to set as the ",
                  "dimnames of dataset \"", name, "\" (one per dimension ",
                  "in \"", name, "\")"))
    if (length(h5dimnames) > ndim)
        stop(wmsg("length of 'h5dimnames' must equal the number of ",
                  "dimensions (", ndim, ") in HDF5 dataset \"", name, "\""))
    for (along in which(!is.na(h5dimnames))) {
        h5dn <- h5dimnames[[along]]
        if (!h5exists(filepath, h5dn))
            stop(wmsg("HDF5 dataset \"", h5dn, "\" does not exist ",
                      "in this HDF5 file"))
    }
    msg <- .validate_h5dimnames_lengths(filepath, name, h5dimnames)
    if (!isTRUE(msg))
        stop(wmsg("invalid 'h5dimnames': ", msg))
}

set_h5dimnames <- function(filepath, name, h5dimnames, dry.run=FALSE)
{
    .check_filepath_and_name(filepath, name)
    .check_h5dimnames(filepath, name, h5dimnames)
    h5setdimscales(filepath, name, dimscales=h5dimnames,
                   scalename="dimnames", dry.run=dry.run)
    invisible(NULL)
}

### For internal use only (not exported).
### FIXME: Make this work on an H5File object!
validate_lengths_of_h5dimnames <- function(filepath, name)
{
    if (is(filepath, "H5File"))
        return(TRUE)  # we cheat for now (fix this!)
    h5dimnames <- get_h5dimnames(filepath, name)
    msg <- .validate_h5dimnames_lengths(filepath, name, h5dimnames)
    if (!isTRUE(msg))
        return(paste0("invalid dimnames found in HDF5 file \"", filepath, "\" ",
                      "for dataset \"", name, "\": ", msg))
    TRUE
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5writeDimnames() / h5readDimnames()
###

.check_dimnames <- function(dimnames, filepath, name)
{
    if (!is.list(dimnames))
        stop(wmsg("'dimnames' must be a list"))
    dim <- h5dim(filepath, name)
    ndim <- length(dim)
    if (length(dimnames) > ndim)
        stop(wmsg("'dimnames' cannot have more list elements than ",
                  "the number of dimensions in dataset \"", name,"\""))
    not_NULL <- !S4Vectors:::sapply_isNULL(dimnames)
    for (along in which(not_NULL)) {
        dn <- dimnames[[along]]
        if (!(is.vector(dn) && is.atomic(dn)))
            stop(wmsg("each list element in the supplied 'dimnames' ",
                      "must an atomic vector or a NULL"))
        if (length(dn) != dim[[along]])
            stop(wmsg("length of 'dimnames[[", along, "]]' ",
                      "(", length(dn), ") must equal the ",
                      "extent of the corresponding dimension in ",
                      "HDF5 dataset \"", name, "\" (", dim[[along]], ")"))
    }
    dimlabels <- names(dimnames)
    if (!is.null(dimlabels) && any(is.na(dimlabels)))
        stop(wmsg("'names(dimnames)' cannot contain NAs"))
    not_NULL
}

.normarg_group <- function(group, name)
{
    if (!isSingleStringOrNA(group))
        stop(wmsg("'group' must be a single string or NA"))
    if (is.na(group)) {
        group <- add_prefix_to_basename(name, prefix=".")
        group <- paste0(group, "_dimnames")
    }
    group
}

.normarg_h5dimnames <- function(h5dimnames, group, not_NULL, filepath, name)
{
    ndim <- length(not_NULL)
    if (is.null(h5dimnames)) {
        ## Generate automatic dataset names.
        digits <- as.integer(log10(ndim + 0.5)) + 1L
        fmt <- paste0("%0", digits, "d")
        h5dimnames <- sprintf(fmt, seq_len(ndim))
    } else {
        if (!is.character(h5dimnames) || length(h5dimnames) != ndim)
            stop(wmsg("'h5dimnames' must be a character vector containing ",
                      "the names of the HDF5 datasets where to write the ",
                      "dimnames of dataset \"", name, "\" (one per dimension ",
                      "in \"", name, "\")"))
        if (any(not_NULL & is.na(h5dimnames)))
            stop(wmsg("'h5dimnames' cannot have NAs associated with ",
                      "list elements in 'dimnames' that are not NULL"))
    }
    if (nzchar(group))
        h5dimnames <- paste0(group, "/", h5dimnames)
    h5dimnames[!not_NULL] <- NA_character_
    for (along in which(not_NULL)) {
        h5dn <- h5dimnames[[along]]
        if (h5exists(filepath, h5dn))
            stop(wmsg("HDF5 dataset \"", h5dn, "\" already exists"))
    }
    h5dimnames
}

### dimnames:   A list (possibly named) with 1 list element per dimension in
###             dataset 'name'.
### name:       The name of the HDF5 dataset on which to set the dimnames.
### group:      The name of the HDF5 group where to write the dimnames.
###             If NA, the group name is automatically generated from 'name'.
###             An empty string ("") means that no group should be used.
###             Otherwise, the names in 'h5dimnames' must be relative to the
###             specified group name.
### h5dimnames: A character vector containing the names of the HDF5 datasets
###             (1 per list element in 'dimnames') where to write the dimnames.
###             Names associated with NULL list elements in 'dimnames' are
###             ignored.
h5writeDimnames <- function(dimnames, filepath, name, group=NA, h5dimnames=NULL)
{
    ## 1. Lots of checks.

    ## Before we start writing to the file we want some guarantees that
    ## the full operation will succeed. The checks we make access the file
    ## in read-only mode.
    .check_filepath_and_name(filepath, name)

    not_NULL <- .check_dimnames(dimnames, filepath, name)

    group <- .normarg_group(group, name)

    h5dimnames <- .normarg_h5dimnames(h5dimnames, group, not_NULL,
                                      filepath, name)

    ## 2. Write to the HDF5 file.

    ## Create group if needed.
    if (!is.na(group) && !h5exists(filepath, group))
        h5createGroup(filepath, group)

    ## Write dimnames.
    for (along in which(not_NULL)) {
        dn <- dimnames[[along]]
        h5dn <- h5dimnames[[along]]
        h5write(dn, filepath, h5dn)
    }

    ## Attach new datasets to dimensions of dataset 'name'.
    set_h5dimnames(filepath, name, h5dimnames)

    ## Set the dimension labels.
    dimlabels <- names(dimnames)
    if (!is.null(dimlabels) && any(nzchar(dimlabels)))
        h5setdimlabels(filepath, name, dimlabels)
}

h5readDimnames <- function(filepath, name, as.character=FALSE)
{
    if (!isTRUEorFALSE(as.character))
        stop(wmsg("'as.character' must be TRUE or FALSE"))
    h5dimnames <- get_h5dimnames(filepath, name)
    dimlabels <- h5getdimlabels(filepath, name)
    if (all(is.na(h5dimnames)) && is.null(dimlabels))
        return(NULL)
    lapply(setNames(h5dimnames, dimlabels),
           function(h5dn) {
               if (is.na(h5dn))
                   return(NULL)
               dn <- h5mread(filepath, h5dn)
               if (as.character) {
                   ## as.character() drops all attributes so no need to
                   ## explicitly drop the "dim" attribute.
                   as.character(dn)
               } else {
                   S4Arrays:::set_dim(dn, NULL)
               }
           })
}

