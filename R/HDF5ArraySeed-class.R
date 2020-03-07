### =========================================================================
### HDF5ArraySeed objects
### -------------------------------------------------------------------------


setClass("HDF5ArraySeed",
    contains="Array",
    representation(
        filepath="character",       # Absolute path to the HDF5 file so the
                                    # object doesn't break when the user
                                    # changes the working directory (e.g.
                                    # with setwd()).
                                    # The path must also be in its canonical
                                    # form so paths from different objects
                                    # can be compared (required by
                                    # quickResaveHDF5SummarizedExperiment()).
        name="character",           # Name of the dataset in the HDF5 file.
        type="character",           # NA or the wanted type.
        dim="integer",
        chunkdim="integer_OR_NULL",
        first_val="ANY"             # First value in the dataset.
    ),
    prototype(
        type=NA_character_
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

### Check that HDF5ArraySeed object 'x' points to an HDF5 dataset that
### is accessible and "as expected".
validate_HDF5ArraySeed_dataset <- function(x)
{
    ## Check that 'x' points to an HDF5 file that is accessible.
    if (!file.exists(x@filepath))
        return(paste0("points to an HDF5 file that does not exist: ",
                      x@filepath))
    if (dir.exists(x@filepath))
        return(paste0("points to a directory ('", x@filepath, "') ",
                      "instead of an HDF5 file"))
    h5_content <- try(h5ls(x@filepath), silent=TRUE)
    if (inherits(h5_content, "try-error"))
        return(paste0("points to an invalid HDF5 file: ", x@filepath))
    if (x@filepath != file_path_as_absolute(x@filepath))
        return(paste0("uses a non-absolute/non-canonical path ",
                      "('", x@filepath, "') to point to the HDF5 file"))

    ## Check that 'x' points to an HDF5 dataset that is accessible.
    h5_dim <- try(h5dim(x@filepath, x@name), silent=TRUE)
    if (inherits(h5_dim, "try-error"))
        return(paste0("points to an HDF5 dataset ('", x@name, "') ",
                      "that does not exist in HDF5 file: ", x@filepath))

    ## Check that 'x' points to an HDF5 dataset that has the
    ## expected dimensions and chunk dimensions.
    if (!identical(h5_dim, x@dim))
        return(paste0("points to an HDF5 dataset ('", x@name, "') ",
                      "in HDF5 file '", x@filepath, "' ",
                      "that does not have the expected dimensions"))
    h5_chunkdim <- h5chunkdim(x@filepath, x@name, adjust=TRUE)
    if (!identical(h5_chunkdim, x@chunkdim))
        return(paste0("points to an HDF5 dataset ('", x@name, "') ",
                      "in HDF5 file '", x@filepath, "' ",
                      "that does not have the expected chunk dimensions"))

    TRUE
}

.validate_HDF5ArraySeed <- function(x)
{
    ## 'filepath' slot.
    x_filepath <- x@filepath
    if (!isSingleString(x_filepath))
        return("'filepath' slot must be a single string")

    ## 'name' slot.
    x_name <- x@name
    if (!isSingleString(x_name))
        return("'name' slot must be a single string")

    ## 'dim' slot.
    msg <- DelayedArray:::validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)

    ## 'chunkdim' slot.
    x_chunkdim <- x@chunkdim
    if (!is.null(x_chunkdim)) {
        msg <- DelayedArray:::validate_dim_slot(x, "chunkdim")
        if (!isTRUE(msg))
            return(msg)
    }

    ## Check that 'x' points to an HDF5 dataset that is accessible
    ## and "as expected".
    msg <- validate_HDF5ArraySeed_dataset(x)
    if (!isTRUE(msg))
        return(paste0("object ", msg))

    ## Check that the dimnames stored in the file are consistent with
    ## the dimensions of the HDF5 dataset.
    msg <- validate_lengths_of_h5dimnames(x_filepath, x_name)
    if (!isTRUE(msg))
        return(msg)

    TRUE
}

setValidity2("HDF5ArraySeed", .validate_HDF5ArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### Does NOT access the file.
setMethod("path", "HDF5ArraySeed", function(object) object@filepath)

normarg_path <- function(path, what1, what2)
{
    if (!isSingleString(path))
        stop(wmsg(what1, " must be a single string specifying the path ",
                  "to the file where the ", what2, " is located"))
    file_path_as_absolute(path)  # return absolute path in canonical form
}

### Return a fake value (of the correct type) if the dataset is empty i.e.
### if at least one of its dimensions is 0.
.read_h5dataset_first_val <- function(filepath, name, dim)
{
    if (any(dim == 0L)) {
        type <- get_h5mread_returned_type(filepath, name)
        val <- vector(type, 1L)  # fake value
    } else {
        index <- rep.int(list(1L), length(dim))
        ans <- h5read2(filepath, name, index)
        stopifnot(length(ans) == 1L)  # sanity check
        val <- ans[[1L]]  # drop any attribute
    }
    val
}

setReplaceMethod("path", "HDF5ArraySeed",
    function(object, value)
    {
        new_filepath <- normarg_path(value, "the supplied path",
                                            "HDF5 dataset")

        ## Check dim compatibility.
        new_dim <- h5dim(new_filepath, object@name)
        object_dim <- object@dim
        if (!identical(new_dim, object_dim)) {
            new_dim_in1string <- paste0(new_dim, collapse=" x ")
            dim_in1string <- paste0(object_dim, collapse=" x ")
            stop(wmsg("dimensions (", new_dim_in1string, ") ",
                      "of HDF5 dataset '", object@name, "' ",
                      "from file '", value, "' are not ",
                      "as expected (", dim_in1string, ")"))
        }

        ## Check first val compatibility.
        new_first_val <- .read_h5dataset_first_val(new_filepath,
                                                   object@name,
                                                   object_dim)
        if (!identical(new_first_val, object@first_val))
            stop(wmsg("first value in HDF5 dataset '", object@name, "' ",
                      "from file '", value, "' is not ",
                      "as expected"))

        ## Set new path.
        object@filepath <- new_filepath
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() getter
###

### Override the default method (defined in the DelayedArray package) with
### a much faster one.
setMethod("type", "HDF5ArraySeed",
    function(x)
    {
        ## Prior to HDF5Array 1.15.6 HDF5ArraySeed objects didn't have
        ## the "type" slot.
        if (!.hasSlot(x, "type"))
            return(type(x@first_val))
        type <- x@type
        if (is.na(type))
            type <- get_h5mread_returned_type(path(x), x@name)
        type
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() getter
###

### Does NOT access the file.
setMethod("dim", "HDF5ArraySeed", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

### If the user requested a specific type when HDF5ArraySeed object 'x'
### was constructed then we must return an array of that type.
.extract_array_from_HDF5ArraySeed <- function(x, index)
{
    ## Prior to HDF5Array 1.15.6 HDF5ArraySeed objects didn't have
    ## the "type" slot.
    if (!.hasSlot(x, "type"))
        return(h5read2(path(x), x@name, index))
    as_int <- !is.na(x@type) && x@type == "integer"
    ans <- h5read2(path(x), x@name, index, as.integer=as_int)
    if (!is.na(x@type) && typeof(ans) != x@type)
        storage.mode(ans) <- x@type
    ans
}

setMethod("extract_array", "HDF5ArraySeed", .extract_array_from_HDF5ArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkdim() getter
###

### Does NOT access the file.
setMethod("chunkdim", "HDF5ArraySeed", function(x) x@chunkdim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Return an HDF5ArraySeed object with NO dimnames!
### FIXME: Investigate the possiblity to store the dimnames in the HDF5 file
### and make dimnames() on the object returned by HDF5ArraySeed() bring them
### back.
HDF5ArraySeed <- function(filepath, name, type=NA)
{
    ## Check 'filepath'.
    filepath <- normarg_path(filepath, "'filepath'", "HDF5 dataset")

    ## Check 'name'.
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying ",
                  "the name of the dataset in the HDF5 file"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))

    ## Check 'type'
    if (!isSingleStringOrNA(type))
        stop(wmsg("'type' must be a single string or NA"))
    if (is.na(type)) {
        type <- as.character(type)
    } else if (type != "list") {
        tmp <- try(vector(type), silent=TRUE)
        if (inherits(tmp, "try-error") || !is.atomic(tmp))
            stop(wmsg("'type' must be an R atomic type ",
                      "(e.g. \"integer\") or \"list\""))
    }

    dim <- h5dim(filepath, name)
    chunkdim <- h5chunkdim(filepath, name, adjust=TRUE)
    first_val <- .read_h5dataset_first_val(filepath, name, dim)

    new2("HDF5ArraySeed", filepath=filepath,
                          name=name,
                          type=type,
                          dim=dim,
                          chunkdim=chunkdim,
                          first_val=first_val)
}

