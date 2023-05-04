### =========================================================================
### HDF5ArraySeed objects
### -------------------------------------------------------------------------


setClass("HDF5ArraySeed",
    contains="Array",
    representation(
      ## ----------------- user supplied slots -----------------

        ## H5File object or **absolute** path to a local HDF5 file so the
        ## object won't break when the user changes the working directory
        ## (e.g. with setwd()). The path must also be in its canonical
        ## form so comparing paths from different objects is meaningful
        ## (required by quickResaveHDF5SummarizedExperiment()).
        filepath="character_OR_H5File",

        ## Name of dataset in the HDF5 file.
        name="character",

        ## Whether the HDF5 dataset should be considered sparse (and treated
        ## as such) or not. Slot added in HDF5Array 1.17.8.
        as_sparse="logical",  # TRUE or FALSE

        ## NA or the desired type. Slot added in HDF5Array 1.15.6.
        type="character",

      ## ------------ automatically populated slots ------------

        dim="integer",
        chunkdim="integer_OR_NULL",
        first_val="ANY"  # first value in the dataset
    ),
    prototype(
        as_sparse=FALSE,
        type=NA_character_
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

### Check that 'x' points to an HDF5 dataset that has the expected dimensions
### and chunk dimensions.
validate_HDF5ArraySeed_dataset_geometry <- function(x, what="object")
{
    h5_dim <- h5dim(x@filepath, x@name)
    if (!identical(h5_dim, x@dim))
        return(paste0(what, " points to an HDF5 dataset (\"", x@name, "\") ",
                      "in HDF5 file \"", x@filepath, "\" ",
                      "that does not have the expected dimensions"))
    h5_chunkdim <- h5chunkdim(x@filepath, x@name, adjust=TRUE)
    if (!identical(h5_chunkdim, x@chunkdim))
        return(paste0(what, " points to an HDF5 dataset (\"", x@name, "\") ",
                      "in HDF5 file \"", x@filepath, "\" ",
                      "that does not have the expected chunk dimensions"))
    TRUE
}

.validate_HDF5ArraySeed <- function(x)
{
    ## 'filepath' and 'name' slots.
    x_filepath <- x@filepath
    x_name <- x@name
    if (is(x_filepath, "H5File")) {
        ## TODO: Implement the H5File case.
        ## Note that using 'validObject(x@filepath)' won't be enough
        ## because a closed H5File object is considered valid. We want to make
        ## sure that the H5File object is opened and has a working file ID.
    } else {
        msg <- validate_h5_absolute_path(x_filepath, "'filepath' slot")
        if (!isTRUE(msg))
            return(msg)
        msg <- validate_h5_dataset_name(x_filepath, x_name, "'name' slot")
        if (!isTRUE(msg))
            return(msg)
    }

    ## 'as_sparse' slot.
    x_as_sparse <- x@as_sparse
    if (!isTRUEorFALSE(x_as_sparse))
        return("'as_sparse' slot must be TRUE or FALSE")

    ## 'dim' slot.
    msg <- S4Arrays:::validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)

    ## 'chunkdim' slot.
    x_chunkdim <- x@chunkdim
    if (!is.null(x_chunkdim)) {
        msg <- S4Arrays:::validate_dim_slot(x, "chunkdim")
        if (!isTRUE(msg))
            return(msg)
    }

    if (!is(x_filepath, "H5File")) {
        ## Check that the dataset has the expected dimensions and
        ## chunk dimensions.
        msg <- validate_HDF5ArraySeed_dataset_geometry(x)
        if (!isTRUE(msg))
            return(msg)
    }

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
setMethod("path", "HDF5ArraySeed",
    function(object)
    {
        filepath <- object@filepath
        if (is(filepath, "H5File"))
            filepath <- path(filepath)
        filepath
    }
)

### Return a fake value (of the correct type) if the dataset is empty i.e.
### if at least one of its dimensions is 0.
.read_h5dataset_first_val <- function(filepath, name, dim)
{
    if (any(dim == 0L)) {
        type <- get_h5mread_returned_type(filepath, name)
        val <- vector(type, 1L)  # fake value
    } else {
        index <- rep.int(list(1L), length(dim))
        ans <- h5mread(filepath, name, index)
        stopifnot(length(ans) == 1L)  # sanity check
        val <- ans[[1L]]  # drop any attribute
    }
    val
}

setReplaceMethod("path", "HDF5ArraySeed",
    function(object, value)
    {
        if (is(value, "H5File")) {
            new_filepath <- value
            value <- path(value)
            ## Check dim compatibility.
            ## TODO: Implement this.
        } else {
            new_filepath <- normarg_h5_filepath(value,
                                                what1="the supplied path",
                                                what2="the HDF5 dataset")
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
        }
        ## Check first val compatibility.
        new_first_val <- .read_h5dataset_first_val(new_filepath,
                                                   object@name,
                                                   object_dim)
        if (!identical(new_first_val, object@first_val))
            stop(wmsg("first value in HDF5 dataset '", object@name, "' ",
                      "from file '", value, "' is not as expected"))

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
            type <- get_h5mread_returned_type(x@filepath, x@name)
        type
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() getter
###

### Does NOT access the file.
setMethod("dim", "HDF5ArraySeed", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dimnames() getter
###

### Does access the file!
setMethod("dimnames", "HDF5ArraySeed",
    function(x) h5readDimnames(x@filepath, x@name, as.character=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

### A thin wrapper around h5mread().
### TODO: Maybe we no longer need this. I mean, this is used in the
### extract_array() and OLD_extract_sparse_array() methods for HDF5ArraySeed
### objects but the 'index' passed to these methods should never contain
### RangeNSBS objects. So it's probably ok to get rid of this and to just
### use h5mread() instead.
.h5mread2 <- function(filepath, name, index=NULL,
                      as.integer=FALSE, as.sparse=FALSE)
{
    if (!is.null(index))
        index <- S4Arrays:::expand_Nindex_RangeNSBS(index)
    h5mread(filepath, name, starts=index,
            as.integer=as.integer, as.sparse=as.sparse)
}


.extract_array_from_HDF5ArraySeed <- function(x, index)
{
    ## Prior to HDF5Array 1.15.6 HDF5ArraySeed objects didn't have
    ## the "type" slot.
    if (!.hasSlot(x, "type"))
        return(.h5mread2(x@filepath, x@name, index))
    ## If the user requested a specific type when HDF5ArraySeed object 'x'
    ## was constructed then we must return an array of that type.
    as_int <- !is.na(x@type) && x@type == "integer"
    ans <- .h5mread2(x@filepath, x@name, index, as.integer=as_int)
    if (!is.na(x@type) && typeof(ans) != x@type)
        storage.mode(ans) <- x@type
    ans
}

setMethod("extract_array", "HDF5ArraySeed", .extract_array_from_HDF5ArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse() and OLD_extract_sparse_array()
###

### Prior to HDF5Array 1.17.8 HDF5ArraySeed objects didn't have the
### "as_sparse" slot.
setMethod("is_sparse", "HDF5ArraySeed",
    function(x) .hasSlot(x, "as_sparse") && x@as_sparse
)

setReplaceMethod("is_sparse", "HDF5ArraySeed",
    function(x, value)
    {
        if (!isTRUEorFALSE(value))
            stop(wmsg("the supplied value must be TRUE or FALSE"))
        if (!.hasSlot(x, "as_sparse"))
            x <- updateObject(x, check=FALSE)
        x@as_sparse <- value
        x
    }
)

.OLD_extract_sparse_array_from_HDF5ArraySeed <- function(x, index)
{
    if (!is_sparse(x))
        stop(wmsg("calling OLD_extract_sparse_array() on an HDF5ArraySeed ",
                  "object is supported only if the object is sparse"))
    ## Prior to HDF5Array 1.15.6 HDF5ArraySeed objects didn't have
    ## the "type" slot.
    if (!.hasSlot(x, "type"))
        return(.h5mread2(x@filepath, x@name, index, as.sparse=TRUE))
    ## If the user requested a specific type when HDF5ArraySeed object 'x'
    ## was constructed then we must return a SparseArraySeed object of
    ## that type.
    as_int <- !is.na(x@type) && x@type == "integer"
    ans <- .h5mread2(x@filepath, x@name, index, as.integer=as_int,
                                         as.sparse=TRUE)
    if (!is.na(x@type) && typeof(ans) != x@type)
        storage.mode(ans@nzdata) <- x@type
    ans
}

setMethod("OLD_extract_sparse_array", "HDF5ArraySeed",
    .OLD_extract_sparse_array_from_HDF5ArraySeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkdim() getter
###

### Does NOT access the file.
setMethod("chunkdim", "HDF5ArraySeed", function(x) x@chunkdim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

HDF5ArraySeed <- function(filepath, name, as.sparse=FALSE, type=NA)
{
    if (!is(filepath, "H5File"))
        filepath <- normarg_h5_filepath(filepath)
    name <- normarg_h5_name(name)

    ## Check 'as.sparse'.
    if (!isTRUEorFALSE(as.sparse))
        stop(wmsg("'as.sparse' must be TRUE or FALSE"))

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
                          as_sparse=as.sparse,
                          type=type,
                          dim=dim,
                          chunkdim=chunkdim,
                          first_val=first_val)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###

setMethod("updateObject", "HDF5ArraySeed",
    function(object, ..., verbose=FALSE)
    {
        ## The "file" slot was renamed "filepath" in HDF5Array 1.7.3 (commit
        ## b30f4d4b).
        if (!.hasSlot(object, "filepath")) {
            return(new2("HDF5ArraySeed", filepath=object@file,
                                         name=object@name,
                                         type=type(object@first_val),
                                         dim=object@dim,
                                         first_val=object@first_val,
                                         check=FALSE))
        }
        ## The "chunkdim" slot was added in HDF5Array 1.7.7 (commit ef4c5b47).
        if (!.hasSlot(object, "chunkdim")) {
            return(new2("HDF5ArraySeed", filepath=object@filepath,
                                         name=object@name,
                                         type=type(object@first_val),
                                         dim=object@dim,
                                         first_val=object@first_val,
                                         check=FALSE))
        }
        ## The "type" slot was added in HDF5Array 1.15.6.
        if (!.hasSlot(object, "type")) {
            return(new2("HDF5ArraySeed", filepath=object@filepath,
                                         name=object@name,
                                         type=type(object@first_val),
                                         dim=object@dim,
                                         chunkdim=object@chunkdim,
                                         first_val=object@first_val,
                                         check=FALSE))
        }
        ## The "as_sparse" slot was added in HDF5Array 1.17.8.
        if (!.hasSlot(object, "as_sparse")) {
            return(new2("HDF5ArraySeed", filepath=object@filepath,
                                         name=object@name,
                                         type=object@type,
                                         dim=object@dim,
                                         chunkdim=object@chunkdim,
                                         first_val=object@first_val,
                                         check=FALSE))
        }
        object
    }
)

