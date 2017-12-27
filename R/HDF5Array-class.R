### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


setClass("HDF5ArraySeed",
    contains="Array",
    representation(
        filepath="character",  # Absolute path to the HDF5 file so the object
                               # doesn't break when the user changes the
                               # working directory (e.g. with setwd()).
        name="character",      # Name of the dataset in the HDF5 file.
        dim="integer",
        first_val="ANY"        # First value in the dataset.
    )
)

setMethod("dim", "HDF5ArraySeed", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.extract_array_from_HDF5ArraySeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)) {
        ans <- x@first_val[0]
        dim(ans) <- ans_dim
    } else {
        ans <- h5read2(x@filepath, x@name, index)
    }
    ans
}

setMethod("extract_array", "HDF5ArraySeed", .extract_array_from_HDF5ArraySeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5ArraySeed constructor
###

.get_h5dataset_dim <- function(filepath, name)
{
    dim <- h5dim(filepath, name)
    if (!is.integer(dim)) {
        if (any(dim > .Machine$integer.max)) {
            dim_in1string <- paste0(dim, collapse=" x ")
            stop(wmsg("The dimensions of HDF5 dataset '", name, "' are: ",
                      dim_in1string, "\n\nThe HDF5Array package only ",
                      "supports datasets with all dimensions <= 2^31-1",
                      " (this is ", .Machine$integer.max, ") at the moment."))
        }
        dim <- as.integer(dim)
    }
    dim
}

### Will fail if the dataset is empty (i.e. if at least one of its
### dimensions is 0).
.read_h5dataset_first_val <- function(filepath, name, ndim)
{
    index <- rep.int(list(1L), ndim)
    ans <- h5read2(filepath, name, index)
    stopifnot(length(ans) == 1L)  # sanity check
    ans[[1L]]  # drop any attribute
}

### Return a HDF5ArraySeed object with NO dimnames!
### FIXME: Investigate the possiblity to store the dimnames in the HDF5 file
### and make dimnames() on the object returned by HDF5ArraySeed() bring them
### back.
HDF5ArraySeed <- function(filepath, name, type=NA)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying the path ",
                  "to the HDF5 file where the dataset is located"))
    filepath <- file_path_as_absolute(filepath)
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the dataset in the HDF5 file"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
    if (!isSingleStringOrNA(type))
        stop("'type' must be a single string or NA")
    dim <- .get_h5dataset_dim(filepath, name)
    if (any(dim == 0L)) {
        if (is.na(type))
            stop(wmsg("This HDF5 dataset is empty! Don't know how to ",
                      "determine the type of an empty HDF5 dataset at the ",
                      "moment. Please use the 'type' argument to help me ",
                      "(see '?HDF5Array' for more information)."))
        first_val <- match.fun(type)(1)  # fake value
        if (!is.atomic(first_val))
            stop(wmsg("invalid type: ", type))
    } else {
        first_val <- .read_h5dataset_first_val(filepath, name, length(dim))
        detected_type <- typeof(first_val)
        if (!(is.na(type) || type == detected_type))
            warning(wmsg("The type specified via the 'type' argument (",
                         type, ") doesn't match the type of this HDF5 ",
                         "dataset (", detected_type, "). Ignoring the ",
                         "former."))
    }
    new2("HDF5ArraySeed", filepath=filepath,
                          name=name,
                          dim=dim,
                          first_val=first_val)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5Array and HDF5Matrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate HDF5Array and HDF5Matrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

setClass("HDF5Array", contains="DelayedArray")

setClass("HDF5Matrix", contains=c("DelayedMatrix", "HDF5Array"))

### Automatic coercion method from HDF5Array to HDF5Matrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("HDF5Array", "HDF5Matrix", function(from) new("HDF5Matrix", from))

### For internal use only.
setMethod("matrixClass", "HDF5Array", function(x) "HDF5Matrix")

.validate_HDF5Array <- function(x)
{
    if (!is(x@seed, "HDF5ArraySeed"))
        return(wmsg("'x@seed' must be a HDF5ArraySeed object"))
    if (!DelayedArray:::is_pristine(x))
        return(wmsg("'x' carries delayed operations"))
    TRUE
}

setValidity2("HDF5Array", .validate_HDF5Array)

setAs("ANY", "HDF5Matrix",
    function(from) as(as(from, "HDF5Array"), "HDF5Matrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "HDF5ArraySeed",
    function(seed) DelayedArray:::new_DelayedArray(seed, Class="HDF5Array")
)

### Works directly on a HDF5ArraySeed object, in which case it must be called
### with a single argument.
HDF5Array <- function(filepath, name, type=NA)
{
    if (is(filepath, "HDF5ArraySeed")) {
        if (!(missing(name) && identical(type, NA)))
            stop(wmsg("HDF5Array() must be called with a single argument ",
                      "when passed a HDF5ArraySeed object"))
        seed <- filepath
    } else {
        seed <- HDF5ArraySeed(filepath, name, type=type)
    }
    DelayedArray(seed)
}

