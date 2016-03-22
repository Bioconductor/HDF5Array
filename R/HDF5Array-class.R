### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


### NOT exported.
setClass("HDF5Dataset",
    representation(
        file="character",   # Single string.
        name="character",   # Dataset name.
        dim="integer",
        first_val="ANY"     # First value in the dataset.
    )
)

setMethod("dim", "HDF5Dataset", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A convenience wrapper to rhdf5::h5read()
###

.quiet_h5read <- function(file, name, index)
{
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    suppressWarnings(
        h5read(file, name, index=index)
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array_from_seed()
###

setMethod("extract_array_from_seed", "HDF5Dataset",
    function(seed, index)
    {
        if (any(lengths(index) == 0L)) {
            ans <- seed@first_val[0]
        } else {
            ans <- .quiet_h5read(seed@file, seed@name, index)
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Construct an HDF5Dataset object from a file or an array
###

.get_h5dataset_dim <- function(file, name)
{
    group <- gsub("(.*/)[^/]*$", "\\1", name)
    name <- gsub(".*/([^/]*)$", "\\1", name)
    f <- H5Fopen(file, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(f))
    g <- H5Gopen(f, group)
    on.exit(H5Gclose(g), add=TRUE)
    d <- H5Dopen(g, name)
    on.exit(H5Dclose(d), add=TRUE)
    H5Sget_simple_extent_dims(H5Dget_space(d))$size
}

### Will fail if the dataset is empty (i.e. if at least one of its dimensions
### is 0).
.read_h5dataset_first_val <- function(file, name, ndim)
{
    index <- rep.int(list(1L), ndim)
    ans <- .quiet_h5read(file, name, index)
    stopifnot(length(ans) == 1L)  # sanity check
    ans[[1L]]  # drop any attribute
}

.new_HDF5Dataset_from_file <- function(file, name, type=NA)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the HDF5 file where the dataset is located"))
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the dataset in the HDF5 file"))
    if (!isSingleStringOrNA(type))
        stop("'type' must be a single string or NA")
    dim <- .get_h5dataset_dim(file, name)
    if (any(dim == 0L)) {
        if (is.na(type))
            stop(wmsg("This HDF5 dataset is empty! Don't know how to ",
                      "determine the type of an empty HDF5 dataset at the ",
                      "moment. Please use the 'type' argument to help me ",
                      "(see '?HDF5Array' for more information)."))
        first_val <- match.fun(type)(1)  # fake value
        if (!is.atomic(first_val))
            stop("invalid type: ", type)
    } else {
        first_val <- .read_h5dataset_first_val(file, name, length(dim))
        detected_type <- typeof(first_val)
        if (!(is.na(type) || type == detected_type))
            warning(wmsg("The type specified via the 'type' argument (",
                         type, ") doesn't match the type of this HDF5 ",
                         "dataset (", detected_type, "). Ignoring the ",
                         "former."))
    }
    new2("HDF5Dataset", file=file,
                        name=name,
                        dim=dim,
                        first_val=first_val)
}

### TODO: Investigate the possiblity to store the dimnames in the HDF5 file
### so a "dimnames" method for HDF5Dataset objects could bring them back.
### Then this coercion would propagate the dimnames.
.new_HDF5Dataset_from_array <- function(a)
{
    if (!is.array(a))
        stop("cannot create an HDF5Dataset object from a ", class(a))
    out_file <- getHDF5OutputFile()
    out_name <- getHDF5OutputName()
    on.exit(setHDF5OutputName())

    h5write(a, out_file, out_name)
    .new_HDF5Dataset_from_file(out_file, out_name, type=type(a))
}

HDF5Dataset <- function(file, name, type=NA)
{
    if (!missing(name))
        return(.new_HDF5Dataset_from_file(file, name, type=type))
    if (!identical(type, NA))
        warning("ignoring supplied 'type'")
    .new_HDF5Dataset_from_array(file)
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

### Overwrite unsafe automatic coercion method that return invalid objects (it
### doesn't validate them).
setAs("HDF5Array", "HDF5Matrix", function(from) new("HDF5Matrix", from))

### For internal use only.
setMethod("matrixClass", "HDF5Array", function(x) "HDF5Matrix")

.validate_HDF5Array <- function(x)
{
    if (!is_pristine(x))
        return(wmsg("'x' carries delayed operations on it"))
    if (!is(x@seeds[[1L]], "HDF5Dataset"))
        return(wmsg("'x@seeds' must be a HDF5Dataset object"))
    TRUE
}

setValidity2("HDF5Array", .validate_HDF5Array)

HDF5Array <- function(file, name, type=NA)
{
    if (is(file, "HDF5Dataset")) {
        seed <- file
    } else {
        seed <- HDF5Dataset(file, name, type=type)
    }
    ans <- new_DelayedArray(seed, Class="HDF5Array")
    ## The dimnames will automatically propagate once we store them in the
    ## HDF5 file and we have a dimnames() getter for HDF5Dataset objects that
    ## knows how to extract them. And so this won't be necessary anymore...
    if (missing(name))
        dimnames(ans) <- dimnames(file)
    ans
}

setAs("DelayedArray", "HDF5Array",
    function(from) stop(wmsg("coercing a ", class(from), " object to an ",
                             "HDF5Array object is not supported yet"))
)

setAs("DelayedArray", "HDF5Matrix",
    function(from) stop(wmsg("Coercing a ", class(from), " object to an ",
                             "HDF5Matrix object is not supported yet. ",
                             "Please coerce to DelayedMatrix instead."))
)

setAs("DelayedMatrix", "HDF5Matrix",
    function(from) stop(wmsg("coercing a ", class(from), " object to an ",
                             "HDF5Matrix object is not supported yet"))
)

