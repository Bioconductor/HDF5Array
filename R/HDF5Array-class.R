### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


### NOT exported.
setClass("HDF5Dataset",
    representation(
        file="character",   # Single string.
        group="character",  # Single string.
        name="character",   # Dataset name.
        dim="integer",
        first_val="ANY"     # First value in the dataset.
    )
)

setMethod("dim", "HDF5Dataset", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A convenience wrapper to rhdf5::h5read()
###

.quiet_h5read <- function(file, group, name, index)
{
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    suppressWarnings(
        h5read(file, paste(group, name, sep="/"), index=index)
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
            ans <- .quiet_h5read(seed@file, seed@group, seed@name, index)
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Construct an HDF5Dataset object from a file or an array
###

.get_h5dataset_dim <- function(file, group, name)
{
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
.read_h5dataset_first_val <- function(file, group, name, ndim)
{
    index <- rep.int(list(1L), ndim)
    ans <- .quiet_h5read(file, group, name, index)
    stopifnot(length(ans) == 1L)  # sanity check
    ans[[1L]]  # drop any attribute
}

.new_HDF5Dataset_from_file <- function(file, group, name, type=NA)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the HDF5 file where the dataset is located"))
    if (!isSingleString(group))
        stop(wmsg("'group' must be a single string specifying the HDF5 group ",
                  "of the dataset, as reported by rhdf5::h5ls()"))
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the dataset, as reported by rhdf5::h5ls()"))
    if (!isSingleStringOrNA(type))
        stop("'type' must be a single string or NA")
    dim <- .get_h5dataset_dim(file, group, name)
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
        first_val <- .read_h5dataset_first_val(file, group, name, length(dim))
        detected_type <- typeof(first_val)
        if (!(is.na(type) || type == detected_type))
            warning(wmsg("The type specified via the 'type' argument (",
                         type, ") doesn't match the type of this HDF5 ",
                         "dataset (", detected_type, "). Ignoring the ",
                         "former."))
    }
    new2("HDF5Dataset", file=file,
                        group=group,
                        name=name,
                        dim=dim,
                        first_val=first_val)
}

.new_HDF5Dataset_from_array <- function(a)
{
    out_file <- getHDF5ArrayOutputFile()
    out_name <- getHDF5ArrayOutputName()
    on.exit(setHDF5ArrayOutputName())

    h5write(a, out_file, out_name)
    .new_HDF5Dataset_from_file(out_file, "/", out_name, type=type(a))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5Array and HDF5Matrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate HDF5Array and HDF5Matrix objects instead of DelayedArray
### and DelayedMatrix objects.

setClass("HDF5Array", contains="DelayedArray")

setAs("HDF5Dataset", "HDF5Array",
    function(from) new_DelayedArray(from, Class="HDF5Array")
)

HDF5Array <- function(file, group, name, type=NA)
{
    as(.new_HDF5Dataset_from_file(file, group, name, type=type), "HDF5Array")
}

setAs("array", "HDF5Array",
    function(from)
    {
        ans <- as(.new_HDF5Dataset_from_array(from), "HDF5Array")
        ## TODO: Investigate the possiblity that .new_HDF5Dataset_from_array()
        ## stores the dimnames in the HDF5 file so .new_HDF5Dataset_from_file()
        ## can bring them back. Then we wouldn't need to explicitely set them
        ## on 'ans' like we do here.
        dimnames(ans) <- dimnames(from)
        ans
    }
)

setClass("HDF5Matrix", contains=c("DelayedMatrix", "HDF5Array"))

setAs("HDF5Array", "HDF5Matrix",
    function(from) as(as(from, "DelayedMatrix"), "HDF5Matrix")
)

HDF5Matrix <- function(file, group, name, type=NA)
{
    as(HDF5Array(file, group, name, type=type), "HDF5Matrix")
}

setAs("matrix", "HDF5Matrix",
    function(from) as(as(from, "HDF5Array"), "HDF5Matrix")
)

setAs("DelayedArray", "HDF5Array",
    function(from) stop(wmsg("coercing a ", class(from), " object to an ",
                             "HDF5Array object is not supported yet"))
)

setAs("DelayedArray", "HDF5Matrix",
    function(from) stop(wmsg("Coercing a ", class(from), " object to an ",
                             "HDF5Matrix object is not supported yet. ",
                             "Please coerce to DelayedMatrix instead."))
)

