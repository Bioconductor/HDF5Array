### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


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
### subset_seed_as_array()
###

.subset_HDF5Dataset_as_array <- function(seed, index)
{
    ans_dim <- get_subscripts_lengths(index, dim(seed))
    if (any(ans_dim == 0L)) {
        ans <- seed@first_val[0]
        dim(ans) <- ans_dim
    } else {
        ans <- h5read2(seed@file, seed@name, index)
    }
    ans
}

setMethod("subset_seed_as_array", "HDF5Dataset", .subset_HDF5Dataset_as_array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5Dataset internal low-level constructor
###

.get_h5dataset_dim <- function(file, name)
{
    if (substr(name, 1L, 1L) != "/")
        name <- paste0("/", name)
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
    ans <- h5read2(file, name, index)
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
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeHDF5Dataset()
###

### Return an HDF5Dataset object pointing to the newly written HDF5 dataset
### on disk.
### TODO: Investigate the possiblity to store the dimnames in the HDF5 file
### so a "dimnames" method for HDF5Dataset objects could bring them back.
writeHDF5Dataset <- function(x, file, name)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the HDF5 file where to write the dataset"))
    if (!isSingleString(name)) 
        stop(wmsg("'name' must be a single string specifying the name of ", 
                  "the HDF5 dataset to write"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))
    if (is(x, "HDF5Dataset")) {
        x <- DelayedArray(x)
    } else if (!(is.array(x) || is(x, "DelayedArray"))) {
        stop(wmsg("writing ", class(x), " object as an HDF5 dataset ",
                  "is not supported"))
    }
    if (!file.exists(file))
        h5createFile(file)
    ans_type <- type(x)
    h5createDataset2(file, name, dim(x), storage.mode=ans_type)
    if (is.array(x)) {
        h5write2(x, file, name)
    } else if (is(x, "DelayedArray")) {
        ## Semantically equivalent to 'h5write2(as.array(x), file, name)'
        ## but uses block-processing so the full DelayedArray object is not
        ## realized at once in memory. Instead the object is first split into
        ## blocks and the blocks are realized and written to disk one at a
        ## time.
        block_APPLY(x, identity, out_file=file, out_name=name)
    } else {
        stop(wmsg("writing ", class(x), " object to HDF5 dataset ",
                  "not supported"))
    }
    invisible(.new_HDF5Dataset_from_file(file, name, type=ans_type))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5Dataset constructor
###

HDF5Dataset <- function(file, name, type=NA)
{
    if (!is.array(file)) {
        if (isSingleString(file))
            return(.new_HDF5Dataset_from_file(file, name, type=type))
        if (!is(file, "DelayedArray"))
            stop(wmsg("cannot create an HDF5Dataset object from this object"))
    }
    if (!(missing(name) && identical(type, NA)))
        stop(wmsg("'name' or 'type' cannot be specified when calling ",
                  "HDF5Dataset() on a DelayedArray or array object"))
    out_file <- getHDF5DumpFile()
    out_name <- getHDF5DumpName()
    on.exit(setHDF5DumpName())
    writeHDF5Dataset(file, out_file, out_name)
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
    if (!is(x@seed, "HDF5Dataset"))
        return(wmsg("'x@seed' must be a HDF5Dataset object"))
    if (!is_pristine(x))
        return(wmsg("'x' carries delayed operations on it"))
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

### Overwrite unsafe automatic coercion methods that return invalid objects
### (they don't validate them).
.from_DelayedArray_to_HDF5Array <- function(from)
{
    use_HDF5Dataset_msg <- wmsg(
        "Coercing a ", class(from), " object to HDF5Array or HDF5Matrix is ",
        "not supported. If you intend to realize the object on disk, you ",
        "should instead call the HDF5Dataset() constructor on it."
    )
    stop(use_HDF5Dataset_msg)
}

setAs("DelayedArray", "HDF5Array", .from_DelayedArray_to_HDF5Array)
setAs("DelayedArray", "HDF5Matrix", .from_DelayedArray_to_HDF5Array)
setAs("DelayedMatrix", "HDF5Matrix", .from_DelayedArray_to_HDF5Array)

