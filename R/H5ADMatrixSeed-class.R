### =========================================================================
### H5ADMatrixSeed objects
### -------------------------------------------------------------------------


setClass("H5ADMatrixSeed", contains="Array", representation("VIRTUAL"))

setClass("Dense_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "HDF5ArraySeed"),
    representation(dimnames="list"),
    prototype(dimnames=list(NULL, NULL))
)
setClass("CSC_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "CSC_H5SparseMatrixSeed")
)
setClass("CSR_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "CSR_H5SparseMatrixSeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dimnames() method for Dense_H5ADMatrixSeed objects
###

### We overwrite the method for HDF5ArraySeed objects with a method that
### accesses the slot, not the file.
setMethod("dimnames", "Dense_H5ADMatrixSeed",
    function(x) DelayedArray:::simplify_NULL_dimnames(x@dimnames)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transposition
###

### S3/S4 combo for t.CSC_H5ADMatrixSeed
t.CSC_H5ADMatrixSeed <- function(x)
{
    new2("CSR_H5ADMatrixSeed", callNextMethod())
}
setMethod("t", "CSC_H5ADMatrixSeed", t.CSC_H5ADMatrixSeed)

### S3/S4 combo for t.CSR_H5ADMatrixSeed
t.CSR_H5ADMatrixSeed <- function(x)
{
    new2("CSC_H5ADMatrixSeed", callNextMethod())
}
setMethod("t", "CSR_H5ADMatrixSeed", t.CSR_H5ADMatrixSeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.load_h5ad_rownames <- function(filepath, name="var")
{
    ok <- try(h5isdataset(filepath, name), silent=TRUE)
    if (isTRUE(ok)) {
        ## Must use rhdf5::h5read() for now, until h5mread() knows how
        ## to read COMPOUND datasets.
        ans <- h5read(filepath, name)$index
        if (!is.null(ans))
            ans <- as.character(ans)
        return(ans)
    }
    ok <- try(h5isgroup(filepath, name), silent=TRUE)
    if (!isTRUE(ok))
        return(NULL)
    ROWNAMES_DATASET <- paste0(name, "/_index")
    ok <- try(h5isdataset(filepath, ROWNAMES_DATASET), silent=TRUE)
    if (!isTRUE(ok))
        return(NULL)
    as.character(h5mread(filepath, ROWNAMES_DATASET))
}

### Must return a list of length 2.
.load_h5ad_dimnames <- function(filepath)
{
    ans_rownames <- .load_h5ad_rownames(filepath)
    ans_colnames <- .load_h5ad_rownames(filepath, name="obs")
    if (is.null(ans_rownames) && is.null(ans_colnames))
        warning(wmsg("could not find dimnames in this h5ad file"))
    list(ans_rownames, ans_colnames)
}

### Returns an H5ADMatrixSeed derivative (can be either a Dense_H5ADMatrixSeed,
### or a CSC_H5SparseMatrixSeed, or a CSR_H5SparseMatrixSeed object).
H5ADMatrixSeed <- function(filepath, layer=NULL)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying the ",
                  "path to the h5ad file"))
    filepath <- file_path_as_absolute(filepath)
    if (is.null(layer)) {
        name <- "/X"
    } else {
        if (!isSingleString(layer) || layer == "")
            stop(wmsg("'layer' must be a single non-empty string"))
        name <- paste0("/layers/", layer)
    }
    if (!h5exists(filepath, name))
        stop(wmsg("HDF5 object \"", name, "\" does not exist ",
                  "in this HDF5 file. Is this a valid h5ad file?"))
    dimnames <- .load_h5ad_dimnames(filepath)

    if (h5isdataset(filepath, name)) {
        ans0 <- HDF5ArraySeed(filepath, name)
        if (length(dim(ans0)) != 2L)
            stop(wmsg("HDF5 dataset \"", name, "\" in file \"", filepath, "\" ",
                      "does not have exactly 2 dimensions. Please consider ",
                      "using the HDF5Array() constructor to access this ",
                      "dataset."))
        ans <- new2("Dense_H5ADMatrixSeed", ans0, dimnames=dimnames)
    } else if (h5isgroup(filepath, name)) {
        ans0 <- H5SparseMatrixSeed(filepath, name)
        if (is(ans0, "CSC_H5SparseMatrixSeed"))
            ans_class <- "CSC_H5ADMatrixSeed"
        else
            ans_class <- "CSR_H5ADMatrixSeed"
        ans <- new2(ans_class, ans0, dimnames=dimnames)
    } else {
        stop(wmsg("HDF5 object \"", name, "\" in file \"", filepath, "\" ",
                  "is neither a dataset or a group. Is this a valid ",
                  "h5ad file?"))
    }
    ans
}

