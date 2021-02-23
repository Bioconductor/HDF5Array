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
### dimnames() method for Dense_H5ADMatrixSeed objects (overwrite method for
### HDF5Array objects)
###

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

### Must return a list of length 2.
.load_h5ad_dimnames <- function(filepath)
{
    ROWNAMES_DATASET <- "/var/_index"
    COLNAMES_DATASET <- "/obs/_index"

    ok1 <- try(h5isdataset(filepath, ROWNAMES_DATASET), silent=TRUE)
    ok2 <- try(h5isdataset(filepath, COLNAMES_DATASET), silent=TRUE)
    if (inherits(ok1, "try-error") || inherits(ok2, "try-error"))
        stop(wmsg("Cannot read ", filepath, ". Is this a valid HDF5 file?"))

    ans_rownames <- ans_colnames <- NULL
    if (ok1 || ok2) {
        if (ok1) {
            ans_rownames <- h5mread(filepath, ROWNAMES_DATASET)
        } else {
            warning(wmsg("could not find rownames (normally stored ",
                         "in dataset '", ROWNAMES_DATASET, "' ",
                         "in this .h5ad file"))
        }
        if (ok2) {
            ans_colnames <- h5mread(filepath, COLNAMES_DATASET)
        } else {
            warning(wmsg("could not find colnames (normally stored ",
                         "in dataset '", COLNAMES_DATASET, "' ",
                         "in this .h5ad file"))
        }
    } else {
        warning(wmsg("could not find dimnames (normally stored ",
                     "in datasets '", ROWNAMES_DATASET, "' ",
                     "and '", COLNAMES_DATASET, "') in this .h5ad file"))
    }
    list(ans_rownames, ans_colnames)
}

### Returns an H5ADMatrixSeed derivative (can be either a Dense_H5ADMatrixSeed,
### or a CSC_H5SparseMatrixSeed, or a CSR_H5SparseMatrixSeed object).
H5ADMatrixSeed <- function(filepath, name="X")
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying the ",
                  "path to the h5ad file where the matrix is located"))
    filepath <- file_path_as_absolute(filepath)
    dimnames <- .load_h5ad_dimnames(filepath)

    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name of ",
                  "the dataset or group in the h5ad file that stores the ",
                  "matrix"))
    if (name == "")
        stop(wmsg("'name' cannot be the empty string"))

    if (h5isdataset(filepath, name)) {
        ans0 <- HDF5ArraySeed(filepath, name)
        if (length(dim(ans0)) != 2L)
            stop(wmsg("Dataset '", name, "' in file '", filepath, "' does ",
                      "not have exactly 2 dimensions."))
        ans <- new2("Dense_H5ADMatrixSeed", ans0, dimnames=dimnames)
    } else if (h5isgroup(filepath, name)) {
        ans0 <- H5SparseMatrixSeed(filepath, name)
        if (is(ans0, "CSC_H5SparseMatrixSeed"))
            ans_class <- "CSC_H5ADMatrixSeed"
        else
            ans_class <- "CSR_H5ADMatrixSeed"
        ans <- new2(ans_class, ans0, dimnames=dimnames)
    } else {
        stop(wmsg("file '", filepath, "' contains no dataset or group ",
                  "named '", name, "'"))
    }
    ans
}

