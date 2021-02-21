### =========================================================================
### H5ADMatrixSeed objects
### -------------------------------------------------------------------------


setClass("H5ADMatrixSeed", contains="Array", representation("VIRTUAL"))

setClass("Dense_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "HDF5ArraySeed")
)
setClass("CSC_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "CSC_H5SparseMatrixSeed")
)
setClass("CSR_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "CSR_H5SparseMatrixSeed")
)

### Returns an H5ADMatrixSeed derivative (can be either a Dense_H5ADMatrixSeed,
### or a CSC_H5SparseMatrixSeed, or a CSR_H5SparseMatrixSeed object).
H5ADMatrixSeed <- function(filepath, name="X")
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying the ",
                  "path to the h5ad file where the matrix is located"))
    filepath <- file_path_as_absolute(filepath)

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
        ans <- new2("Dense_H5ADMatrixSeed", ans0)
    } else if (h5isgroup(filepath, name)) {
        ans0 <- H5SparseMatrixSeed(filepath, name)
        if (is(ans0, "CSC_H5SparseMatrixSeed"))
            ans_class <- "CSC_H5ADMatrixSeed"
        else
            ans_class <- "CSR_H5ADMatrixSeed"
        ans <- new2(ans_class, ans0)
    } else {
        stop(wmsg("file '", filepath, "' contains no dataset or group ",
                  "named '", name, "'"))
    }
    ans
}

