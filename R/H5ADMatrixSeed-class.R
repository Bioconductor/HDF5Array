### =========================================================================
### H5ADMatrixSeed objects
### -------------------------------------------------------------------------


setClass("H5ADMatrixSeed", contains="Array", representation("VIRTUAL"))

setClass("Dense_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "HDF5ArraySeed")
)
setClass("Sparse_H5ADMatrixSeed",
    contains=c("H5ADMatrixSeed", "H5SparseMatrixSeed")
)

### Returns either a Dense_H5ADMatrixSeed or Sparse_H5ADMatrixSeed object.
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
        ans <- HDF5ArraySeed(filepath, name)
        if (length(dim(ans)) != 2L)
            stop(wmsg("Dataset '", name, "' in file '", filepath, "' does ",
                      "not have exactly 2 dimensions."))
        ans <- new2("Dense_H5ADMatrixSeed", ans)
    } else if (h5isgroup(filepath, name)) {
        ans <- H5SparseMatrixSeed(filepath, name)
        ans <- new2("Sparse_H5ADMatrixSeed", ans)
    } else {
        stop(wmsg("file '", filepath, "' contains no dataset or group ",
                  "named '", name, "'"))
    }
    ans
}

