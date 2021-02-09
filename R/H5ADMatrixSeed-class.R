### =========================================================================
### H5ADMatrixSeed() constructor
### -------------------------------------------------------------------------


setClass("Dense_H5ADMatrixSeed", contains="HDF5ArraySeed")
setClass("Sparse_H5ADMatrixSeed", contains="H5SparseMatrixSeed")

setClassUnion(
    "H5ADMatrixSeed", c("Dense_H5ADMatrixSeed", "Sparse_H5ADMatrixSeed")
)


### Returns either a Dense_H5ADMatrixSeed or Sparse_H5ADMatrixSeed object.
H5ADMatrixSeed <- function(filepath)
{
    if (!isSingleString(filepath)) 
        stop(wmsg("'filepath' must be a single string specifying ", 
                  "the path to the h5ad file"))
    filepath <- file_path_as_absolute(filepath)

    if (!h5exists(filepath, "X"))
        stop(wmsg("File '", filepath, "' contains no dataset or group ",
                  "named 'X'. Is this a valid h5ad file?"))

    h5attrs <- h5readAttributes(filepath, "X")
    if (length(h5attrs) == 0L) {
        ans <- HDF5ArraySeed(filepath, "X")
        if (length(dim(ans)) != 2L)
            stop(wmsg("Dataset 'X' in file '", filepath, "' does not have ",
                      "exactly 2 dimensions. Is this a valid h5ad file?"))
        ans <- new2("Dense_H5ADMatrixSeed", ans)
    } else {
        ans <- H5SparseMatrixSeed(filepath, "X")
        ans <- new2("Sparse_H5ADMatrixSeed", ans)
    }
    ans
}

