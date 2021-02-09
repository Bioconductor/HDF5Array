### =========================================================================
### H5ADMatrix objects
### -------------------------------------------------------------------------
###


setClass("H5ADMatrix",
    contains="DelayedMatrix",
    representation(seed="H5ADMatrixSeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### We must define the 2 methods here. If we define a single method
### for the union class (H5ADMatrixSeed), calling DelayedArray() on
### a Dense_H5ADMatrixSeed object will return an HDF5Matrix object.
### That's because method dispatch will choose the DelayedArray() method
### for HDF5MatrixSeed objects over the method for H5ADMatrixSeed objects.
setMethod("DelayedArray", "Dense_H5ADMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="H5ADMatrix")
)
setMethod("DelayedArray", "Sparse_H5ADMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="H5ADMatrix")
)

### Works directly on a Dense_H5ADMatrixSeed or Sparse_H5ADMatrixSeed object,
### in which case it must be called with a single argument.
H5ADMatrix <- function(filepath)
{
    if (is(filepath, "H5ADMatrixSeed")) {
        #if (!missing(group))
        #    stop(wmsg("H5ADMatrix() must be called with a single argument ",
        #              "when passed a Dense_H5ADMatrixSeed or ",
        #              "Sparse_H5ADMatrixSeed object"))
        seed <- filepath
    } else {
        seed <- H5ADMatrixSeed(filepath)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Taking advantage of sparsity
###

### Will work only if the seed is a Sparse_H5ADMatrixSeed object.
setMethod("sparsity", "H5ADMatrix", function(x) sparsity(x@seed))

setMethod("read_sparse_block", "H5ADMatrix",
    function(x, viewport) read_sparse_block(x@seed, viewport)
)

### Will work only if the seed is a Sparse_H5ADMatrixSeed object.
setMethod("extractNonzeroDataByCol", "H5ADMatrix",
    function(x, j) extractNonzeroDataByCol(x@seed, j)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to dgCMatrix
###

.from_H5ADMatrix_to_dgCMatrix <- function(from) as(from@seed, "dgCMatrix")
setAs("H5ADMatrix", "dgCMatrix", .from_H5ADMatrix_to_dgCMatrix)
setAs("H5ADMatrix", "sparseMatrix", .from_H5ADMatrix_to_dgCMatrix)

