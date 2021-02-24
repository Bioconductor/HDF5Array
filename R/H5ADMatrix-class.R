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

setMethod("DelayedArray", "H5ADMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="H5ADMatrix")
)

### Works directly on an H5ADMatrixSeed derivative, in which case it must
### be called with a single argument.
H5ADMatrix <- function(filepath, name="X")
{
    if (is(filepath, "H5ADMatrixSeed")) {
        if (!identical(name, "X"))
            stop(wmsg("H5ADMatrix() must be called with a single argument ",
                      "when passed an H5ADMatrixSeed derivative"))
        seed <- filepath
    } else {
        seed <- H5ADMatrixSeed(filepath, name=name)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Taking advantage of sparsity
###

### Will work only if the seed is an H5SparseMatrixSeed derivative, that is,
### if it's a CSC_H5ADMatrixSeed or CSR_H5ADMatrixSeed object.
setMethod("sparsity", "H5ADMatrix", function(x) sparsity(x@seed))

setMethod("read_sparse_block", "H5ADMatrix",
    function(x, viewport) read_sparse_block(x@seed, viewport)
)

### Will work only if the seed is a CSC_H5ADMatrixSeed object.
setMethod("extractNonzeroDataByCol", "H5ADMatrix",
    function(x, j) extractNonzeroDataByCol(x@seed, j)
)

### Will work only if the seed is a CSR_H5ADMatrixSeed object.
setMethod("extractNonzeroDataByRow", "H5ADMatrix",
    function(x, i) extractNonzeroDataByCol(x@seed, i)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to dgCMatrix
###

.from_H5ADMatrix_to_dgCMatrix <- function(from) as(from@seed, "dgCMatrix")
setAs("H5ADMatrix", "dgCMatrix", .from_H5ADMatrix_to_dgCMatrix)
setAs("H5ADMatrix", "sparseMatrix", .from_H5ADMatrix_to_dgCMatrix)

