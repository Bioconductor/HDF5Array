### =========================================================================
### H5SparseMatrix objects
### -------------------------------------------------------------------------
###


setClass("H5SparseMatrix",
    contains="DelayedMatrix",
    representation(seed="H5SparseMatrixSeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "H5SparseMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="H5SparseMatrix")
)

### Works directly on an H5SparseMatrixSeed derivative, in which case it must
### be called with a single argument.
H5SparseMatrix <- function(filepath, group)
{
    if (is(filepath, "H5SparseMatrixSeed")) {
        if (!missing(group))
            stop(wmsg("H5SparseMatrix() must be called with a single argument ",
                      "when passed an H5SparseMatrixSeed object"))
        seed <- filepath
    } else {
        seed <- H5SparseMatrixSeed(filepath, group)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Taking advantage of sparsity
###

setMethod("sparsity", "H5SparseMatrix", function(x) sparsity(x@seed))

setMethod("read_sparse_block", "H5SparseMatrix",
    function(x, viewport) read_sparse_block(x@seed, viewport)
)

setMethod("extractNonzeroDataByCol", "H5SparseMatrix",
    function(x, j) extractNonzeroDataByCol(x@seed, j)
)

setMethod("extractNonzeroDataByRow", "H5SparseMatrix",
    function(x, i) extractNonzeroDataByCol(x@seed, i)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to dgCMatrix
###

.from_H5SparseMatrix_to_dgCMatrix <- function(from) as(from@seed, "dgCMatrix")
setAs("H5SparseMatrix", "dgCMatrix", .from_H5SparseMatrix_to_dgCMatrix)
setAs("H5SparseMatrix", "sparseMatrix", .from_H5SparseMatrix_to_dgCMatrix)

