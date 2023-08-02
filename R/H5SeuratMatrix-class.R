### =========================================================================
### H5SeuratMatrix objects
### -------------------------------------------------------------------------
###
### Note that we could just wrap a H5SeuratMatrixSeed object in a DelayedMatrix
### object to represent and manipulate a 10x Genomics dataset as a
### DelayedMatrix object. So, strictly speaking, we don't really need the
### H5SeuratMatrix class. However, we define this class mostly for cosmetic
### reasons, that is, to hide the DelayedMatrix class from the user.
### So the user will see and manipulate H5SeuratMatrix objects instead of
### DelayedMatrix objects.
###


setClass("H5SeuratMatrix",
    contains="DelayedMatrix",
    representation(seed="H5SeuratMatrixSeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "H5SeuratMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="H5SeuratMatrix")
)

### Works directly on a TENxMatrixSeed object, in which case it must be
### called with a single argument.
H5SeuratMatrix <- function(filepath, group="matrix")
{
    if (is(filepath, "H5SeuratMatrixSeed")) {
        if (!missing(group))
            stop(wmsg("H5SeuratMatrix() must be called with a single argument ",
                      "when passed a H5SeuratMatrixSeed object"))
        seed <- filepath
    } else {
        seed <- H5SeuratMatrixSeed(filepath, group)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Taking advantage of sparsity
###

setMethod("sparsity", "H5SeuratMatrix", function(x) sparsity(x@seed))

setMethod("read_sparse_block", "H5SeuratMatrix",
    function(x, viewport) read_sparse_block(x@seed, viewport)
)

setMethod("extractNonzeroDataByCol", "H5SeuratMatrix",
    function(x, j) extractNonzeroDataByCol(x@seed, j)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to dgCMatrix
###

.from_H5SeuratMatrix_to_dgCMatrix <- function(from) as(from@seed, "dgCMatrix")
setAs("H5SeuratMatrix", "dgCMatrix", .from_H5SeuratMatrix_to_dgCMatrix)
setAs("H5SeuratMatrix", "sparseMatrix", .from_H5SeuratMatrix_to_dgCMatrix)

