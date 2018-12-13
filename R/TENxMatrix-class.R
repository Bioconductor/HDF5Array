### =========================================================================
### TENxMatrix objects
### -------------------------------------------------------------------------
###
### Note that we could just wrap a TENxMatrixSeed object in a DelayedMatrix
### object to represent and manipulate a 10x Genomics dataset as a
### DelayedMatrix object. So, strictly speaking, we don't really need the
### TENxMatrix class. However, we define this class mostly for cosmetic
### reasons, that is, to hide the DelayedMatrix class from the user.
### So the user will see and manipulate TENxMatrix objects instead of
### DelayedMatrix objects.
###


setClass("TENxMatrix",
    contains="DelayedMatrix",
    representation(seed="TENxMatrixSeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "TENxMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="TENxMatrix")
)

### Works directly on a TENxMatrixSeed object, in which case it must be
### called with a single argument.
TENxMatrix <- function(filepath, group="mm10")
{
    if (is(filepath, "TENxMatrixSeed")) {
        if (!missing(group))
            stop(wmsg("TENxMatrix() must be called with a single argument ",
                      "when passed a TENxMatrixSeed object"))
        seed <- filepath
    } else {
        seed <- TENxMatrixSeed(filepath, group)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Taking advantage of sparsity
###

setMethod("sparsity", "TENxMatrix", function(x) sparsity(x@seed))

setMethod("read_sparse_block", "TENxMatrix",
    function(x, viewport) read_sparse_block(x@seed, viewport)
)

setMethod("extractNonzeroDataByCol", "TENxMatrix",
    function(x, j) extractNonzeroDataByCol(x@seed, j)
)

