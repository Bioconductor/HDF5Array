### =========================================================================
### ReshapedHDF5Array objects
### -------------------------------------------------------------------------
###


setClass("ReshapedHDF5Array",
    contains="HDF5Array",
    representation(seed="ReshapedHDF5ArraySeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "ReshapedHDF5ArraySeed",
    function(seed) new_DelayedArray(seed, Class="ReshapedHDF5Array")
)

### Works directly on a ReshapedHDF5ArraySeed object, in which case it must
### be called with a single argument.
ReshapedHDF5Array <- function(filepath, name, dim, type=NA)
{
    if (is(filepath, "ReshapedHDF5ArraySeed")) {
        if (!(missing(name) && missing(dim) && identical(type, NA)))
            stop(wmsg("ReshapedHDF5Array() must be called with a single ",
                      "argument when passed a ReshapedHDF5ArraySeed object"))
        seed <- filepath
    } else {
        seed <- ReshapedHDF5ArraySeed(filepath, name, dim, type=type)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ReshapedHDF5Matrix objects
###

setClass("ReshapedHDF5Matrix", contains=c("ReshapedHDF5Array", "HDF5Matrix"))

### Required for DelayedArray internal business.
setMethod("matrixClass", "ReshapedHDF5Array", function(x) "ReshapedHDF5Matrix")

### Automatic coercion method from ReshapedHDF5Array to ReshapedHDF5Matrix
### silently returns a broken object (unfortunately these dummy automatic
### coercion methods don't bother to validate the object they return).
### So we overwrite it.
setAs("ReshapedHDF5Array", "ReshapedHDF5Matrix",
    function(from) new("ReshapedHDF5Matrix", from)
)

### The user should not be able to degrade a ReshapedHDF5Matrix object to
### a ReshapedHDF5Array object so 'as(x, "ReshapedHDF5Array", strict=TRUE)'
### should fail or be a no-op when 'x' is a ReshapedHDF5Matrix object.
### Making this coercion a no-op seems to be the easiest (and safest) way
### to go.
setAs("ReshapedHDF5Matrix", "ReshapedHDF5Array", function(from) from)  # no-op

setAs("ANY", "ReshapedHDF5Matrix",
    function(from) as(as(from, "ReshapedHDF5Array"), "ReshapedHDF5Matrix")
)

