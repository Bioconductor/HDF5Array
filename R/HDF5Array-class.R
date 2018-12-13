### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------
###
### Note that we could just wrap an HDF5ArraySeed object in a DelayedArray
### object to represent and manipulate an HDF5 dataset as a DelayedArray
### object. So, strictly speaking, we don't really need the HDF5Array and
### HDF5Matrix classes. However, we define these classes mostly for cosmetic
### reasons, that is, to hide the DelayedArray and DelayedMatrix classes
### from the user. So the user will see and manipulate HDF5Array and
### HDF5Matrix objects instead of DelayedArray and DelayedMatrix objects.
###


setClass("HDF5Array",
    contains="DelayedArray",
    representation(seed="HDF5ArraySeed")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "HDF5ArraySeed",
    function(seed) new_DelayedArray(seed, Class="HDF5Array")
)

### Works directly on an HDF5ArraySeed object, in which case it must be
### called with a single argument.
HDF5Array <- function(filepath, name, type=NA)
{
    if (is(filepath, "HDF5ArraySeed")) {
        if (!(missing(name) && identical(type, NA)))
            stop(wmsg("HDF5Array() must be called with a single argument ",
                      "when passed an HDF5ArraySeed object"))
        seed <- filepath
    } else {
        seed <- HDF5ArraySeed(filepath, name, type=type)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5Matrix objects
###

setClass("HDF5Matrix", contains=c("HDF5Array", "DelayedMatrix"))

### Required for DelayedArray internal business.
setMethod("matrixClass", "HDF5Array", function(x) "HDF5Matrix")

### Automatic coercion method from HDF5Array to HDF5Matrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods
### don't bother to validate the object they return). So we overwrite it.
setAs("HDF5Array", "HDF5Matrix", function(from) new("HDF5Matrix", from))

### The user should not be able to degrade an HDF5Matrix object to
### an HDF5Array object so 'as(x, "HDF5Array", strict=TRUE)' should
### fail or be a no-op when 'x' is an HDF5Matrix object. Making this
### coercion a no-op seems to be the easiest (and safest) way to go.
setAs("HDF5Matrix", "HDF5Array", function(from) from)  # no-op

setAs("ANY", "HDF5Matrix",
    function(from) as(as(from, "HDF5Array"), "HDF5Matrix")
)

