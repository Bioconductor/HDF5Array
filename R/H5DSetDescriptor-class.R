### =========================================================================
### H5DSetDescriptor objects
### -------------------------------------------------------------------------


setClass("H5DSetDescriptor",
    representation(
        xp="externalptr"
    )
)

.destroy_H5DSetDescriptor_xp <- function(xp)
{
    .Call2("C_destroy_H5DSetDescriptor_xp", xp, PACKAGE="HDF5Array")
}

H5DSetDescriptor <- function(filepath, name, as.integer=FALSE)
{
    xp <- .Call2("C_new_H5DSetDescriptor_xp", filepath, name, as.integer,
                                              PACKAGE="HDF5Array")
    reg.finalizer(xp, .destroy_H5DSetDescriptor_xp, onexit=TRUE)
    new2("H5DSetDescriptor", xp=xp)
}

destroy_H5DSetDescriptor <- function(x)
{
    invisible(.destroy_H5DSetDescriptor_xp(x@xp))
}

setMethod("show", "H5DSetDescriptor",
    function(object)
        .Call2("C_show_H5DSetDescriptor_xp", object@xp, PACKAGE="HDF5Array")
)

### The R type returned by h5mread() is determined by arguments 'filepath',
### 'name', and 'as.integer'.
get_h5mread_returned_type <- function(filepath, name, as.integer=FALSE)
{
    .Call2("C_get_h5mread_returned_type", filepath, name, as.integer,
           PACKAGE="HDF5Array")
}

