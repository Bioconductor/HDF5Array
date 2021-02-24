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
    if (!is(filepath, "H5File")) {
        filepath <- H5File(filepath)
        on.exit(close(filepath))
    }
    name <- normarg_h5_name(name)

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

