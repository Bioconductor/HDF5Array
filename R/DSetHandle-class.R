### =========================================================================
### DSetHandle objects
### -------------------------------------------------------------------------


setClass("DSetHandle",
    representation(
        xp="externalptr"
    )
)

.destroy_DSetHandle_xp <- function(xp)
{
    .Call("C_destroy_DSetHandle_xp", xp, PACKAGE="HDF5Array")
}

DSetHandle <- function(filepath, name, as.integer=FALSE)
{
    xp <- .Call("C_create_DSetHandle_xp", filepath, name, as.integer,
                                          PACKAGE="HDF5Array")
    reg.finalizer(xp, .destroy_DSetHandle_xp, onexit=TRUE)
    new2("DSetHandle", xp=xp)
}

destroy_DSetHandle <- function(x)
{
    invisible(.destroy_DSetHandle_xp(x@xp))
}

setMethod("show", "DSetHandle",
    function(object)
        .Call("C_show_DSetHandle_xp", object@xp, PACKAGE="HDF5Array")
)

### The R type returned by h5mread() is determined by arguments 'filepath',
### 'name', and 'as.integer'.
get_h5mread_returned_type <- function(filepath, name, as.integer=FALSE)
{
    .Call("C_get_h5mread_returned_type", filepath, name, as.integer,
          PACKAGE="HDF5Array")
}

