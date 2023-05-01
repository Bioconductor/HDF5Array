### =========================================================================
### H5File objects
### -------------------------------------------------------------------------


.ID_to_H5IdComponent <- function(ID) new("H5IdComponent", ID=ID, native=FALSE)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .h5openfile() / .h5closefile()
###
### The .h5open*file() functions return an h5 ID (hid_t value) as a string.
###
### The .h5closefile() function returns NULL.
###

.h5openlocalfile <- function(filepath, readonly=TRUE, use.rhdf5=FALSE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying ",
                  "the path to an HDF5 file"))
    if (!isTRUEorFALSE(readonly))
        stop(wmsg("'readonly' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(use.rhdf5))
        stop(wmsg("'use.rhdf5' must be TRUE or FALSE"))

    filepath <- file_path_as_absolute(filepath)
    if (use.rhdf5) {
        flags <- if (readonly) "H5F_ACC_RDONLY" else "H5F_ACC_RDWR"
        rhdf5::H5Fopen(filepath, flags)@ID
    } else {
        .Call2("C_h5openlocalfile", filepath, readonly, PACKAGE="HDF5Array")
    }
}

.h5openS3file <- function(filepath, s3credentials=NULL, use.rhdf5=FALSE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying ",
                  "the URL to an HDF5 file located on S3"))
    if (is.null(s3credentials)) {
        auth <- FALSE
        aws_region <- secret_id <- secret_key <- ""
    } else if (is.list(s3credentials) && length(s3credentials) == 3L) {
        auth <- TRUE
        aws_region <- s3credentials[[1L]]
        secret_id <- s3credentials[[2L]]
        secret_key <- s3credentials[[3L]]
    } else {
        stop(wmsg("'s3credentials' must be NULL or a list of 3 strings: ",
                  "(1) aws_region, (2) secret_id, (3) secret_key"))
    }
    if (!isTRUEorFALSE(use.rhdf5))
        stop(wmsg("'use.rhdf5' must be TRUE or FALSE"))

    if (use.rhdf5) {
        fapl_id <- rhdf5::H5Pcreate("H5P_FILE_ACCESS")
        on.exit(rhdf5::H5Pclose(fapl_id))
        rhdf5::H5Pset_fapl_ros3(fapl_id, s3credentials)
        loc <- rhdf5:::h5checktypeOrOpenLocS3(filepath, readonly=TRUE,
                                              fapl=fapl_id, native=FALSE)
        loc$H5Identifier@ID
    } else {
        .Call2("C_h5openS3file", filepath, auth,
                                 aws_region, secret_id, secret_key,
                                 PACKAGE="HDF5Array")
    }
}

.h5openfile <- function(filepath, s3=FALSE, s3credentials=NULL, use.rhdf5=FALSE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying ",
                  "the path or URL to an HDF5 file"))
    if (!isTRUEorFALSE(s3))
        stop(wmsg("'s3' must be TRUE or FALSE"))

    if (s3) {
        ID <- .h5openS3file(filepath, s3credentials=s3credentials,
                                      use.rhdf5=use.rhdf5)
    } else {
        ID <- .h5openlocalfile(filepath, readonly=TRUE, use.rhdf5=use.rhdf5)
    }
    ID
}

.h5closefile <- function(ID, use.rhdf5=FALSE)
{
    if (!isTRUEorFALSE(use.rhdf5))
        stop(wmsg("'use.rhdf5' must be TRUE or FALSE"))

    if (use.rhdf5) {
        fid <- .ID_to_H5IdComponent(ID)
        rhdf5::H5Fclose(fid)
    } else {
        .Call2("C_h5closefile", ID, PACKAGE="HDF5Array")
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level manipulation of the H5FileID external pointer and ID behind it
###

### 'ID' must be a string or NA_character_.
.set_H5FileID_xp_ID <- function(xp, ID)
{
    .Call2("C_set_H5FileID_xp_ID", xp, ID, PACKAGE="HDF5Array")
}

### Return a string, or NA_character_, or NULL.
### Note that NULL is returned if the H5FileID object was not properly
### initialized e.g. if it was constructed with 'new("H5FileID")' instead
### of with the H5FileID() constructor function.
.get_H5FileID_xp_ID <- function(xp)
{
    .Call2("C_get_H5FileID_xp_ID", xp, PACKAGE="HDF5Array")
}

### 'ID' must be a string or NA_character_.
.new_H5FileID_xp <- function(ID)
{
    .Call2("C_new_H5FileID_xp", ID, PACKAGE="HDF5Array")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### H5FileID objects
###

.ID_is_closed <- function(ID) { is.null(ID) || is.na(ID) }

.open_H5FileID_xp <- function(xp, filepath, s3=FALSE, s3credentials=NULL,
                                            use.rhdf5=FALSE)
{
    ID <- .get_H5FileID_xp_ID(xp)
    if (!.ID_is_closed(ID)) {
        ## H5FileID object is already open.
        return(FALSE)
    }
    ID <- .h5openfile(filepath, s3=s3, s3credentials=s3credentials,
                                use.rhdf5=use.rhdf5)
    .set_H5FileID_xp_ID(xp, ID)
    TRUE
}

.close_H5FileID_xp <- function(xp, use.rhdf5=FALSE)
{
    ID <- .get_H5FileID_xp_ID(xp)
    if (.ID_is_closed(ID)) {
        ## H5FileID object is already closed.
        return(FALSE)
    }
    .h5closefile(ID, use.rhdf5=use.rhdf5)
    .set_H5FileID_xp_ID(xp, NA_character_)
    TRUE
}

setClass("H5FileID",
    representation(
        xp="externalptr",
        from_rhdf5="logical"
    ),
    prototype(
        #xp=.new_H5FileID_xp(NA_character_),  # cannot be called at load time!
        from_rhdf5=FALSE
    )
)

open.H5FileID <- function(con, ...)
{
    .open_H5FileID_xp(con@xp, ..., use.rhdf5=con@from_rhdf5)
}

close.H5FileID <- function(con, ...)
{
    .close_H5FileID_xp(con@xp, ..., use.rhdf5=con@from_rhdf5)
}

H5FileID <- function(filepath, s3=FALSE, s3credentials=NULL, use.rhdf5=FALSE)
{
    ID <- .h5openfile(filepath, s3=s3, s3credentials=s3credentials,
                                use.rhdf5=use.rhdf5)
    xp <- .new_H5FileID_xp(ID)
    reg.finalizer(xp,
                  function(e) .close_H5FileID_xp(e, use.rhdf5=use.rhdf5),
                  onexit=TRUE)
    new2("H5FileID", xp=xp, from_rhdf5=use.rhdf5)
}

setMethod("show", "H5FileID",
    function(object)
    {
        ID <- .get_H5FileID_xp_ID(object@xp)
        explain <- ""
        if (!is.null(ID)) {
            if (is.na(ID)) {
                explain <- " (was closed)"
            } else if (object@from_rhdf5) {
                explain <- " (from rhdf5)"
            }
        }
        cat("H5FileID: ", ID, explain, "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### H5File objects
###
### Unfortunately, because HDF5Array and rhdf5 are both **statically** linked
### to the hdf5 library (libhdf5.a in the Rhdf5lib package), h5 IDs returned
### by calls to H5Fopen() in HDF5Array's C code cannot be used in rhdf5's
### calls to the hdf5 lib and vice versa. This is a huge bummer!
###
### We work around this with a crazy hack: we store **two** h5 IDs in an
### H5File object:
###   - The 1st h5 ID is compatible with HDF5Array, that is, it's made by
###     HDF5Array for HDF5Array (so it can be used in HDF5Array's calls to
###     the hdf5 lib).
###   - The 2nd h5 ID is compatible with rhdf5, that is, it's made by rhdf5
###     for rhdf5 (so it can be used in hdf5's calls to the hdf5 lib).
###     Note that this ID is only needed by HDF5Array functions that call
###     rhdf5's C code (only HDF5Array::h5ls() at the moment).
###
### One drawback of this hack is that a call to H5File(..., s3=TRUE, ...)
### needs to authenticate twice so takes twice longer.
###

setClass("H5File",
    representation(
        filepath="character",       # URL or **absolute** path to a local HDF5
                                    # file so the object won't break when the
                                    # user changes the working directory (e.g.
                                    # with setwd()).
        s3="logical",
        HDF5Array_h5id="H5FileID",  # compatible with HDF5Array
        no_rhdf5_h5id="logical",    # TRUE or FALSE
        rhdf5_h5id="H5FileID"       # compatible with rhdf5
    )
)

setMethod("path", "H5File", function(object) object@filepath)

open.H5File <- function(con, ...)
{
    if (!con@no_rhdf5_h5id)
        open(con@rhdf5_h5id, con@filepath, s3=con@s3, ...)
    open(con@HDF5Array_h5id, con@filepath, s3=con@s3, ...)
}

close.H5File <- function(con, ...)
{
    if (!con@no_rhdf5_h5id)
        close(con@rhdf5_h5id, ...)
    close(con@HDF5Array_h5id, ...)
}

H5File <- function(filepath, s3=FALSE, s3credentials=NULL, .no_rhdf5_h5id=FALSE)
{
    if (!isTRUEorFALSE(.no_rhdf5_h5id))
        stop(wmsg("'.no_rhdf5_h5id' must be TRUE or FALSE"))

    HDF5Array_h5id <- H5FileID(filepath, s3=s3, s3credentials=s3credentials)
    if (.no_rhdf5_h5id) {
        rhdf5_h5id <- new("H5FileID")
    } else {
        rhdf5_h5id <- H5FileID(filepath, s3=s3, s3credentials=s3credentials,
                                         use.rhdf5=TRUE)
    }
    if (!s3)
        filepath <- file_path_as_absolute(filepath)
    new2("H5File", filepath=filepath, s3=s3,
                   HDF5Array_h5id=HDF5Array_h5id,
                   no_rhdf5_h5id=.no_rhdf5_h5id,
                   rhdf5_h5id=rhdf5_h5id)
}

setMethod("show", "H5File",
    function(object)
    {
        ID <- .get_H5FileID_xp_ID(object@HDF5Array_h5id@xp)
        if (.ID_is_closed(ID))
            cat("CLOSED ")
        cat("H5File object", if (object@s3) " for a file on S3" else "",
            ":\n", sep="")
        cat("  path: ", path(object), "\n", sep="")
    }
)

setAs("H5File", "H5IdComponent",
    function(from)
    {
        if (from@no_rhdf5_h5id)
            stop(wmsg("This H5File object is not compatible with rhdf5. ",
                      "Maybe it was created with ",
                      "'H5File(..., .no_rhdf5_h5id=TRUE)'?"))
        h5id <- from@rhdf5_h5id
        ID <- .get_H5FileID_xp_ID(h5id@xp)
        if (.ID_is_closed(ID))
            stop(wmsg("H5File object is closed"))
        .ID_to_H5IdComponent(ID)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Union class: character_OR_H5File
###

setClassUnion("character_OR_H5File", c("character", "H5File"))

