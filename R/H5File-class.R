### =========================================================================
### H5File objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .h5openfile() / .h5closefile()
###
### The .h5open*file() functions return an h5 id (hid_t value) as a string.
###
### The .h5closefile() function returns NULL.
###

.h5openlocalfile <- function(filepath, readonly=TRUE)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying ",
                  "the path to an HDF5 file"))
    if (!isTRUEorFALSE(readonly))
        stop(wmsg("'readonly' must be TRUE or FALSE"))
    filepath <- file_path_as_absolute(filepath)
    .Call2("C_h5openlocalfile", filepath, readonly, PACKAGE="HDF5Array")
}

.h5openS3file <- function(filepath, s3credentials=NULL)
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
    .Call2("C_h5openS3file", filepath, auth, aws_region, secret_id, secret_key,
                             PACKAGE="HDF5Array")
}

.h5openfile <- function(filepath, s3=FALSE, s3credentials=NULL)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string specifying ",
                  "the path or URL to an HDF5 file"))
    if (!isTRUEorFALSE(s3))
        stop(wmsg("'s3' must be TRUE or FALSE"))
    if (s3) {
        ID <- .h5openS3file(filepath, s3credentials=s3credentials)
    } else {
        ID <- .h5openlocalfile(filepath, readonly=TRUE)
    }
    ID
}

.h5closefile <- function(ID)
{
    .Call2("C_h5closefile", ID, PACKAGE="HDF5Array")
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

.open_H5FileID_xp <- function(xp, filepath, s3=FALSE, s3credentials=NULL)
{
    ID <- .get_H5FileID_xp_ID(xp)
    if (!(is.null(ID) || is.na(ID))) {
        ## H5FileID object is already open and needs to be closed first.
        return(FALSE)
    }
    ID <- .h5openfile(filepath, s3=s3, s3credentials=s3credentials)
    .set_H5FileID_xp_ID(xp, ID)
    TRUE
}

.close_H5FileID_xp <- function(xp)
{
    ID <- .get_H5FileID_xp_ID(xp)
    if (is.null(ID) || is.na(ID)) {
        ## H5FileID object is already closed and needs to be opened first.
        return(FALSE)
    }
    .h5closefile(ID)
    .set_H5FileID_xp_ID(xp, NA_character_)
    TRUE
}

setClass("H5FileID",
    representation(
        xp="externalptr"
    )
    #prototype(
    #    xp=.new_H5FileID_xp(NA_character_)  # cannot be called at load time!
    #)
)

open.H5FileID <- function(con, ...)
{
    .open_H5FileID_xp(con@xp, ...)
}

close.H5FileID <- function(con, ...)
{
    .close_H5FileID_xp(con@xp, ...)
}

H5FileID <- function(filepath, s3=FALSE, s3credentials=NULL)
{
    ID <- .h5openfile(filepath, s3=s3, s3credentials=s3credentials)
    xp <- .new_H5FileID_xp(ID)
    reg.finalizer(xp, .close_H5FileID_xp, onexit=TRUE)
    new2("H5FileID", xp=xp)
}

setMethod("show", "H5FileID",
    function(object)
        cat("H5FileID: ", .get_H5FileID_xp_ID(object@xp), "\n", sep="")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### H5File objects
###

setClass("H5File",
    representation(
        h5fid="H5FileID",
        filepath="character",
        s3="logical"
    )
)

open.H5File <- function(con, ...)
{
    open(con@h5fid, con@filepath, s3=con@s3, ...)
}

close.H5File <- function(con, ...)
{
    close(con@h5fid, ...)
}

H5File <- function(filepath, s3=FALSE, s3credentials=NULL)
{
    h5fid <- H5FileID(filepath, s3=s3, s3credentials=s3credentials)
    if (!s3)
        filepath <- file_path_as_absolute(filepath)
    new2("H5File", h5fid=h5fid, filepath=filepath, s3=s3)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Union class: character_OR_H5File
###

setClassUnion("character_OR_H5File", c("character", "H5File"))

