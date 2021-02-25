### =========================================================================
### A wrapper to rhdf5::h5ls() that works on H5File objects
### -------------------------------------------------------------------------
###

h5ls <- function(file, recursive=TRUE, all=FALSE, datasetinfo=TRUE,
                 index_type=h5default("H5_INDEX"), order=h5default("H5_ITER"),
                 s3=FALSE, s3credentials=NULL, native=FALSE)
{
    if (is(file, "H5File")) {
        if (!(identical(s3, FALSE)
              && is.null(s3credentials)
              && identical(native, FALSE)))
            stop(wmsg("arguments 's3', 's3credentials', and 'native', ",
                      "cannot be used when 'file' is an H5File object"))
        if (file@no_rhdf5_h5id)
            stop(wmsg("this H5File object is not compatible with h5ls()"))
        file <- as(file, "H5IdComponent")
    }
    rhdf5::h5ls(file, recursive=recursive, all=all, datasetinfo=datasetinfo,
                      index_type=index_type, order=order,
                      s3=s3, s3credentials=s3credentials, native=native)
}

