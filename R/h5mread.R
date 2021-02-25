### =========================================================================
### h5mread()
### -------------------------------------------------------------------------
###
### An alternative to rhdf5::h5read() -- STILL EXPERIMENTAL!
###

### The R type returned by h5mread() is determined by arguments 'filepath',
### 'name', and 'as.integer'.
get_h5mread_returned_type <- function(filepath, name, as.integer=FALSE)
{
    if (!is(filepath, "H5File")) {
        filepath <- H5File(filepath, .no_rhdf5_h5id=TRUE)
        on.exit(close(filepath))
    }
    name <- normarg_h5_name(name)

    .Call2("C_get_h5mread_returned_type", filepath, name, as.integer,
           PACKAGE="HDF5Array")
}

### When both 'starts' and 'counts' are specified, the selection must be
### strictly ascending along each dimension.
### By default the user-supplied selection is checked and reduced (if it
### can be).
### Set 'noreduce' to TRUE to skip the reduction step.
### Set 'as.integer' to TRUE to force returning the result as an integer array.
h5mread <- function(filepath, name, starts=NULL, counts=NULL, noreduce=FALSE,
                    as.integer=FALSE, as.sparse=FALSE,
                    method=0L, use.H5Dread_chunk=FALSE)
{
    if (!is(filepath, "H5File")) {
        filepath <- H5File(filepath, .no_rhdf5_h5id=TRUE)
        on.exit(close(filepath))
    }
    name <- normarg_h5_name(name)

    if (!isTRUEorFALSE(as.sparse))
        stop(wmsg("'as.sparse' must be TRUE or FALSE"))
    if (is.null(starts)) {
        if (!is.null(counts))
            stop(wmsg("'counts' must be NULL when 'starts' is NULL"))
    } else if (is.list(starts)) {
        order_starts <- is.null(counts) &&
                        !all(S4Vectors:::sapply_isNULL(starts))
        if (order_starts) {
            ## Round the 'starts'.
            starts0 <- lapply(starts,
                function(start) {
                    if (is.null(start))
                        return(NULL)
                    if (!is.numeric(start))
                        stop(wmsg("each list element in 'starts' must ",
                                  "be NULL or a numeric vector"))
                    if (!is.integer(start))
                        start <- round(start)
                    start
                })
            ok <- vapply(starts0,
                function(start0) is.null(start0) || isStrictlySorted(start0),
                logical(1))
            order_starts <- !all(ok)
            if (order_starts) {
                starts <- lapply(seq_along(starts0),
                    function(i) {
                        start0 <- starts0[[i]]
                        if (ok[[i]])
                            return(start0)
                        start0 <- sort(start0)
                        start <- unique(start0)
                        if (as.sparse && length(start) != length(start0))
                            stop(wmsg("when 'as.sparse' is TRUE, list ",
                                      "elements in 'starts' are not allowed ",
                                      "to contain duplicates"))
                        start
                    })
            } else {
                starts <- starts0
            }
        }
    } else {
        stop(wmsg("'starts' must be a list (or NULL)"))
    }
    ## C_h5mread() will return an ordinary array if 'as.sparse' is FALSE,
    ## or 'list(ans_dim, nzindex, nzdata)' if it's TRUE.
    ans <- .Call2("C_h5mread", filepath, name, starts, counts, noreduce,
                               as.integer, as.sparse,
                               method, use.H5Dread_chunk,
                               PACKAGE="HDF5Array")
    if (as.sparse)
        ans <- SparseArraySeed(ans[[1L]], ans[[2L]], ans[[3L]], check=FALSE)
    if (is.null(starts) || !order_starts)
        return(ans)
    index <- lapply(seq_along(starts0),
        function(i) {
            if (ok[[i]])
                return(NULL)
            match(starts0[[i]], starts[[i]])
        })
    if (as.sparse) {
        extract_sparse_array(ans, index)
    } else {
        extract_array(ans, index)
    }
}

