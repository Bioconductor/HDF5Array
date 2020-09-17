### =========================================================================
### h5mread()
### -------------------------------------------------------------------------
###
### An alternative to rhdf5::h5read() -- STILL EXPERIMENTAL!
###

### NOT exported!
check_selection <- function(dim, starts=NULL, counts=NULL)
{
    .Call2("C_check_selection", dim, starts, counts, PACKAGE="HDF5Array")
}

### NOT exported!
check_ordered_selection <- function(dim, starts=NULL, counts=NULL)
{
    .Call2("C_check_ordered_selection", dim, starts, counts,
           PACKAGE="HDF5Array")
}

### The selection must be strictly ascending along each dimension.
### Return NULL if the selection could not be reduced.
### Typical usage:
###     reduced <- reduce_selection(dim, starts, counts)
###     if (!is.null(reduced)) {
###         starts <- reduced[[1L]]
###         counts <- reduced[[2L]]
###     }
### NOT exported!
reduce_selection <- function(dim, starts=NULL, counts=NULL)
{
    .Call2("C_reduce_selection", dim, starts, counts, PACKAGE="HDF5Array")
}

### NOT exported!
map_starts_to_chunks <- function(starts, dim, chunkdim)
{
    .Call2("C_map_starts_to_chunks", starts, dim, chunkdim, PACKAGE="HDF5Array")
}

### When both 'starts' and 'counts' are specified, the selection must be
### strictly ascending along each dimension.
### By default the supplied selection is checked and reduced (if it can be).
### Set 'noreduce' to TRUE to skip the reduction step.
### Set 'as.integer' to TRUE to force returning the result as an integer array.
h5mread <- function(filepath, name, starts=NULL, counts=NULL, noreduce=FALSE,
                    as.integer=FALSE, as.sparse=FALSE, method=0L)
{
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
    ## C_h5mread will return an ordinary array if 'as.sparse' is FALSE,
    ## or a list of length 2 if it's TRUE. If the latter:
    ## - 'ans[[1L]]' will be an atomic vector of length 0 with the type of
    ##   the sparse data.
    ## - 'ans[[2L]]' will be a list of length the total nb of chunks touched
    ##   by the selection where each list element is itself a list of length 2.
    ##   This list of length 2 is a sparse representation of the data loaded
    ##   from the corresponding touched chunk i.e. its 1st element is the
    ##   "nzindex" (matrix) and its 2nd element the "nzdata" (atomic vector).
    ## - 'ans[[3L]]' will be an integer vector containing the dimensions of
    ##   the returned sparse data.
    ans <- .Call2("C_h5mread", filepath, name, starts, counts, noreduce,
                               as.integer, as.sparse, method,
                               PACKAGE="HDF5Array")
    if (as.sparse) {
        sparse_data_list <- ans[[2L]]
        if (length(sparse_data_list) == 0L) {
            nzindex <- NULL
            nzdata <- ans[[1L]]
        } else {
            ## Concatenate all the sparse data in 'sparse_data_list'.
            nzindex <- do.call(rbind, lapply(sparse_data_list, `[[`, 1L))
            nzdata <- unlist(lapply(sparse_data_list, `[[`, 2L))
        }
        ans <- SparseArraySeed(ans[[3L]], nzindex, nzdata, check=FALSE)
    }
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

