### =========================================================================
### h5mread()
### -------------------------------------------------------------------------
###
### An alternative to rhdf5::h5read() -- STILL EXPERIMENTAL!
###

### NOT exported!
check_selection <- function(dim, starts=NULL, counts=NULL)
{
    .Call("C_check_selection", dim, starts, counts, PACKAGE="HDF5Array")
}

### NOT exported!
check_ordered_selection <- function(dim, starts=NULL, counts=NULL)
{
    .Call("C_check_ordered_selection", dim, starts, counts, PACKAGE="HDF5Array")
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
    .Call("C_reduce_selection", dim, starts, counts, PACKAGE="HDF5Array")
}

### NOT exported!
map_starts_to_chunks <- function(starts, dim, chunkdim)
{
    .Call("C_map_starts_to_chunks", starts, dim, chunkdim, PACKAGE="HDF5Array")
}

### When both 'starts' and 'counts' are specified, the selection must be
### strictly ascending along each dimension.
### By default the supplied selection is checked and reduced (if it can be).
### Set 'noreduce' to TRUE to skip the reduction step.
### Set 'as.integer' to TRUE to force returning the result as an integer array.
h5mread <- function(filepath, name, starts=NULL, counts=NULL, noreduce=FALSE,
                    as.integer=FALSE, method=0L)
{
    if (is.null(starts)) {
        if (!is.null(counts))
            stop(wmsg("'counts' must be NULL when 'starts' is NULL"))
    } else if (is.list(starts)) {
        if (is.null(counts)) {
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
                function(start) is.null(start) || isStrictlySorted(start),
                logical(1))
            starts <- lapply(seq_along(starts0),
                function(i) {
                    start <- starts0[[i]]
                    if (ok[[i]])
                        return(start)
                    unique(sort(start))
                })
        }
    } else {
        stop(wmsg("'starts' must be a list (or NULL)"))
    }
    ans <- .Call("C_h5mread", filepath, name, starts, counts, noreduce,
                              as.integer, method,
                              PACKAGE="HDF5Array")
    if (is.null(starts) || !is.null(counts))
        return(ans)
    Nindex <- lapply(seq_along(starts0),
        function(i) {
            if (ok[[i]])
                return(NULL)
            match(starts0[[i]], starts[[i]])
        })
    DelayedArray:::subset_by_Nindex(ans, Nindex, drop=FALSE)
}

