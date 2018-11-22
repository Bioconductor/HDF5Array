### =========================================================================
### h5mread()
### -------------------------------------------------------------------------
###
### An alternative to rhdf5::h5read() -- STILL EXPERIMENTAL!
###

### The selection must be strictly ascending along each dimension.
### If 'dim' is provided, it must also be within the extend of each dimension.
### Return NULL if the selection could not be reduced.
### Typical usage:
###     reduced <- reduce_selection(starts, counts)
###     if (!is.null(reduced)) {
###         starts <- reduced[[1L]]
###         counts <- counts[[1L]]
###     }
reduce_selection <- function(starts, counts=NULL, dim=NULL)
{
    .Call("C_reduce_selection", starts, counts, dim, PACKAGE="HDF5Array")
}

### NOT exported!
map_starts_to_chunks <- function(starts, dim, chunkdim)
{
    .Call("C_map_starts_to_chunks", starts, dim, chunkdim, PACKAGE="HDF5Array")
}

### The selection must be strictly ascending along each dimension.
### By default the supplied selection is checked and reduced (if it can be).
### Set 'noreduce' to TRUE to skip the reduction step.
### Set 'as.integer' to TRUE to force returning the result as an integer array.
h5mread <- function(filepath, name, starts, counts=NULL, noreduce=FALSE,
                    method=1L, as.integer=FALSE)
{
    .Call("C_h5mread", filepath, name, starts, counts, noreduce,
                       method, as.integer,
                       PACKAGE="HDF5Array")
}

