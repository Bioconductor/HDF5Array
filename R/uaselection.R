### =========================================================================
### Manipulation of a user-supplied array selection
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

check_uaselection <- function(dim, starts=NULL, counts=NULL)
{
    .Call2("C_check_uaselection", dim, starts, counts, PACKAGE="HDF5Array")
}

check_ordered_uaselection <- function(dim, starts=NULL, counts=NULL)
{
    .Call2("C_check_ordered_uaselection", dim, starts, counts,
           PACKAGE="HDF5Array")
}

### The selection to reduce must be strictly ascending along each dimension.
### Return NULL if the selection could not be reduced.
### Typical usage:
###     reduced <- reduce_uaselection(dim, starts, counts)
###     if (!is.null(reduced)) {
###         starts <- reduced[[1L]]
###         counts <- reduced[[2L]]
###     }
reduce_uaselection <- function(dim, starts=NULL, counts=NULL)
{
    .Call2("C_reduce_uaselection", dim, starts, counts, PACKAGE="HDF5Array")
}

map_starts_to_chunks <- function(starts, dim, chunkdim)
{
    .Call2("C_map_starts_to_chunks", starts, dim, chunkdim, PACKAGE="HDF5Array")
}

