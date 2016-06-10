### =========================================================================
### Bind DelayedArray objects along their rows or columns
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ArrayBinder objects
###
### This class is for internal use only and is not exported.
###
 
setClass("ArrayBinder",
    representation(
        seeds="list",    # List of array-like objects to bind. Each object
                         # is expected to satisfy the "seed contract" i.e.
                         # to support dim(), dimnames(), and
                         # subset_seed_as_array().

        along="integer"  # Single integer indicating the dimension along
                         # which to bind the seeds.
    ),
    prototype(
        seeds=list(new("array")),
        along=1L
    )
)

.validate_ArrayBinder <- function(x)
{
    if (length(x@seeds) == 0L)
        return(wmsg("'x@seeds' cannot be empty"))
    if (!(isSingleInteger(x@along) && x@along > 0L))
        return(wmsg("'x@along' must be a single positive integer"))
    dims <- IRanges:::get_dims_to_bind(x@seeds, x@along)
    if (is.character(dims))
        return(wmsg(dims))
    TRUE
}

setValidity2("ArrayBinder", .validate_ArrayBinder)

.new_ArrayBinder <- function(seeds, along)
    new2("ArrayBinder", seeds=seeds, along=along)

### Implement the "seed contract" i.e. dim(), dimnames(), and
### subset_seed_as_array().

.get_ArrayBinder_dim <- function(x)
{
    dims <- IRanges:::get_dims_to_bind(x@seeds, x@along)
    IRanges:::combine_dims_along(dims, x@along)
}

setMethod("dim", "ArrayBinder", .get_ArrayBinder_dim)

.get_ArrayBinder_dimnames <- function(x)
{
    dims <- IRanges:::get_dims_to_bind(x@seeds, x@along)
    IRanges:::combine_dimnames_along(x@seeds, dims, x@along)
}

setMethod("dimnames", "ArrayBinder", .get_ArrayBinder_dimnames)

.subset_ArrayBinder_as_array <- function(x, index)
{
    dims <- IRanges:::get_dims_to_bind(x@seeds, x@along)
    breakpoints <- cumsum(dims[x@along, ])
    part_idx <- get_part_index(index[[x@along]], breakpoints)
    split_part_idx <- split_part_index(part_idx, length(breakpoints))
    FUN <- function(i) {
        index[[x@along]] <- split_part_idx[[i]]
        subset_seed_as_array(x@seeds[[i]], index)
    }
    tmp <- lapply(seq_along(x@seeds), FUN)

    ## Bind the ordinary arrays in 'tmp'.
    ans <- do.call(IRanges:::simple_abind, c(tmp, list(along=x@along)))

    ## Reorder the rows or columns in 'ans'.
    subscript <- rep.int(alist(foo=), length(index))
    subscript[[x@along]] <- get_rev_index(part_idx)
    subset_array_like_by_list(ans, subscript)
}

setMethod("subset_seed_as_array", "ArrayBinder",
    function(seed, index) .subset_ArrayBinder_as_array(seed, index)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### arbind() and acbind()
###

.DelayedArray_arbind <- function(...)
{
    objects <- unname(list(...))
    dims <- IRanges:::get_dims_to_bind(objects, 1L)
    if (is.character(dims))
        stop(wmsg(dims))
    DelayedArray(.new_ArrayBinder(objects, 1L))
}

.DelayedArray_acbind <- function(...)
{
    objects <- unname(list(...))
    dims <- IRanges:::get_dims_to_bind(objects, 2L)
    if (is.character(dims))
        stop(wmsg(dims))
    DelayedArray(.new_ArrayBinder(objects, 2L))
}

setMethod("arbind", "DelayedArray", .DelayedArray_arbind)
setMethod("acbind", "DelayedArray", .DelayedArray_acbind)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rbind() and cbind()
###

setMethod("rbind", "DelayedMatrix", .DelayedArray_arbind)
setMethod("cbind", "DelayedMatrix", .DelayedArray_acbind)

.as_DelayedMatrix_objects <- function(objects)
{
    lapply(objects,
        function(object) {
            if (length(dim(object)) != 2L)
                stop(wmsg("cbind() and rbind() are not supported on ",
                          "DelayedArray objects that don't have exactly ",
                          "2 dimensions. Please use acbind() or arnind() ",
                          "instead."))
            as(object, "DelayedMatrix")
        })
}

.DelayedArray_rbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("rbind", objects)
}

.DelayedArray_cbind <- function(...)
{
    objects <- .as_DelayedMatrix_objects(list(...))
    do.call("cbind", objects)
}

setMethod("rbind", "DelayedArray", .DelayedArray_rbind)
setMethod("cbind", "DelayedArray", .DelayedArray_cbind)

