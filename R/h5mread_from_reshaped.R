### =========================================================================
### h5mread_from_reshaped()
### -------------------------------------------------------------------------
###
### An h5mread() wrapper that reads data from a virtually reshaped
### HDF5 dataset.
###

.INVALID_RESHAPING_MSG <- c(
    "Reshaping only supports reducing the number of dimensions by ",
    "collapsing a group of consecutive dimensions into a single ",
    "dimension (e.g. reshaping a 10 x 3 x 5 x 1000 array as a ",
    "10 x 15 x 1000 array or as a 150 x 1000 matrix)."
)

find_dims_to_collapse <- function(dim, dim0)
{
    if (prod(dim) != prod(dim0))
        stop(wmsg("Reshaping must preserve the length of the HDF5 dataset. ",
                  "More precisely, 'prod(dim)' must be equal to 'prod(dim0)' ",
                  "where 'dim0' is the vector of physical dimensions ",
                  "of the HDF5 dataset."))
    ndim <- length(dim)
    ndim0 <- length(dim0)
    if (ndim > ndim0)
        stop(wmsg("Trying to set ", ndim, " dimensions on an HDF5 dataset ",
                  "with ", ndim0, " dimensions. Reshaping doesn't support ",
                  "increasing the number of dimensions at the moment."))
    if (ndim == ndim0) {
        if (all(dim == dim0))
            return(NULL)  # no reshaping
        stop(wmsg(.INVALID_RESHAPING_MSG))
    }
    idx <- which(dim != head(dim0, n=ndim))
    if (length(idx) == 0L)
        return(c(ndim, ndim0))
    along1 <- idx[[1L]]
    along2 <- along1 + ndim0 - ndim
    if (!all(tail(dim, n=ndim-along1) == tail(dim0, n=ndim0-along2)))
        stop(wmsg(.INVALID_RESHAPING_MSG))
    c(along1, along2)
}

collapse_dims <- function(dim0, collapse_along)
{
    if (is.null(collapse_along))
        return(dim0)
    along1 <- collapse_along[[1L]]
    along2 <- collapse_along[[2L]]
    c(dim0[seq_len(along1 - 1L)],
      prod(dim0[along1:along2]),
      dim0[seq_len(length(dim0) - along2) + along2])
}

.h5mread_and_collapse_dims <- function(filepath, name, starts, collapse_along,
                                       noreduce=FALSE, as.integer=FALSE,
                                       method=0L)
{
    ans <- h5mread(filepath, name, starts, noreduce=noreduce,
                   as.integer=as.integer, method=method)
    dim(ans) <- collapse_dims(dim(ans), collapse_along)
    ans
}

.Mindex_as_index_list <- function(Mindex)
{
    args <- c(lapply(2:ncol(Mindex), function(j) Mindex[ , j]), list(sep=","))
    rle <- Rle(do.call(paste, args))
    skeleton <- PartitioningByWidth(runLength(rle), names=runValue(rle))
    tmp1 <- relist(Mindex[ , 1L], skeleton)
    tmp2 <- strsplit(names(tmp1), ",", fixed=TRUE)
    lapply(seq_along(tmp1),
        function(i) c(list(tmp1[[i]]), as.list(as.integer(tmp2[[i]]))))
}

### 'dim' specifies how to reshape the HDF5 dataset to read from. Note that
### we do NOT support arbitrary reshaping (see .INVALID_RESHAPING_MSG at the
### top of this file for the kind of reshaping that is currently supported).
### 'starts' must be a multidimensional subsetting index with respect to
### the reshaped dataset. Note that the user gets the illusion that the
### reshaping happens **before** the data is read even though the dataset
### in the HDF5 dataset is not touched (it's treated as read-only).
h5mread_from_reshaped <- function(filepath, name, dim, starts, noreduce=FALSE,
                                  as.integer=FALSE, method=0L)
{
    dim <- DelayedArray:::normarg_dim(dim)
    dim0 <- h5dim(filepath, name)
    collapse_along <- find_dims_to_collapse(dim, dim0)
    ndim <- length(dim)
    stopifnot(is.list(starts), length(starts) == ndim)
    if (is.null(collapse_along)) {
        ## No reshaping.
        ans <- h5mread(filepath, name, starts, noreduce=noreduce,
                       as.integer=as.integer, method=method)
        return(ans)
    }
    along1 <- collapse_along[[1L]]
    along2 <- collapse_along[[2L]]
    idx0 <- along1:along2
    Lstarts <- starts[seq_len(along1 - 1L)]
    Rstarts <- starts[seq_len(ndim - along1) + along1]
    starts0 <- c(Lstarts, vector("list", length=length(idx0)), Rstarts)
    start1 <- starts[[along1]]
    if (is.null(start1)) {
        ans <- .h5mread_and_collapse_dims(filepath, name,
                                          starts0, collapse_along,
                                          noreduce=noreduce,
                                          as.integer=as.integer,
                                          method=method)
        return(ans)
    }
    ## Other 'starts' list elements will be checked by h5mread().
    if (!is.numeric(start1))
        stop(wmsg("each list element in 'starts' must ", 
                  "be NULL or a numeric vector"))
    if (length(start1) == 0L) {
        starts0[idx0] <- rep.int(list(integer(0)), length(idx0))
        ans <- .h5mread_and_collapse_dims(filepath, name,
                                          starts0, collapse_along,
                                          noreduce=noreduce,
                                          as.integer=as.integer,
                                          method=method)
        return(ans)
    }

    Mindex <- Lindex2Mindex(start1, dim0[idx0])
    index_list <- .Mindex_as_index_list(Mindex)
    tmp <- lapply(seq_along(index_list),
        function(i) {
            starts0[idx0] <- index_list[[i]]
            .h5mread_and_collapse_dims(filepath, name,
                                       starts0, collapse_along,
                                       noreduce=noreduce,
                                       as.integer=as.integer,
                                       method=method)
        })
    do.call(DelayedArray:::simple_abind, c(tmp, list(along=along1)))
}

