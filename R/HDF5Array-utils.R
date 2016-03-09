### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A low-level utility for putting HDF5Array object in a "straight" form
###
### Untranspose the HDF5Array object and put its rows and columns in their
### "native" order. The goal is to put the matrix elements in their "native"
### order (i.e. in the same order as on disk) so as.vector() is faster on the
### resulting object is faster than on the original object.
###

.straighten_index <- function(i)
{
    i_len <- length(i)
    if (i_len == 0L)
        return(i)
    i_max <- max(i)
    ## Threshold is a rough estimate obtained empirically.
    ## TODO: Refine this.
    if (i_max <= 2L * i_len * log(i_len))
        which(as.logical(tabulate(i, nbins=i_max)))
    else
        sort(unique(i))
}

.straighten <- function(x, untranspose=FALSE, straighten.index=FALSE)
{
    if (is.array(x))
        return(x)
    if (untranspose)
        x@is_transposed <- FALSE
    if (!straighten.index)
        return(x)
    x_index <- index(x)
    index_was_touched <- FALSE
    for (n in seq_along(x_index)) {
        if (isStrictlySorted(x_index[[n]]))
            next
        x_index[[n]] <- .straighten_index(x_index[[n]])
        index_was_touched <- TRUE
    }
    if (index_was_touched)
        index(x) <- x_index
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

.HDF5Array_block_anyNA <- function(x, recursive=FALSE)
{
    x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
    blocks <- ArrayBlocks(dim(x), get_block_length(type(x)))
    for (i in seq_along(blocks)) {
        subarray <- extract_array_block(x, blocks, i)
        if (anyNA(as.vector(subarray)))
            return(TRUE)
    }
    FALSE
}

setMethod("anyNA", "HDF5Array", .HDF5Array_block_anyNA)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summary group generic
###
### max(), min(), range(), sum(), prod(), any(), all().
###

setMethod("Summary", "HDF5Array",
    function(x, ..., na.rm=FALSE)
    {
        if (missing(x)) {
            objects <- list(...)
        } else {
            objects <- list(x, ...)
        }
        ans <- suppressWarnings(callGeneric(NULL))  # init value
        for (x in objects) {
            if (.Generic %in% c("sum", "prod")) {
                x <- .straighten(x, untranspose=TRUE)
            } else {
                x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
            }
            blocks <- ArrayBlocks(dim(x), get_block_length(type(x)))
            for (i in seq_along(blocks)) {
                subarray <- extract_array_block(x, blocks, i)
                subarray_ans <- callGeneric(as.vector(subarray), na.rm=na.rm)
                ## Early bailout for any() and all().
                if (.Generic == "any") {
                    if (identical(subarray_ans, TRUE))
                        return(TRUE)
                } else if (.Generic == "all") {
                    if (identical(subarray_ans, FALSE))
                        return(FALSE)
                }
                ans <- callGeneric(ans, subarray_ans)
            }
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### Same arguments as base::mean.default().
.HDF5Array_block_mean <- function(x, trim=0, na.rm=FALSE)
{
    if (!identical(trim, 0))
        stop("\"mean\" method for HDF5Array objects ",
             "does not support the 'trim' argument yet")
    x <- .straighten(x, untranspose=TRUE)
    blocks <- ArrayBlocks(dim(x), get_block_length(type(x)))
    sum <- nval <- 0
    for (i in seq_along(blocks)) {
        subarray <- extract_array_block(x, blocks, i)
        tmp <- as.vector(subarray, mode="numeric")
        subarray_sum <- sum(tmp, na.rm=na.rm)
        if (is.na(subarray_sum))
            return(subarray_sum)
        subarray_nval <- if (na.rm) sum(!is.na(tmp)) else length(tmp)
        sum <- sum + subarray_sum
        nval <- nval + subarray_nval
    }
    sum / nval
}

### S3/S4 combo for mean.HDF5Array
mean.HDF5Array <- function(x, trim=0, na.rm=FALSE, ...)
    .HDF5Array_block_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "HDF5Array", .HDF5Array_block_mean)

