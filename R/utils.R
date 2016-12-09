### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Simple wrappers around rhdf5:::h5read() and rhdf5:::h5write()
###

.replace_missing_subscripts_by_NULLS <- function(subscripts)
{
    ans <- vector("list", length(subscripts))
    not_null_idx <- which(vapply(subscripts, class, character(1)) != "name")
    ans[not_null_idx] <- subscripts[not_null_idx]
    ans
}

h5read2 <- function(file, name, index=NULL)
{
    if (!is.null(index))
        index <- .replace_missing_subscripts_by_NULLS(index)
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    suppressWarnings(h5read(file, name, index=index))
}

h5write2 <- function(obj, file, name, index=NULL)
{
    if (!is.null(index))
        index <- .replace_missing_subscripts_by_NULLS(index)
    h5write(obj, file, name, index=index)
}

