### =========================================================================
### HDF5Matrix objects
### -------------------------------------------------------------------------


### TODO: Try to also inherit from DataTable. Benefits?
setClass("HDF5Matrix",
    contains="HDF5Array",
    representation(
        ## x@rowdim and x@coldim are 2 integers such that
        ##     1 <= x@rowdim < x@coldim <= length(x@index)
        rowdim="integer",  # single integer
        coldim="integer"   # single integer
    )
)

.is_valid_rowdim <- function(rowdim, ndim)
{
    is.integer(rowdim) && isSingleNumber(rowdim) &&
        rowdim >= 1L && rowdim <= ndim
}

.validate_HDF5Matrix <- function(x)
{
    if (!.is_valid_rowdim(x@rowdim, length(x@index))
     || !.is_valid_rowdim(x@coldim, length(x@index)))
        return(wmsg("'x@rowdim' and 'x@coldim' must be single integers ",
                    ">= 1 and <= 'length(x@index)'"))
    if (x@rowdim >= x@coldim)
        return("'x@rowdim' must be < 'x@coldim'")
    array_dim <- get_HDF5Array_dim(x)
    if (!all(array_dim[-c(x@rowdim, x@coldim)] == 1L))
        return(wmsg("'x@rowdim' and 'x@coldim' are incompatible with the ",
                    "dimensions of the underlying HDF5Array object"))
    TRUE
}

setValidity2("HDF5Matrix", .validate_HDF5Matrix)

.get_HDF5Matrix_dim <- function(x)
{
    array_dim <- get_HDF5Array_dim(x)
    array_dim[c(x@rowdim, x@coldim)]
}

setMethod("dim", "HDF5Matrix", .get_HDF5Matrix_dim)

.get_HDF5Matrix_dimnames <- function(x)
{
    ans <- lapply(x@index[c(x@rowdim, x@coldim)], names)
    if (is.null(unlist(ans)))
        return(NULL)
    ans
}

setMethod("dimnames", "HDF5Matrix", .get_HDF5Matrix_dimnames)

.set_HDF5Matrix_dimnames <- function(x, value)
{
    value <- normalize_dimnames_replacement_value(value, 2L)
    ## 'x@index' can be big so avoid copies when possible. With this trick
    ## no-op dimnames(x) <- dimnames(x) is instantaneous.
    if (!identical(names(x@index[[x@rowdim]]), value[[1L]]))
        names(x@index[[x@rowdim]]) <- value[[1L]]
    if (!identical(names(x@index[[x@coldim]]), value[[2L]]))
        names(x@index[[x@coldim]]) <- value[[2L]]
    x
}

setReplaceMethod("dimnames", "HDF5Matrix", .set_HDF5Matrix_dimnames)

.from_HDF5Array_to_HDF5Matrix <- function(from)
{
    array_dim <- get_HDF5Array_dim(from)
    if (length(array_dim) < 2L)
        stop(wmsg(class(from), " object with less than 2 dimensions ",
                  "cannot be coerced to a HDF5Matrix object at the moment"))
    idx <- which(array_dim != 1L)
    if (length(idx) > 2L)
        stop(wmsg(class(from), " object with more than 2 effective dimensions ",
                  "cannot be coerced to a HDF5Matrix object. ", slicing_tip))
    if (length(idx) == 2L) {
        ans_rowdim <- idx[[1L]]
        ans_coldim <- idx[[2L]]
    } else if (length(idx) == 0L) {
        ans_rowdim <- 1L
        ans_coldim <- 2L
    } else {
        ## length(idx) == 1L
        ans_rowdim <- idx[[1L]]
        if (ans_rowdim == length(array_dim))
            stop(wmsg("A ", class(from), " object where the only effective ",
                      "dimension is its last dimension cannot be coerced ",
                      "to a HDF5Matrix object at the moment"))
        ans_coldim <- ans_rowdim + 1L
    }
    new2("HDF5Matrix", from, rowdim=ans_rowdim, coldim=ans_coldim)
}

setAs("HDF5Array", "HDF5Matrix", .from_HDF5Array_to_HDF5Matrix)

### Constructor.
HDF5Matrix <- function(file, group, name)
{
    hdf5array <- HDF5Array(file, group, name)
    as(hdf5array, "HDF5Matrix")
}

.extract_HDF5Matrix_subset <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")

    ## Check dimensionality of subsetting i.e whether it's 1D- or 2D- or
    ## 3D- etc... subsetting.
    ndim <- nargs() - 1L
    if (!missing(drop))
        ndim <- ndim - 1L
    if (ndim == 1L && missing(i))
        ndim <- 0L
    if (ndim != 0L && ndim != 2L) {
        if (ndim == 1L)
            stop("1D-subsetting is not supported")
        stop("incorrect number of dimensions")
    }

    ## Perform the subsetting.
    if (!missing(i))
        x@index[[x@rowdim]] <- extractROWS(x@index[[x@rowdim]], i)
    if (!missing(j))
        x@index[[x@coldim]] <- extractROWS(x@index[[x@coldim]], j)
    x
}

setMethod("[", "HDF5Matrix", .extract_HDF5Matrix_subset)

