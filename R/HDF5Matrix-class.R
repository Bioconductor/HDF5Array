### =========================================================================
### HDF5Matrix objects
### -------------------------------------------------------------------------


### TODO: Try to also inherit from DataTable. Benefits?
setClass("HDF5Matrix",
    contains="HDF5Array",
    representation(
        ## x@N1 and x@N2 must be 2 integers such that
        ##     1 <= x@N1 < x@N2 <= length(x@index)
        N1="integer",  # single integer
        N2="integer"   # single integer
    )
)

.is_valid_N1 <- function(N1, ndim)
{
    isSingleInteger(N1) && N1 >= 1L && N1 <= ndim
}

.validate_HDF5Matrix <- function(x)
{
    if (!.is_valid_N1(x@N1, length(x@index))
     || !.is_valid_N1(x@N2, length(x@index)))
        return(wmsg("'x@N1' and 'x@N2' must be single integers ",
                    ">= 1 and <= 'length(x@index)'"))
    if (x@N1 >= x@N2)
        return("'x@N1' must be < 'x@N2'")
    array_dim <- get_HDF5Array_dim(x)
    if (!all(array_dim[-c(x@N1, x@N2)] == 1L))
        return(wmsg("'x@N1' and 'x@N2' are incompatible with the ",
                    "dimensions of the underlying HDF5Array object"))
    TRUE
}

setValidity2("HDF5Matrix", .validate_HDF5Matrix)

.get_HDF5Matrix_dim <- function(x)
{
    ans <- lengths(x@index)[c(x@N1, x@N2)]
    if (x@transpose)
        ans <- rev(ans)
    ans
}

setMethod("dim", "HDF5Matrix", .get_HDF5Matrix_dim)

.get_HDF5Matrix_dimnames <- function(x)
{
    ans <- lapply(x@index[c(x@N1, x@N2)], names)
    if (is.null(unlist(ans)))
        return(NULL)
    if (x@transpose)
        ans <- rev(ans)
    ans
}

setMethod("dimnames", "HDF5Matrix", .get_HDF5Matrix_dimnames)

.set_HDF5Matrix_dimnames <- function(x, value)
{
    value <- normalize_dimnames_replacement_value(value, 2L)
    if (x@transpose)
        value <- rev(value)
    ## 'x@index' can be big so avoid copies when possible. With this trick
    ## no-op dimnames(x) <- dimnames(x) is instantaneous.
    if (!identical(names(x@index[[x@N1]]), value[[1L]]))
        names(x@index[[x@N1]]) <- value[[1L]]
    if (!identical(names(x@index[[x@N2]]), value[[2L]]))
        names(x@index[[x@N2]]) <- value[[2L]]
    x
}

setReplaceMethod("dimnames", "HDF5Matrix", .set_HDF5Matrix_dimnames)

.from_HDF5Array_to_HDF5Matrix <- function(from)
{
    dim0 <- lengths(from@index)
    if (length(dim0) < 2L)
        stop(wmsg(class(from), " object with less than 2 dimensions ",
                  "cannot be coerced to an HDF5Matrix object at the moment"))
    idx <- which(dim0 != 1L)
    if (length(idx) > 2L)
        stop(wmsg(class(from), " object with more than 2 effective dimensions ",
                  "cannot be coerced to an HDF5Matrix object. ", slicing_tip))
    if (length(idx) == 2L) {
        N1 <- idx[[1L]]
        N2 <- idx[[2L]]
    } else if (length(idx) == 0L) {
        N1 <- 1L
        N2 <- 2L
    } else {
        ## length(idx) == 1L
        N1 <- idx[[1L]]
        if (N1 == length(dim0))
            stop(wmsg("A ", class(from), " object where the only effective ",
                      "dimension is its last dimension cannot be coerced ",
                      "to a HDF5Matrix object at the moment"))
        N2 <- N1 + 1L
    }
    new2("HDF5Matrix", from, N1=N1, N2=N2)
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
            stop("1D-style subsetting is not supported")
        stop("incorrect number of dimensions")
    }

    ## Perform the subsetting.
    if (!missing(i)) {
        n <- if (x@transpose) x@N2 else x@N1
        x@index[[n]] <- extractROWS(x@index[[n]], i)
    }
    if (!missing(j)) {
        n <- if (x@transpose) x@N1 else x@N2
        x@index[[n]] <- extractROWS(x@index[[n]], j)
    }
    x
}

setMethod("[", "HDF5Matrix", .extract_HDF5Matrix_subset)

