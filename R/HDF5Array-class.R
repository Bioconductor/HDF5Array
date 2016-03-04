### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


setClass("HDF5Array",
    representation(
        file="character",    # single string
        group="character",   # single string
        name="character",    # dataset name
        index="list",        # list of integer vectors
        transpose="logical"  # TRUE or FALSE
    )
)

get_HDF5Array_dim <- function(x) lengths(x@index)

setMethod("dim", "HDF5Array", get_HDF5Array_dim)

.get_HDF5Array_dimnames <- function(x)
{
    ans <- lapply(x@index, names)
    if (is.null(unlist(ans)))
        return(NULL)
    ans
}

setMethod("dimnames", "HDF5Array", .get_HDF5Array_dimnames)

normalize_dimnames_replacement_value <- function(value, ndim)
{
    if (is.null(value))
        return(vector("list", length=ndim))
    if (!is.list(value))
        stop("the supplied dimnames must be a list")
    if (length(value) > ndim)
        stop(wmsg("the supplied dimnames is longer ",
                  "than the number of dimensions"))
    if (length(value) <- ndim)
        length(value) <- ndim
    value
}

.set_HDF5Array_dimnames <- function(x, value)
{
    value <- normalize_dimnames_replacement_value(value, length(x@index))
    for (n in seq_along(x@index))
        ## 'x@index' can be big so avoid copies when possible. With this trick
        ## no-op dimnames(x) <- dimnames(x) is instantaneous.
        if (!identical(names(x@index[[n]]), value[[n]]))
            names(x@index[[n]]) <- value[[n]]
    x
}

setReplaceMethod("dimnames", "HDF5Array", .set_HDF5Array_dimnames)

.get_HDF5dataset_dim <- function(file, group, name)
{
    f <- H5Fopen(file, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(f))
    g <- H5Gopen(f, group)
    on.exit(H5Gclose(g), add=TRUE)
    d <- H5Dopen(g, name)
    on.exit(H5Dclose(d), add=TRUE)
    H5Sget_simple_extent_dims(H5Dget_space(d))$size
}

### Constructor.
HDF5Array <- function(file, group, name)
{
    dim0 <- .get_HDF5dataset_dim(file, group, name)
    index <- lapply(dim0, seq_len)
    new2("HDF5Array", file=file, group=group, name=name, index=index)
}

.reduce_array_dimensions <- function(x)
{
    x_dim <- dim(x)
    x_dimnames <- dimnames(x)
    edim_idx <- which(x_dim != 1L)  # index of effective dimensions
    if (length(edim_idx) >= 2L) {
        dim(x) <- x_dim[edim_idx]
        dimnames(x) <- x_dimnames[edim_idx]
    } else {
        dim(x) <- NULL
        if (length(edim_idx) == 1L)
            names(x) <- x_dimnames[[edim_idx]]
    }
    x
}

.from_HDF5Array_to_array <- function(x, drop=FALSE)
{
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    ans <- h5read(x@file, paste(x@group, x@name, sep="/"), index=x@index)
    dimnames(ans) <- .get_HDF5Array_dimnames(x)
    if (drop)
        ans <- .reduce_array_dimensions(ans)
    ans
}

setMethod("as.array", "HDF5Array", .from_HDF5Array_to_array)

slicing_tip <- c(
    "Consider reducing its number of effective dimensions by slicing it ",
    "first (e.g. x[8, 30, , 2, ]). Make sure that all the indices used for ",
    "the slicing have length 1 except at most 2 of them which can be of ",
    "arbitrary length or missing."
)

.from_HDF5Array_to_matrix <- function(x)
{
    if (sum(dim(x) != 1L) > 2L)
        stop(wmsg(class(x), " object with more than 2 effective dimensions ",
                  "cannot be coerced to a matrix. ", slicing_tip))
    ans <- as.array(x, drop=TRUE)
    as.matrix(ans)
}

setMethod("as.matrix", "HDF5Array", .from_HDF5Array_to_matrix)

setMethod("show", "HDF5Array",
    function(object)
    {
        cat(paste0(dim(object), collapse=" x "), class(object), "object\n")
    }
)

.extract_HDF5Array_subset <- function(x, i, j, ..., drop=TRUE)
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
    if (ndim != 0L && ndim != length(x@index)) {
        if (ndim == 1L)
            stop("1D-subsetting is not supported")
        stop("incorrect number of dimensions")
    }

    ## Perform the subsetting.
    if (!missing(i))
        x@index[[1L]] <- extractROWS(x@index[[1L]], i)
    if (!missing(j))
        x@index[[2L]] <- extractROWS(x@index[[2L]], j)
    ## Hack: missing values in '...' are "name" objects.
    xargs <- substitute(...())  # list of non-evaluated args
    is_missing <- sapply(xargs,
        function(xarg) { is.name(xarg) && as.character(xarg) == "" }
    )
    for (n in seq_along(xargs)) {
        if (is_missing[[n]])
            next
        xarg <- eval(xargs[[n]], envir=parent.frame(2L))
        x@index[[2L + n]] <- extractROWS(x@index[[2L + n]], xarg)
    }
    x
}

setMethod("[", "HDF5Array", .extract_HDF5Array_subset)

