### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


setClass("HDF5Array",
    representation(
        file="character",    # single string
        group="character",   # single string
        name="character",    # dataset name
        type="character",    # single string
        index="list",        # list of N integer vectors, one per dimension
        transpose="logical"  # TRUE or FALSE
    ),
    prototype(
        transpose=FALSE
    )
)

### The transposition of the data is delayed i.e. it will be realized on the
### fly only when as.array() (or as.matrix()) is called on 'x'.
setMethod("t", "HDF5Array",
    function(x)
    {
        x@transpose <- !x@transpose
        x
    }
)

get_HDF5Array_dim <- function(x)
{
    ans <- lengths(x@index)
    if (x@transpose)
        ans <- rev(ans)
    ans
}

setMethod("dim", "HDF5Array", get_HDF5Array_dim)

.get_HDF5Array_dimnames_before_transpose <- function(x)
{
    ans <- lapply(x@index, names)
    if (is.null(unlist(ans)))
        return(NULL)
    ans
}

.get_HDF5Array_dimnames <- function(x)
{
    ans <- .get_HDF5Array_dimnames_before_transpose(x)
    if (x@transpose)
        ans <- rev(ans)
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
    if (x@transpose)
        value <- rev(value)
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

.get_HDF5dataset_type <- function(file, group, name, ndim)
{
    ## Get the 1st value in the dataset. Will fail if the dataset is empty.
    ## (There must be a better way to obtain the type of an HDF5 dataset.)
    index <- rep.int(list(1L), ndim)
    typeof(h5read(file, paste(group, name, sep="/"), index=index))
}

### Constructor.
HDF5Array <- function(file, group, name)
{
    dim0 <- .get_HDF5dataset_dim(file, group, name)
    type <- .get_HDF5dataset_type(file, group, name, length(dim0))
    index <- lapply(dim0, seq_len)
    new2("HDF5Array", file=file, group=group, name=name,
                      type=type, index=index)
}

.reduce_array_dimensions <- function(x)
{
    x_dim <- dim(x)
    x_dimnames <- dimnames(x)
    effdim_idx <- which(x_dim != 1L)  # index of effective dimensions
    if (length(effdim_idx) >= 2L) {
        dim(x) <- x_dim[effdim_idx]
        dimnames(x) <- x_dimnames[effdim_idx]
    } else {
        dim(x) <- NULL
        if (length(effdim_idx) == 1L)
            names(x) <- x_dimnames[[effdim_idx]]
    }
    x
}

.from_HDF5Array_to_array <- function(x, drop=FALSE)
{
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    ans <- h5read(x@file, paste(x@group, x@name, sep="/"), index=x@index)
    dimnames(ans) <- .get_HDF5Array_dimnames_before_transpose(x)
    if (drop)
        ans <- .reduce_array_dimensions(ans)
    ## Base R doesn't support transposition of an array of arbitrary dimension
    ## (generalized transposition) so the call to t() below will fail if 'ans'
    ## has more than 2 dimensions. If we want as.array() to work on a
    ## transposed HDF5Array object of arbitrary dimension, we need to implement
    ## our own generalized transposition of an ordinary array.
    if (x@transpose) {
        if (length(dim(ans)) > 2L)
            stop("can't do as.array() on this object, sorry")
        ans <- t(ans)
    }
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

.extract_HDF5Array_subset <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")

    ## Check the dimensionality of the user call i.e whether the function was
    ## called with 1D-style, or 2D-style, or 3D-style etc... subsetting.
    ndim <- nargs() - 1L
    if (!missing(drop))
        ndim <- ndim - 1L
    if (ndim == 1L && missing(i))
        ndim <- 0L
    if (ndim != 0L && ndim != length(x@index)) {
        if (ndim == 1L)
            stop("1D-style subsetting is not supported")
        stop("incorrect number of dimensions")
    }

    ## Perform the subsetting.
    if (!missing(i)) {
        n <- if (x@transpose) length(x@index) else 1L
        x@index[[n]] <- extractROWS(x@index[[n]], i)
    }
    if (!missing(j)) {
        n <- if (x@transpose) length(x@index) - 1L else 2L
        x@index[[n]] <- extractROWS(x@index[[n]], j)
    }
    ## Hack: missing values in '...' are "name" objects.
    xargs <- substitute(...())  # list of non-evaluated args
    is_missing <- sapply(xargs,
        function(k) { is.name(k) && as.character(k) == "" }
    )
    for (n2 in seq_along(xargs)) {
        if (is_missing[[n2]])
            next
        k <- eval(xargs[[n2]], envir=parent.frame(2L))
        n <- if (x@transpose) length(x@index) - 1L - n2 else 2L + n2
        x@index[[n]] <- extractROWS(x@index[[n]], k)
    }
    x
}

setMethod("[", "HDF5Array", .extract_HDF5Array_subset)

show_HDF5Array_topline <- function(x)
{
    x_dim <- dim(x)
    cat(class(x), " object of ", paste0(x_dim, collapse=" x "),
        " ", x@type, ifelse(any(x_dim > 1L), "s", ""), sep="")
}

setMethod("show", "HDF5Array",
    function(object)
    {
        show_HDF5Array_topline(object)
        cat("\n")
    }
)

