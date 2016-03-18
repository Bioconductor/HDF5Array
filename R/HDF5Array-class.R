### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


setClass("HDF5Array",
    representation(
        file="character",           # Single string
        group="character",          # Single string
        name="character",           # Dataset name
        h5dataset_first_val="ANY",  # First value in the HDF5 dataset
        h5index="list",             # List of N integer vectors, one per
                                    # dimension in the HDF5 dataset.
        is_transposed="logical",    # Is it transposed with respect to the
                                    # HDF5 dataset?
        delayed_ops="list"          # List of delayed operations. See below for
                                    # the details.
    ),
    prototype(
        is_transposed=FALSE
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Management of delayed operations
###
### The 'delayed_ops' slot represents the list of delayed operations op1, op2,
### etc... Each delayed operation is itself represented by a list of length 4:
###   1) The name of the function to call (e.g. "+" or "log").
###   2) The list of "left arguments" i.e. the list of arguments to place
###      before the array in the function call.
###   3) The list of "right arguments" i.e. the list of arguments to place
###      after the array in the function call.
###   4) A single logical. Indicates the dimension along which the (left or
###      right) argument of the function call needs to be recycled at
###      evaluation time (evaluation is performed by .apply_delayed_ops()
###      which is called by as.array()). FALSE: along the 1st dim; TRUE: along
###      the last dim; NA: no recycling. Recycling is only supported for
###      function calls with 2 arguments (i.e. the array and the recycled
###      argument) at the moment.
###      
### Each operation must return an array of the same dimensions as the original
### array.
###

register_delayed_op <- function(x, FUN, Largs=list(), Rargs=list(),
                                        recycle_along_last_dim=NA)
{
    if (isTRUEorFALSE(recycle_along_last_dim)) {
        nLargs <- length(Largs)
        nRargs <- length(Rargs)
        ## Recycling is only supported for function calls with 2 arguments
        ## (i.e. the array and the recycled argument) at the moment.
        stopifnot(nLargs + nRargs == 1L)
        partially_recycled_arg <- if (nLargs == 1L) Largs[[1L]] else Rargs[[1L]]
        stopifnot(length(partially_recycled_arg) == nrow(x))
    }
    delayed_op <- list(FUN, Largs, Rargs, recycle_along_last_dim)
    x@delayed_ops <- c(x@delayed_ops, list(delayed_op))
    x
}

.subset_delayed_op_args <- function(delayed_op, i, subset_along_last_dim)
{
    recycle_along_last_dim <- delayed_op[[4L]]
    if (is.na(recycle_along_last_dim)
     || recycle_along_last_dim != subset_along_last_dim)
        return(delayed_op)
    Largs <- delayed_op[[2L]]
    Rargs <- delayed_op[[3L]]
    nLargs <- length(Largs)
    nRargs <- length(Rargs)
    stopifnot(nLargs + nRargs == 1L)
    if (nLargs == 1L) {
        delayed_op[[2L]] <- list(extractROWS(Largs[[1L]], i))
    } else {
        delayed_op[[3L]] <- list(extractROWS(Rargs[[1L]], i))
    }
    delayed_op
}

.subset_delayed_ops_args <- function(delayed_ops, i, subset_along_last_dim)
    lapply(delayed_ops, .subset_delayed_op_args, i, subset_along_last_dim)

### 'a' is the array returned by .read_HDF5Array_as_unprocessed_array().
.apply_delayed_ops <- function(a, delayed_ops)
{
    a_dim <- dim(a)
    first_dim <- a_dim[[1L]]
    last_dim <- a_dim[[length(a_dim)]]
    a_len <- length(a)
    if (a_len == 0L) {
        p1 <- p2 <- 0L
    } else {
        p1 <- a_len / first_dim
        p2 <- a_len / last_dim
    }

    recycle_arg <- function(partially_recycled_arg, recycle_along_last_dim) {
        if (recycle_along_last_dim) {
            stopifnot(length(partially_recycled_arg) == last_dim)
            rep(partially_recycled_arg, each=p2)
        } else {
            stopifnot(length(partially_recycled_arg) == first_dim)
            rep.int(partially_recycled_arg, p1)
        }
    }

    prepare_call_args <- function(a, delayed_op) {
        Largs <- delayed_op[[2L]]
        Rargs <- delayed_op[[3L]]
        recycle_along_last_dim <- delayed_op[[4L]]
        if (isTRUEorFALSE(recycle_along_last_dim)) {
            nLargs <- length(Largs)
            nRargs <- length(Rargs)
            stopifnot(nLargs + nRargs == 1L)
            if (nLargs == 1L) {
                Largs <- list(recycle_arg(Largs[[1L]], recycle_along_last_dim))
            } else {
                Rargs <- list(recycle_arg(Rargs[[1L]], recycle_along_last_dim))
            }
        }
        c(Largs, list(a), Rargs)
    }

    for (delayed_op in delayed_ops) {
        FUN <- delayed_op[[1L]]
        call_args <- prepare_call_args(a, delayed_op)

        ## Perform the delayed operation.
        a <- do.call(FUN, call_args)

        ## Some vectorized operations on an ordinary array can drop the dim
        ## attribute (e.g. comparing a zero-col matrix with an atomic vector).
        a_new_dim <- dim(a)
        if (is.null(a_new_dim)) {
            ## Restore the dim attribute.
            dim(a) <- a_dim
        } else {
            ## Sanity check.
            stopifnot(identical(a_dim, a_new_dim))
        }
    }
    a
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transpose
###

### The actual transposition of the data is delayed i.e. it will be realized
### on the fly only when as.array() (or as.vector() or as.matrix()) is called
### on 'x'.
setMethod("t", "HDF5Array",
    function(x)
    {
        x@is_transposed <- !x@is_transposed
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###
### None of these accessors actually needs to access the HDF5 file.
###

### The index() getter and setter are for internal use only.

setGeneric("index",                                # NOT exported
    function(x) standardGeneric("index")
)
setGeneric("index<-", signature="x",               # NOT exported
    function(x, value) standardGeneric("index<-")
)
setMethod("index", "HDF5Array",
    function(x) x@h5index
)
setReplaceMethod("index", "HDF5Array",
    function(x, value) { x@h5index <- value; x }
)

### Even though prod() always returns a double, it seems that the length()
### primitive function automatically turns this double into an integer if
### it's <= .Machine$integer.max
setMethod("length", "HDF5Array", function(x) prod(lengths(index(x))))

setMethod("isEmpty", "HDF5Array", function(x) any(lengths(index(x)) == 0L))

.get_HDF5Array_dim_before_transpose <- function(x)
{
    lengths(index(x))
}
.get_HDF5Array_dimnames_before_transpose <- function(x)
{
    ans <- lapply(index(x), names)
    if (is.null(unlist(ans)))
        return(NULL)
    ans
}

.get_HDF5Array_dim <- function(x)
{
    ans <- .get_HDF5Array_dim_before_transpose(x)
    if (x@is_transposed)
        ans <- rev(ans)
    ans
}
.get_HDF5Array_dimnames <- function(x)
{
    ans <- .get_HDF5Array_dimnames_before_transpose(x)
    if (x@is_transposed)
        ans <- rev(ans)
    ans
}

setMethod("dim", "HDF5Array", .get_HDF5Array_dim)
setMethod("dimnames", "HDF5Array", .get_HDF5Array_dimnames)

.normalize_dimnames_replacement_value <- function(value, ndim)
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
    x_index <- index(x)
    value <- .normalize_dimnames_replacement_value(value, length(x_index))
    if (x@is_transposed)
        value <- rev(value)
    index_was_touched <- FALSE
    for (n in seq_along(x_index)) {
        ## 'x_index' can be big so avoid copies when possible. With this
        ## trick no-op 'dimnames(x) <- dimnames(x)' is very cheap (i.e. almost
        ## instantaneous).
        if (identical(names(x_index[[n]]), value[[n]]))
            next
        names(x_index[[n]]) <- value[[n]]
        index_was_touched <- TRUE
    }
    if (index_was_touched)
        index(x) <- x_index
    x
}

setReplaceMethod("dimnames", "HDF5Array", .set_HDF5Array_dimnames)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers to read stuff from the HDF5 file
###

.quiet_h5read <- function(file, group, name, index)
{
    ## h5read() emits an annoying warning when it loads integer values that
    ## cannot be represented in R (and thus are converted to NAs).
    suppressWarnings(
        h5read(file, paste(group, name, sep="/"), index=index)
    )
}

### "unprocessed" here means untransposed and before applying the delayed
### operations.
.read_HDF5Array_as_unprocessed_array <- function(x)
{
    dim_before_transpose <- .get_HDF5Array_dim_before_transpose(x)
    is_empty <- any(dim_before_transpose == 0L)
    if (is_empty) {
        ans <- x@h5dataset_first_val[0]
    } else {
        ans <- .quiet_h5read(x@file, x@group, x@name, x@h5index)
    }
    dim(ans) <- dim_before_transpose
    dimnames(ans) <- .get_HDF5Array_dimnames_before_transpose(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.get_h5dataset_dim <- function(file, group, name)
{
    f <- H5Fopen(file, flags="H5F_ACC_RDONLY")
    on.exit(H5Fclose(f))
    g <- H5Gopen(f, group)
    on.exit(H5Gclose(g), add=TRUE)
    d <- H5Dopen(g, name)
    on.exit(H5Dclose(d), add=TRUE)
    H5Sget_simple_extent_dims(H5Dget_space(d))$size
}

### Will fail if the dataset is empty (i.e. if at least one of its dimensions
### is 0).
.read_h5dataset_first_val <- function(file, group, name, ndim)
{
    index <- rep.int(list(1L), ndim)
    ans <- .quiet_h5read(file, group, name, index)
    stopifnot(length(ans) == 1L)  # sanity check
    ans[[1L]]  # drop any attribute
}

HDF5Array <- function(file, group, name, type=NA)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the HDF5 file where the dataset is located"))
    if (!isSingleString(group))
        stop(wmsg("'group' must be a single string specifying the HDF5 group ",
                  "of the dataset, as reported by rhdf5::h5ls()"))
    if (!isSingleString(name))
        stop(wmsg("'name' must be a single string specifying the name ",
                  "of the dataset, as reported by rhdf5::h5ls()"))
    if (!isSingleStringOrNA(type))
        stop("'type' must be a single string or NA")
    h5dataset_dim <- .get_h5dataset_dim(file, group, name)
    if (any(h5dataset_dim == 0L)) {
        if (is.na(type))
            stop(wmsg("This HDF5 dataset is empty! Don't know how to ",
                      "determine the type of an empty HDF5 dataset at the ",
                      "moment. Please use the 'type' argument to help me ",
                      "(see '?HDF5Array' for more information)."))
        h5dataset_first_val <- match.fun(type)(1)  # fake value
        if (!is.atomic(h5dataset_first_val))
            stop("invalid type: ", type)
    } else {
        h5dataset_first_val <- .read_h5dataset_first_val(file, group, name,
                                                         length(h5dataset_dim))
        detected_type <- typeof(h5dataset_first_val)
        if (!(is.na(type) || type == detected_type))
            warning(wmsg("The type specified via the 'type' argument (",
                         type, ") doesn't match the type of this HDF5 ",
                         "dataset (", detected_type, "). Ignoring the ",
                         "former."))
    }
    h5index <- lapply(h5dataset_dim, seq_len)
    new2("HDF5Array", file=file,
                      group=group,
                      name=name,
                      h5dataset_first_val=h5dataset_first_val,
                      h5index=h5index)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

### HDF5Array -> array

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
    ans <- .read_HDF5Array_as_unprocessed_array(x)
    ans <- .apply_delayed_ops(ans, x@delayed_ops)
    if (drop)
        ans <- .reduce_array_dimensions(ans)
    ## Base R doesn't support transposition of an array of arbitrary dimension
    ## (generalized transposition) so the call to t() below will fail if 'ans'
    ## has more than 2 dimensions. If we want as.array() to work on a
    ## transposed HDF5Array object of arbitrary dimension, we need to implement
    ## our own generalized transposition of an ordinary array.
    if (x@is_transposed) {
        if (length(dim(ans)) > 2L)
            stop("can't do as.array() on this object, sorry")
        ans <- t(ans)
    }
    ans
}

setMethod("as.array", "HDF5Array", .from_HDF5Array_to_array)

### HDF5Array -> vector

.from_HDF5Array_to_vector <- function(x, mode="any")
{
    ans <- as.array(x, drop=TRUE)
    as.vector(ans, mode=mode)
}

setMethod("as.vector", "HDF5Array", .from_HDF5Array_to_vector)

### HDF5Array -> matrix

slicing_tip <- c(
    "Consider reducing its number of effective dimensions by slicing it ",
    "first (e.g. x[8, 30, , 2, ]). Make sure that all the indices used for ",
    "the slicing have length 1 except at most 2 of them which can be of ",
    "arbitrary length or missing."
)

.from_HDF5Array_to_matrix <- function(x)
{
    x_dim <- dim(x)
    if (sum(x_dim != 1L) > 2L)
        stop(wmsg(class(x), " object with more than 2 effective dimensions ",
                  "cannot be coerced to a matrix. ", slicing_tip))
    ans <- as.array(x, drop=TRUE)
    if (length(x_dim) == 2L) {
        dim(ans) <- x_dim
        dimnames(ans) <- dimnames(x)
    } else {
        as.matrix(ans)
    }
    ans
}

setMethod("as.matrix", "HDF5Array", .from_HDF5Array_to_matrix)

### array -> HDF5Array

### Write a new HDF5 dataset to disk. Return an HDF5Array object that points
### to this new dataset.
.from_array_to_HDF5Array <- function(from)
{
    out_file <- getHDF5ArrayOutputFile()
    out_name <- getHDF5ArrayOutputName()
    on.exit(setHDF5ArrayOutputName())

    h5write(from, out_file, out_name)
    ## TODO: Investigate the possiblity to store the dimnames in the HDF5 file
    ## so the HDF5Array() constructor can bring them back. Then we wouldn't
    ## need to explicitely set them on 'ans' like we do below.
    ans <- HDF5Array(out_file, "/", out_name, type=type(from))
    dimnames(ans) <- dimnames(from)
    ans
}

setAs("array", "HDF5Array", .from_array_to_HDF5Array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

### 'subscript' must be a multidimensional subscript i.e. a list with one
### subscript per dimension in 'x'. Missing subscripts are represented by
### list elements of class "name".
.extract_subarray_from_HDF5Array <- function(x, subscript)
{
    x_index <- index(x)
    x_ndim <- length(x_index)
    x_delayed_ops <- x@delayed_ops
    index_was_touched <- FALSE
    for (n in seq_along(subscript)) {
        h5n <- if (x@is_transposed) x_ndim - n + 1L else n
        k <- subscript[[n]]
        if (missing(k))
            next
        x_index[[h5n]] <- extractROWS(x_index[[h5n]], k)
        index_was_touched <- TRUE
        if (h5n == 1L)
            x_delayed_ops <- .subset_delayed_ops_args(x_delayed_ops, k, FALSE)
        if (h5n == x_ndim)
            x_delayed_ops <- .subset_delayed_ops_args(x_delayed_ops, k, TRUE)
    }
    if (index_was_touched) {
        index(x) <- x_index
        if (!identical(x@delayed_ops, x_delayed_ops))
            x@delayed_ops <- x_delayed_ops
    }
    x
}

.extract_subarray <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")

    ## Check the dimensionality of the user call i.e whether the function was
    ## called with 1D-style, or 2D-style, or 3D-style etc... subsetting.
    ndim <- nargs() - 1L
    x_ndim <- length(dim(x))
    if (!missing(drop))
        ndim <- ndim - 1L
    if (ndim == 1L && missing(i))
        ndim <- 0L
    if (ndim != 0L && ndim != x_ndim) {
        if (ndim == 1L)
            stop("1D-style subsetting is not supported")
        stop("incorrect number of dimensions")
    }

    ## Prepare the multidimensional subscript.
    subscript <- rep.int(alist(foo=), x_ndim)
    if (!missing(i))
        subscript[[1L]] <- i
    if (!missing(j))
        subscript[[2L]] <- j
    dots <- substitute(...())  # list of non-evaluated args
    for (n2 in seq_along(dots)) {
        k <- dots[[n2]]
        if (!missing(k))
            subscript[[2L + n2]] <- eval(k, envir=parent.frame(2L))
    }

    ## Performs the subsetting.
    .extract_subarray_from_HDF5Array(x, subscript)
}

setMethod("[", "HDF5Array", .extract_subarray)

.get_HDF5Array_element <- function(x, i)
{
    i <- S4Vectors:::normalizeDoubleBracketSubscript(i, x)
    subscript <- as.integer(arrayInd(i, dim(x)))
    as.vector(.extract_subarray_from_HDF5Array(x, subscript))
}

### Only support linear subscripting at the moment.
### TODO: Support multidimensional subscripting e.g. x[[1, 5]].
setMethod("[[", "HDF5Array",
    function(x, i, j, ...)
    {
        dots <- list(...)
        if (length(dots) > 0L)
            dots <- dots[names(dots) != "exact"]
        if (!missing(j) || length(dots) > 0L)
            stop("incorrect number of subscripts")
        .get_HDF5Array_element(x, i)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type()
###

### If 'x' is an HDF5Array object, 'type(x)' must always return the same
### as 'typeof(as.array(x))'. For internal use only.
setMethod("type", "HDF5Array",
    function(x)
    {
        subscript <- as.list(integer(length(dim(x))))
        x0 <- .extract_subarray_from_HDF5Array(x, subscript)  # x[0, ..., 0]
        typeof(as.array(x0, drop=TRUE))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

show_HDF5Array_topline <- function(x)
{
    x_dim <- dim(x)
    cat(class(x), " object of ", paste0(dim(x), collapse=" x "),
        " ", type(x), ifelse(any(x_dim >= 2L), "s", ""), sep="")
}

setMethod("show", "HDF5Array",
    function(object)
    {
        show_HDF5Array_topline(object)
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###
### Combining arrays with c() is NOT an endomorphism!
###

### 'objects' must be a list of array-like objects that support as.vector().
combine_array_objects <- function(objects)
{
    if (!is.list(objects))
        stop("'objects' must be a list")

    NULL_idx <- which(S4Vectors:::sapply_isNULL(objects))
    if (length(NULL_idx) != 0L)
        objects <- objects[-NULL_idx]
    if (length(objects) == 0L)
        return(NULL)

    unlist(lapply(objects, as.vector), recursive=FALSE, use.names=FALSE)
}

setMethod("c", "HDF5Array",
    function (x, ..., recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("\"c\" method for HDF5Array objects ",
                 "does not support the 'recursive' argument")
        if (missing(x)) {
            objects <- list(...)
        } else {
            objects <- list(x, ...)
        }
        combine_array_objects(objects)
    }
)

