### =========================================================================
### DelayedArray objects
### -------------------------------------------------------------------------


setClass("DelayedArray",
    representation(
        seeds="list",              # List of n conformable array-like objects.
        index="list",              # List of (possibly named) integer vectors,
                                   # one per dimension in the seeds.
        COMBINING_OP="character",  # n-ary operator to combine the seeds
                                   # after they all went thru
                                   # extract_array_from_seed(seed, index).
        Rargs="list",
        delayed_ops="list",        # List of delayed operations. See below
                                   # for the details.
        is_transposed="logical"    # Is the object considered to be transposed
                                   # with respect to the seeds?
    ),
    prototype(
        seeds=list(new("array")),
        index=list(integer(0)),
        COMBINING_OP="identity",
        is_transposed=FALSE
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.objects_are_conformable_arrays <- function(objects)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    first_ndim <- ndims[[1L]]
    if (!all(ndims == first_ndim))
        return(FALSE)
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        return(FALSE)
    dims <- matrix(tmp, nrow=first_ndim)
    first_dim <- dims[ , 1L]
    all(dims == first_dim)
}

.validate_DelayedArray <- function(x)
{
    ## 'seeds' slot.
    if (length(x@seeds) == 0L)
        return(wmsg("'x@seeds' cannot be empty"))
    if (!.objects_are_conformable_arrays(x@seeds))
        return(wmsg("'x@seeds' must be a list of conformable ",
                    "array-like objects"))
    ## 'index' slot.
    if (length(x@index) != length(dim(x@seeds[[1L]])))
        return(wmsg("'x@index' must have one element per dimension ",
                    "in 'x@seeds[[1L]]'"))
    ## 'COMBINING_OP' slot.
    if (!isSingleString(x@COMBINING_OP))
        return(wmsg("'x@COMBINING_OP' must be a single string"))
    OP <- try(match.fun(x@COMBINING_OP), silent=TRUE)
    if (is(OP, "try-error"))
        return(wmsg("the name in 'x@COMBINING_OP' (\"", x@COMBINING_OP,
                    "\") must refer to a known function"))
    ## 'is_transposed' slot.
    if (!isTRUEorFALSE(x@is_transposed))
        return(wmsg("'x@is_transposed' must be TRUE or FALSE"))
    TRUE
}

setValidity2("DelayedArray", .validate_DelayedArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### For internal use only.
new_DelayedArray <- function(a=new("array"),
                             ..., COMBINING_OP="identity", Rargs=list(),
                             Class="DelayedArray")
{
    seeds <- list(a, ...)
    index <- lapply(dim(a), seq_len)
    new2(Class, seeds=seeds, COMBINING_OP=COMBINING_OP, Rargs=Rargs,
                index=index)
}

setAs("ANY", "DelayedArray",
    function(from)
    {
        ans <- new_DelayedArray(from)
        dimnames(ans) <- dimnames(from)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pristine objects
###
### A pristine DelayedArray object is an object that does not carry any
### delayed operation on it. In other words, it's in sync with (i.e. reflects
### the true content of) its unique seed.
###

### Note that false negatives happen when 'x' carries delayed operations that
### do nothing, but that's ok.
is_pristine <- function(x)
{
    if (length(x@seeds) != 1L)
        return(FALSE)
    ## 'x' should not carry any delayed operation on it, that is, all the
    ## DelayedArray slots must be in their original state.
    x1 <- new_DelayedArray(x@seeds[[1L]])
    x2 <- as(x, "DelayedArray", strict=TRUE)
    dimnames(x2) <- NULL
    if (!identical(x1, x2))
        return(FALSE)
    if (!is(x, "DelayedMatrix"))
        return(TRUE)
    length(x@index) == 2L && x@N1 == 1L && x@N2 == 2L
}

### When a pristine DelayedArray derived object (i.e. an HDF5Array object) is
### about to be touched, we first need to downgrade it to a DelayedArray or
### DelayedMatrix *instance*.
downgrade_to_DelayedArray_or_DelayedMatrix <- function(x)
{
    if (is(x, "DelayedMatrix"))
        return(as(x, "DelayedMatrix", strict=TRUE))
    as(x, "DelayedArray", strict=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

### The index() getter and setter are for internal use only.

setGeneric("index",                                # NOT exported
    function(x) standardGeneric("index")
)
setGeneric("index<-", signature="x",               # NOT exported
    function(x, value) standardGeneric("index<-")
)
setMethod("index", "DelayedArray",
    function(x) x@index
)
setReplaceMethod("index", "DelayedArray",
    function(x, value) { x@index <- value; x }
)

.get_DelayedArray_dim_before_transpose <- function(x)
{
    lengths(index(x))
}
.get_DelayedArray_dim <- function(x)
{
    ans <- .get_DelayedArray_dim_before_transpose(x)
    if (x@is_transposed)
        ans <- rev(ans)
    ans
}

setMethod("dim", "DelayedArray", .get_DelayedArray_dim)

### Even though prod() always returns a double, it seems that the length()
### primitive function automatically turns this double into an integer if
### it's <= .Machine$integer.max
setMethod("length", "DelayedArray", function(x) prod(dim(x)))

setMethod("isEmpty", "DelayedArray", function(x) any(dim(x) == 0L))

.get_DelayedArray_dimnames_before_transpose <- function(x)
{
    ans <- lapply(index(x), names)
    if (is.null(unlist(ans)))
        return(NULL)
    ans
}
.get_DelayedArray_dimnames <- function(x)
{
    ans <- .get_DelayedArray_dimnames_before_transpose(x)
    if (x@is_transposed)
        ans <- rev(ans)
    ans
}

setMethod("dimnames", "DelayedArray", .get_DelayedArray_dimnames)

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

.set_DelayedArray_dimnames <- function(x, value)
{
    x_index <- index(x)
    x_ndim <- length(x_index)
    value <- .normalize_dimnames_replacement_value(value, x_ndim)
    if (x@is_transposed)
        value <- rev(value)
    x_index_was_touched <- FALSE
    for (n in seq_along(x_index)) {
        ## 'x_index' can be big so avoid copies when possible. With this
        ## trick no-op 'dimnames(x) <- dimnames(x)' is very cheap (i.e. almost
        ## instantaneous).
        if (identical(names(x_index[[n]]), value[[n]]))
            next
        names(x_index[[n]]) <- value[[n]]
        x_index_was_touched <- TRUE
    }
    if (x_index_was_touched)
        index(x) <- x_index
    x
}

setReplaceMethod("dimnames", "DelayedArray", .set_DelayedArray_dimnames)

.get_DelayedArray_names <- function(x)
{
    if (length(dim(x)) != 1L)
        return(NULL)
    dimnames(x)[[1L]]
}

setMethod("names", "DelayedArray", .get_DelayedArray_names)

.set_DelayedArray_names <- function(x, value)
{
    if (length(dim(x)) != 1L) {
        if (!is.null(value))
            stop("setting the names of a ", class(x), " object with more ",
                 "than 1 dimension is not supported")
        return(x)
    }
    dimnames(x)[[1L]] <- value
    x
}

setReplaceMethod("names", "DelayedArray", .set_DelayedArray_names)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [
###

### 'subscript' must be a multidimensional subscript i.e. a list with one
### subscript per dimension in 'x'. Missing subscripts are represented by
### list elements of class "name".
.extract_subarray_from_DelayedArray <- function(x, subscript)
{
    x_index <- x_index0 <- index(x)
    x_ndim <- length(x_index)
    x_delayed_ops <- x@delayed_ops
    for (n in seq_along(subscript)) {
        n0 <- if (x@is_transposed) x_ndim - n + 1L else n
        k <- subscript[[n]]
        if (missing(k))
            next
        x_index[[n0]] <- extractROWS(x_index[[n0]], k)
        if (n0 == 1L)
            x_delayed_ops <- .subset_delayed_ops_args(x_delayed_ops, k, FALSE)
        if (n0 == x_ndim)
            x_delayed_ops <- .subset_delayed_ops_args(x_delayed_ops, k, TRUE)
    }
    if (!identical(x_index0, x_index)) {
        x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
        index(x) <- x_index
        if (!identical(x@delayed_ops, x_delayed_ops))
            x@delayed_ops <- x_delayed_ops
    }
    x
}

.extract_DelayedArray_subset <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop("'x' is missing")

    ## Check the dimensionality of the user call i.e whether the function was
    ## called with 1D-style, or 2D-style, or 3D-style etc... subsetting.
    ndim <- nargs() - 1L
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
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

    ## Perform the subsetting.
    .extract_subarray_from_DelayedArray(x, subscript)
}

setMethod("[", "DelayedArray", .extract_DelayedArray_subset)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array_from_seed()
###
### 'index' is expected to be a list of integer vectors. There is no
### need to support anything else.
### Must return an ordinary array. No need to propagate the dimnames.
###

setGeneric("extract_array_from_seed", signature="seed",
    function(seed, index) standardGeneric("extract_array_from_seed")
)

setMethod("extract_array_from_seed", "ANY",
    function(seed, index)
        as.array(do.call(`[`, c(list(seed), index, drop=FALSE)))
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
    x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
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

### 'a' is the ordinary array returned by the "combining" operator.
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
setMethod("t", "DelayedArray",
    function(x)
    {
        x <- downgrade_to_DelayedArray_or_DelayedMatrix(x)
        x@is_transposed <- !x@is_transposed
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### as.array()
###

### TODO: Not sure we need this. Using drop() should do it.
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

.from_DelayedArray_to_array <- function(x, drop=FALSE)
{
    if (!isTRUEorFALSE(drop))
        stop("'drop' must be TRUE or FALSE")
    dim_before_transpose <- .get_DelayedArray_dim_before_transpose(x)
    arrays <- lapply(x@seeds,
        function(seed) {
            a <- extract_array_from_seed(seed, x@index)
            dim(a) <- dim_before_transpose
            a
        })
    ans <- do.call(x@COMBINING_OP, c(arrays, x@Rargs))
    ans <- .apply_delayed_ops(ans, x@delayed_ops)
    dimnames(ans) <- .get_DelayedArray_dimnames_before_transpose(x)
    if (drop)
        ans <- .reduce_array_dimensions(ans)
    ## Base R doesn't support transposition of an array of arbitrary dimension
    ## (generalized transposition) so the call to t() below will fail if 'ans'
    ## has more than 2 dimensions. If we want as.array() to work on a
    ## transposed DelayedArray object of arbitrary dimension, we need to
    ## implement our own generalized transposition of an ordinary array.
    if (x@is_transposed) {
        if (length(dim(ans)) > 2L)
            stop("can't do as.array() on this object, sorry")
        ans <- t(ans)
    }
    ans
}

### S3/S4 combo for as.array.DelayedArray
as.array.DelayedArray <- function(x, ...)
    .from_DelayedArray_to_array(x, ...)
setMethod("as.array", "DelayedArray", .from_DelayedArray_to_array)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other coercions based on as.array()
###

### S3/S4 combo for as.vector.DelayedArray
as.vector.DelayedArray <- function(x, mode="any")
{
    ans <- as.array(x, drop=TRUE)
    as.vector(ans, mode=mode)
}
setMethod("as.vector", "DelayedArray", as.vector.DelayedArray)

slicing_tip <- c(
    "Consider reducing its number of effective dimensions by slicing it ",
    "first (e.g. x[8, 30, , 2, ]). Make sure that all the indices used for ",
    "the slicing have length 1 except at most 2 of them which can be of ",
    "arbitrary length or missing."
)

.from_DelayedArray_to_matrix <- function(x)
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

### S3/S4 combo for as.matrix.DelayedArray
as.matrix.DelayedArray <- function(x, ...) .from_DelayedArray_to_matrix(x, ...)
setMethod("as.matrix", "DelayedArray", .from_DelayedArray_to_matrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type()
###
### For internal use only.
###

setGeneric("type", function(x) standardGeneric("type"))

setMethod("type", "array", function(x) typeof(x))

### If 'x' is a DelayedArray object, 'type(x)' must always return the same
### as 'typeof(as.array(x))'.
setMethod("type", "DelayedArray",
    function(x)
    {
        subscript <- as.list(integer(length(dim(x))))
        ## x0 <- x[0, ..., 0]
        x0 <- .extract_subarray_from_DelayedArray(x, subscript)
        typeof(as.array(x0, drop=TRUE))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### [[
###

.get_DelayedArray_element <- function(x, i)
{
    i <- S4Vectors:::normalizeDoubleBracketSubscript(i, x)
    subscript <- as.integer(arrayInd(i, dim(x)))
    as.vector(.extract_subarray_from_DelayedArray(x, subscript))
}

### Only support linear subscripting at the moment.
### TODO: Support multidimensional subscripting e.g. x[[5, 15, 2]] or
### x[["E", 15, "b"]].
setMethod("[[", "DelayedArray",
    function(x, i, j, ...)
    {
        dots <- list(...)
        if (length(dots) > 0L)
            dots <- dots[names(dots) != "exact"]
        if (!missing(j) || length(dots) > 0L)
            stop("incorrect number of subscripts")
        .get_DelayedArray_element(x, i)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

.split_names <- function(x_names, idx1, idx2)
{
    make_elt_indices <- function(i) {
        if (length(i) == 0L)
            return(character(0))
        paste0("[", i, "]", sep="")
    }
    if (is.null(x_names)) {
        s1 <- make_elt_indices(idx1)
        s2 <- make_elt_indices(idx2)
    } else {
        s1 <- x_names[idx1]
        s2 <- x_names[idx2]
    }
    format(c(s1, ".", s2), justify="right")
}

.prepare_1D_sample <- function(x, n1, n2)
{
    x_len <- length(x)
    x_names <- names(x)
    if (x_len <= n1 + n2 + 1L) {
        ans <- format(as.vector(x))
        idx1 <- seq_len(x_len)
        idx2 <- integer(0)
        names(ans) <- .split_names(x_names, idx1, idx2)[idx1]
    } else {
        idx1 <- seq_len(n1)
        idx2 <- seq(to=x_len, by=1L, length.out=n2)
        ans1 <- format(as.vector(x[idx1]))
        ans2 <- format(as.vector(x[idx2]))
        ans <- c(ans1, ".", ans2)
        names(ans) <- .split_names(x_names, idx1, idx2)
    }
    ans
}

.print_1D_sample <- function(x, n1, n2)
{
    stopifnot(length(dim(x)) == 1L)
    out <- .prepare_1D_sample(x, n1, n2)
    print(out, quote=FALSE, right=TRUE, max=length(out))
}

.split_rownames <- function(x_rownames, idx1, idx2)
{
    make_row_indices <- function(i) {
        if (length(i) == 0L)
            return(character(0))
        paste0("[", i, ",]", sep="")
    }
    if (is.null(x_rownames)) {
       s1 <- make_row_indices(idx1)
       s2 <- make_row_indices(idx2)
    } else {
       s1 <- x_rownames[idx1]
       s2 <- x_rownames[idx2]
    }
    max_width <- max(nchar(s1, type="width"), nchar(s2, type="width"))
    if (max_width <= 1L) {
        ellipsis <- "."
    } else if (max_width == 2L) {
        ellipsis <- ".."
    } else {
        ellipsis <- "..."
    }
    format(c(s1, ellipsis, s2), justify="right")
}

.split_colnames <- function(x_colnames, idx1, idx2)
{
    make_col_indices <- function(j) {
        if (length(j) == 0L)
            return(character(0))
        paste0("[,", j, "]", sep="")
    }
    if (is.null(x_colnames)) {
        s1 <- make_col_indices(idx1)
        s2 <- make_col_indices(idx2)
    } else {
        s1 <- x_colnames[idx1]
        s2 <- x_colnames[idx2]
    }
    format(c(s1, ".", s2), justify="right")
}

.rsplit_2D_DelayedArray <- function(x, m1, m2)
{
    x_nrow <- nrow(x)
    x_rownames <- rownames(x)
    idx1 <- seq_len(m1)
    idx2 <- seq(to=x_nrow, by=1L, length.out=m2)

    ans1 <- format(as.matrix(x[idx1, , drop=FALSE]))
    ans2 <- format(as.matrix(x[idx2, , drop=FALSE]))
    dots <- rep.int(".", ncol(ans1))
    ans <- rbind(ans1, matrix(dots, nrow=1L), ans2)

    rownames(ans) <- .split_rownames(x_rownames, idx1, idx2)
    ans
}

.csplit_2D_DelayedArray <- function(x, n1, n2)
{
    x_ncol <- ncol(x)
    x_colnames <- colnames(x)
    idx1 <- seq_len(n1)
    idx2 <- seq(to=x_ncol, by=1L, length.out=n2)

    ans1 <- format(as.matrix(x[ , idx1, drop=FALSE]))
    ans2 <- format(as.matrix(x[ , idx2, drop=FALSE]))
    dots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(dots, ncol=1L), ans2)

    colnames(ans) <- .split_colnames(x_colnames, idx1, idx2)
    ans
}

.split_2D_DelayedArray <- function(x, m1, m2, n1, n2)
{
    x_ncol <- ncol(x)
    x_colnames <- colnames(x)
    idx1 <- seq_len(n1)
    idx2 <- seq(to=x_ncol, by=1L, length.out=n2)

    x1 <- x[ , idx1, drop=FALSE]
    x2 <- x[ , idx2, drop=FALSE]
    ans1 <- .rsplit_2D_DelayedArray(x1, m1, m2)
    ans2 <- .rsplit_2D_DelayedArray(x2, m1, m2)
    dots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(dots, ncol=1L), ans2)

    colnames(ans) <- .split_colnames(x_colnames, idx1, idx2)
    ans
}

.prepare_2D_sample <- function(x, m1, m2, n1, n2)
{
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_nrow <= m1 + m2 + 1L) {
        if (x_ncol <= n1 + n2 + 1L) {
            ans <- format(as.matrix(x))
            ## Only needed because of this bug in base::print.default:
            ##   https://stat.ethz.ch/pipermail/r-devel/2016-March/072479.html
            ## TODO: Remove when the bug is fixed.
            if (is.null(colnames(ans))) {
                idx1 <- seq_len(ncol(ans))
                idx2 <- integer(0)
                colnames(ans) <- .split_colnames(NULL, idx1, idx2)[idx1]
            }
        } else {
            ans <- .csplit_2D_DelayedArray(x, n1, n2)
        }
    } else {
        if (x_ncol <= n1 + n2 + 1L) {
            ans <- .rsplit_2D_DelayedArray(x, m1, m2)
            ## Only needed because of this bug in base::print.default:
            ##   https://stat.ethz.ch/pipermail/r-devel/2016-March/072479.html
            ## TODO: Remove when the bug is fixed.
            if (is.null(colnames(ans))) {
                idx1 <- seq_len(ncol(ans))
                idx2 <- integer(0)
                colnames(ans) <- .split_colnames(NULL, idx1, idx2)[idx1]
            }
        } else {
            ans <- .split_2D_DelayedArray(x, m1, m2, n1, n2)
        }
    }
    ans
}

.print_2D_sample <- function(x, m1, m2, n1, n2)
{
    stopifnot(length(dim(x)) == 2L)
    out <- .prepare_2D_sample(x, m1, m2, n1, n2)
    print(out, quote=FALSE, right=TRUE, max=length(out))
}

.print_2D_slices <- function(x, m1, m2, n1, n2, blocks, idx)
{
    subscript2string <- function(subscript, dimnames=NULL) {
        s <- as.character(subscript)
        if (!is.null(dimnames)) {
            usename_idx <- which(nzchar(s) &
                                 lengths(subscript) == 1L &
                                 lengths(dimnames) != 0L)
            s[usename_idx] <- mapply(`[`, dimnames[usename_idx],
                                     subscript[usename_idx])
        }
        paste0(s, collapse=", ")
    }
    x_dimnames <- dimnames(x)
    for (i in idx) {
        subscript <- get_array_block_subscript(blocks, i,
                                               expand.RangeNSBS=TRUE)
        cat(subscript2string(subscript, x_dimnames), "\n", sep="")
        slice <- as(extract_array_block1(x, subscript), "DelayedMatrix")
        .print_2D_sample(slice, m1, m2, n1, n2)
        cat("\n")
    }
}

.print_array_data <- function(x, n1, n2)
{
    x_dim <- dim(x)
    if (length(x_dim) == 1L)
        return(.print_1D_sample(x, n1, n2))
    if (length(x_dim) == 2L) {
        nhead <- get_showHeadLines()
        ntail <- get_showTailLines()
        return(.print_2D_sample(x, nhead, ntail, n1, n2))
    }
    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]
    if (x_ncol <= 5L) {
        if (x_nrow <= 3L) {
            m1 <- m2 <- 3L  # print all rows of each slice
            z1 <- z2 <- 3L  # print first 3 and last 3 slices
        } else {
            m1 <- m2 <- 2L  # print first 2 and last 2 rows of each slice
            z1 <- z2 <- 1L  # print only first and last slices
        }
    } else {
        if (x_nrow <= 3L) {
            m1 <- m2 <- 2L  # print first 2 and last 2 rows of each slice
            z1 <- z2 <- 2L  # print first 2 and last 2 slices
        } else {
            m1 <- m2 <- 2L  # print first 2 and last 2 rows of each slice
            z1 <- z2 <- 1L  # print only first and last slices
        }
    }
    blocks <- ArrayBlocks(x_dim, prod(x_dim[1:2]))
    nblock <- length(blocks)
    if (nblock <= z1 + z2 + 1L) {
        idx <- seq_len(nblock)
        .print_2D_slices(x, m1, m2, n1, n2, blocks, idx)
    } else {
        idx1 <- seq_len(z1)
        idx2 <- seq(to=nblock, by=1L, length.out=z2)
        .print_2D_slices(x, m1, m2, n1, n2, blocks, idx1)
        cat("...\n\n")
        .print_2D_slices(x, m1, m2, n2, n2, blocks, idx2)
    }
}

setMethod("show", "DelayedArray",
    function(object)
    {
        object_class <- class(object)
        object_dim <- dim(object)
        dim_in1string <- paste0(object_dim, collapse=" x ")
        object_type <- type(object)
        if (any(object_dim == 0L)) {
            cat(sprintf("<%s> %s object of type \"%s\"\n",
                        dim_in1string, object_class, object_type))
        } else {
            cat(sprintf("%s object of %s %s%s:\n",
                        object_class, dim_in1string, object_type,
                        ifelse(any(object_dim >= 2L), "s", "")))
            .print_array_data(object, 4L, 4L)
        }
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

setMethod("c", "DelayedArray",
    function (x, ..., recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("\"c\" method for DelayedArray objects ",
                 "does not support the 'recursive' argument")
        if (missing(x)) {
            objects <- list(...)
        } else {
            objects <- list(x, ...)
        }
        combine_array_objects(objects)
    }
)

