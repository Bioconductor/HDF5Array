### =========================================================================
### Common operations on DelayedArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various unary operators + the "Math" group generic
###
### All these operations return a DelayedArray object of the same dimensions
### as 'x'.
###

setMethod("is.na", "DelayedArray", function(x) register_delayed_op(x, "is.na"))

setMethod("!", "DelayedArray", function(x) register_delayed_op(x, "!"))

setMethod("Math", "DelayedArray", function(x) register_delayed_op(x, .Generic))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Ops" group generics
###
### Arith members: "+", "-", "*", "/", "^", "%%", "%/%"
### Compare members: ==, !=, <=, >=, <, >
### Logic members: &, |
###

### Return a DelayedArray object of the same dimensions as 'e1'.
.DelayedArray_Ops_with_right_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e1, "DelayedArray"))
    e1_class <- class(e1)
    e2_class <- class(e2)
    if (!is.vector(e2))
        e2 <- as.vector(e2)
    if (!is.atomic(e2))
        stop(wmsg("`", .Generic, "` between ", e1_class, " and ",
                  e2_class, " objects is not supported"))
    e2_len <- length(e2)
    if (e2_len == 1L)
        return(register_delayed_op(e1, .Generic, Rargs=list(e2)))
    e1_len <- length(e1)
    if (e2_len > e1_len)
        stop(wmsg("right object is longer than left object"))
    e1_nrow <- nrow(e1)
    if (e1_nrow != 0L) {
        if (e2_len == 0L || e1_nrow %% e2_len != 0L)
            stop(wmsg("length of right object is not a divisor ",
                      "of number of rows of left object"))
        e2 <- rep(e2, length.out=nrow(e1))
    }
    register_delayed_op(e1, .Generic, Rargs=list(e2),
                                      recycle_along_last_dim=e1@is_transposed)
}

### Return a DelayedArray object of the same dimensions as 'e2'.
.DelayedArray_Ops_with_left_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e2, "DelayedArray"))
    e1_class <- class(e1)
    e2_class <- class(e2)
    if (!is.vector(e1))
        e1 <- as.vector(e1)
    if (!is.atomic(e1))
        stop(wmsg("`", .Generic, "` between ", e1_class, " and ",
                  e2_class, " objects is not supported"))
    e1_len <- length(e1)
    if (e1_len == 1L)
        return(register_delayed_op(e2, .Generic, Largs=list(e1)))
    e2_len <- length(e2)
    if (e1_len > e2_len)
        stop(wmsg("left object is longer than right object"))
    e2_nrow <- nrow(e2)
    if (e2_nrow != 0L) {
        if (e1_len == 0L || e2_nrow %% e1_len != 0L)
            stop(wmsg("length of left object is not a divisor ",
                      "of number of rows of right object"))
        e1 <- rep(e1, length.out=nrow(e2))
    }
    register_delayed_op(e2, .Generic, Largs=list(e1),
                                      recycle_along_last_dim=e2@is_transposed)
}

.combine_dimnames <- function(e1, e2)
{
    ans_rownames <- rownames(e1)
    if (is.null(ans_rownames))
        ans_rownames <- rownames(e2)
    ans_colnames <- colnames(e1)
    if (is.null(ans_colnames))
        ans_colnames <- colnames(e2)
    if (is.null(ans_rownames) && is.null(ans_colnames)) {
        ans_dimnames <- NULL
    } else {
        ans_dimnames <- list(ans_rownames, ans_colnames)
    }
    ans_dimnames
}

### Return a DelayedArray object of the same dimensions as 'e1' and 'e2'.
.DelayedArray_Ops_COMBINE_seeds <- function(.Generic, e1, e2)
{
    if (!identical(dim(e1), dim(e2)))
        stop("non-conformable arrays")
    ans <- new_DelayedArray(e1, e2, COMBINING_OP=.Generic)
    dimnames(ans) <- .combine_dimnames(e1, e2)
    ans
}

.DelayedArray_Ops <- function(.Generic, e1, e2)
{
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (identical(e1_dim, e2_dim))
        return(.DelayedArray_Ops_COMBINE_seeds(.Generic, e1, e2))
    ## Effective dimensions.
    effdim_idx1 <- which(e1_dim != 1L)
    effdim_idx2 <- which(e2_dim != 1L)
    if ((length(effdim_idx1) == 1L) == (length(effdim_idx2) == 1L))
        stop("non-conformable arrays")
    if (length(effdim_idx1) == 1L) {
        .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
    } else {
        .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
    }
}

setMethod("Ops", c("DelayedArray", "vector"),
    function(e1, e2)
        .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
)

setMethod("Ops", c("vector", "DelayedArray"),
    function(e1, e2)
        .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
)

setMethod("Ops", c("DelayedArray", "DelayedArray"),
    function(e1, e2)
        .DelayedArray_Ops(.Generic, e1, e2)
)

### Support unary operators "+" and "-".
setMethod("+", c("DelayedArray", "missing"),
    function(e1, e2) register_delayed_op(e1, .Generic, Largs=list(0L))
)
setMethod("-", c("DelayedArray", "missing"),
    function(e1, e2) register_delayed_op(e1, .Generic, Largs=list(0L))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmax2() and pmin2()
###
### We treat them like the binary operators of the "Ops" group generics.
###

setGeneric("pmax2", function(e1, e2) standardGeneric("pmax2"))
setGeneric("pmin2", function(e1, e2) standardGeneric("pmin2"))

### Mimicking how the "Ops" members combine the "dim", "names", and "dimnames"
### attributes of the 2 operands.
.check_and_combine_dims <- function(e1, e2)
{
    dim1 <- dim(e1)
    dim2 <- dim(e2)
    if (is.null(dim1))
        return(dim2)
    if (is.null(dim2))
        return(dim1)
    if (!identical(dim1, dim2))
        stop("non-conformable arrays")
    dim1
}

.combine_names <- function(e1, e2)
{
    len1 <- length(e1)
    len2 <- length(e2)
    names1 <- names(e1)
    if (len1 > len2)
        return(names1)
    names2 <- names(e2)
    if (len2 > len1 || is.null(names1))
        return(names2)
    names1
}

setMethod("pmax2", c("ANY", "ANY"),
    function(e1, e2)
    {
        ans_dim <- .check_and_combine_dims(e1, e2)
        ans <- pmax(e1, e2)
        if (is.null(ans_dim)) {
            names(ans) <- .combine_names(e1, e2)
        } else {
            dim(ans) <- ans_dim
            dimnames(ans) <- .combine_dimnames(e1, e2)
        }
        ans
    }
)

setMethod("pmin2", c("ANY", "ANY"),
    function(e1, e2)
    {
        ans_dim <- .check_and_combine_dims(e1, e2)
        ans <- pmin(e1, e2)
        if (is.null(ans_dim)) {
            names(ans) <- .combine_names(e1, e2)
        } else {
            dim(ans) <- ans_dim
            dimnames(ans) <- .combine_dimnames(e1, e2)
        }
        ans
    }
)

for (.Generic in c("pmax2", "pmin2")) {
    setMethod(.Generic, c("DelayedArray", "vector"),
        function(e1, e2)
            .DelayedArray_Ops_with_right_vector(.Generic, e1, e2)
    )
    setMethod(.Generic, c("vector", "DelayedArray"),
        function(e1, e2)
            .DelayedArray_Ops_with_left_vector(.Generic, e1, e2)
    )
    setMethod(.Generic, c("DelayedArray", "DelayedArray"),
        function(e1, e2)
            .DelayedArray_Ops(.Generic, e1, e2)
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A low-level utility for putting DelayedArray object in a "straight" form
###
### Untranspose the DelayedArray object and put its rows and columns in their
### "native" order. The result is a DelayedArray object where the array
### elements are in the same order as in the seeds. This makes block-processing
### faster if the seeds are on-disk objects where the 1st dimension is the fast
### changing dimension (e.g. 5x faster if the seeds are HDF5Dataset objects).
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
    x_index <- x@index
    for (N in x@subindex) {
        if (isStrictlySorted(x_index[[N]]))
            next
        x_index[[N]] <- .straighten_index(x_index[[N]])
    }
    if (!identical(x@index, x_index))
        x@index <- x_index
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

.DelayedArray_block_anyNA <- function(x, recursive=FALSE)
{
    REDUCE <- anyNA
    COMBINE <- `||`
    init <- FALSE
    BREAKIF <- identity

    x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
    block_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init, BREAKIF)
}

setMethod("anyNA", "DelayedArray", .DelayedArray_block_anyNA)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Summary" group generic
###
### Members: max, min, range, sum, prod, any, all
###

.collect_objects <- function(x, ...)
{
    if (missing(x)) {
        objects <- unname(list(...))
    } else {
        objects <- unname(list(x, ...))
    }
    NULL_idx <- which(S4Vectors:::sapply_isNULL(objects))
    if (length(NULL_idx) != 0L)
        objects <- objects[-NULL_idx]
    is_array_like <- function(x) is(x, "DelayedArray") || is.array(x)
    if (!all(vapply(objects, is_array_like, logical(1))))
        stop("the objects to combine must be array-like objects (or NULLs)")
    objects
}

.DelayedArray_block_Summary <- function(.Generic, x, ..., na.rm=FALSE)
{
    objects <- .collect_objects(x, ...)

    GENERIC <- match.fun(.Generic)
    REDUCE <- function(subarray) {
        ## We get a warning if 'subarray' is empty (which can't happen, blocks
        ## can't be empty) or if 'na.rm' is TRUE and 'subarray' contains only
        ## NA's or NaN's.
        reduced <- tryCatch(GENERIC(subarray, na.rm=na.rm), warning=identity)
        if (is(reduced, "warning"))
            return(NULL)
        reduced
    }
    COMBINE <- function(init, reduced) {
        if (is.null(init) && is.null(reduced))
            return(NULL)
        GENERIC(init, reduced)
    }
    init <- NULL
    BREAKIF <- function(init) {
        if (is.null(init))
            return(FALSE)
        switch(.Generic,
            max=         is.na(init) || init == Inf,
            min=         is.na(init) || init == -Inf,
            range=       is.na(init[[1L]]) || all(init == c(-Inf, Inf)),
            sum=, prod=  is.na(init),
            any=         identical(init, TRUE),
            all=         identical(init, FALSE),
            FALSE)  # fallback (actually not needed)
    }

    for (x in objects) {
        if (.Generic %in% c("sum", "prod")) {
            x <- .straighten(x, untranspose=TRUE)
        } else {
            x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
        }
        init <- block_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init, BREAKIF)
    }
    if (is.null(init))
        init <- GENERIC()
    init
}

setMethod("Summary", "DelayedArray",
    function(x, ..., na.rm=FALSE)
        .DelayedArray_block_Summary(.Generic, x, ..., na.rm=na.rm)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### Same arguments as base::mean.default().
.DelayedArray_block_mean <- function(x, trim=0, na.rm=FALSE)
{
    if (!identical(trim, 0))
        stop("\"mean\" method for DelayedArray objects ",
             "does not support the 'trim' argument yet")

    REDUCE <- function(subarray) {
        tmp <- as.vector(subarray, mode="numeric")
        subarray_sum <- sum(tmp, na.rm=na.rm)
        subarray_nval <- length(tmp)
        if (na.rm)
            subarray_nval <- subarray_nval - sum(is.na(tmp))
        c(subarray_sum, subarray_nval)
    }
    COMBINE <- `+`
    init <- numeric(2)  # sum and nval
    BREAKIF <- function(init) is.na(init[[1L]])

    x <- .straighten(x, untranspose=TRUE)
    ans <- block_REDUCE_and_COMBINE(x, REDUCE, COMBINE, init, BREAKIF)
    ans[[1L]] / ans[[2L]]
}

### S3/S4 combo for mean.DelayedArray
mean.DelayedArray <- function(x, trim=0, na.rm=FALSE, ...)
    .DelayedArray_block_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "DelayedArray", .DelayedArray_block_mean)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### apply()
###

setGeneric("apply", signature="X")

.simplify_apply_answer <- function(ans)
{
    if (!all(vapply(ans, is.atomic, logical(1), USE.NAMES=FALSE)))
        return(ans)  # won't simplify

    ans_lens <- lengths(ans, use.names=FALSE)
    mat_nrow <- ans_lens[[1L]]
    if (!all(ans_lens == mat_nrow))
        return(ans)  # won't simplify

    mat_data <- unlist(unname(ans))
    if (mat_nrow == 0L)
        return(mat_data)  # zero-length atomic vector

    mat_colnames <- names(ans)
    if (mat_nrow == 1L)
        return(setNames(mat_data, mat_colnames))  # atomic vector parallel
                                                  # to 'ans'

    ## Simplify as matrix.
    mat_data_names <- names(mat_data)  # comes from the 'ans' inner names
    if (is.null(mat_data_names)) {
        mat_rownames <- NULL
    } else {
        mat_rownames <- head(mat_data_names, n=mat_nrow)
        if (!all(mat_data_names == mat_rownames))
            mat_rownames <- NULL
    }
    if (is.null(mat_rownames) && is.null(mat_colnames)) {
        mat_dimnames <- NULL
    } else {
        mat_dimnames <- list(mat_rownames, mat_colnames)
    }
    matrix(mat_data, ncol=length(ans), dimnames=mat_dimnames)
}

### MARGIN must be a single integer.
.DelayedArray_apply <- function(X, MARGIN, FUN, ...)
{
    FUN <- match.fun(FUN)
    X_dim <- dim(X)
    if (!isSingleNumber(MARGIN))
        stop("'MARGIN' must be a single integer")
    if (!is.integer(MARGIN))
        MARGIN <- as.integer(MARGIN)
    if (MARGIN < 1L || MARGIN > length(X_dim))
        stop("'MARGIN' must be >= 1 and <= length(dim(X))")

    if (X_dim[[MARGIN]] == 0L) {
        ## base::apply seems to be doing something like that!
        ans <- FUN(X, ...)
        return(as.vector(ans[0L]))
    }

    ## TODO: Try using sapply() instead of lapply(). Maybe we're lucky
    ## and it achieves the kind of simplification that we're doing with
    ## .simplify_apply_answer() so we can get rid of .simplify_apply_answer().
    ans_names <-  dimnames(X)[[MARGIN]]
    ans <- lapply(setNames(seq_len(X_dim[[MARGIN]]), ans_names),
        function(i) {
            subscript <- rep.int(alist(foo=), length(X_dim))
            subscript[[MARGIN]] <- i
            args <- c(list(X), subscript)
            slice <- do.call(`[`, args)
            dim(slice) <- dim(slice)[-MARGIN]
            FUN(slice, ...)
        })

    ## Try to simplify the answer.
    .simplify_apply_answer(ans)
}

setMethod("apply", "DelayedArray", .DelayedArray_apply)

