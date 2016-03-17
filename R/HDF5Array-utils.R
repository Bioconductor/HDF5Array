### =========================================================================
### Common operations on HDF5Array objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various unary operators + the "Math" group generic
###

### All these operations return an HDF5Array object of the same dimensions
### as 'x'. This object points to the same HDF5 dataset than 'x' but has the
### operation stored in it (in the delayed_ops slot) as a delayed operation.

setMethod("is.na", "HDF5Array", function(x) register_delayed_op(x, "is.na"))

setMethod("!", "HDF5Array", function(x) register_delayed_op(x, "!"))

setMethod("Math", "HDF5Array", function(x) register_delayed_op(x, .Generic))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "Ops" group generics
###
### Arith members: "+", "-", "*", "/", "^", "%%", "%/%"
### Compare members: ==, !=, <=, >=, <, >
### Logic members: &, |
###

### Return an HDF5Array object of the same dimensions as 'e1'. This object
### points to the same HDF5 dataset than 'e1' but has the operation stored
### in it (in the delayed_ops slot) as a delayed operation.
.HDF5Array_delayed_Ops_with_right_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e1, "HDF5Array"))
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

### Return an HDF5Array object of the same dimensions as 'e2'. This object
### points to the same HDF5 dataset than 'e2' but has the operation stored
### in it (in the delayed_ops slot) as a delayed operation.
.HDF5Array_delayed_Ops_with_left_vector <- function(.Generic, e1, e2)
{
    stopifnot(is(e2, "HDF5Array"))
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

### Return an HDF5Array object. This object points to its own HDF5 dataset
### stored in a new file.
.HDF5Array_block_Ops <- function(.Generic, e1, e2)
{
    if (!identical(dim(e1), dim(e2)))
        stop("non-conformable arrays")
    GENERIC <- match.fun(.Generic)

    out_file <- paste0(tempfile(), ".h5")
    out_name <- sprintf("%s %s %s", "e1", "Ops", "e2")
    ans_type <- typeof(GENERIC(match.fun(type(e1))(1), match.fun(type(e2))(1)))
    h5createFile(out_file)
    h5createDataset(out_file, out_name, dim(e1), storage.mode=ans_type)
    
    block_MAPPLY(GENERIC, e1, e2,
        out_file=out_file,
        out_name=out_name
    )

    ans <- HDF5Array(out_file, "/", out_name)
    if (is(e1, "HDF5Matrix") || is(e2, "HDF5Matrix"))
        ans <- as(ans, "HDF5Matrix")
    ans
}

.HDF5Array_Ops <- function(.Generic, e1, e2)
{
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (identical(e1_dim, e2_dim))
        return(.HDF5Array_block_Ops(.Generic, e1, e2))
    ## Effective dimensions.
    effdim_idx1 <- which(e1_dim != 1L)
    effdim_idx2 <- which(e2_dim != 1L)
    if ((length(effdim_idx1) == 1L) == (length(effdim_idx2) == 1L))
        stop("non-conformable arrays")
    if (length(effdim_idx1) == 1L) {
        .HDF5Array_delayed_Ops_with_left_vector(.Generic, e1, e2)
    } else {
        .HDF5Array_delayed_Ops_with_right_vector(.Generic, e1, e2)
    }
}

setMethod("Ops", c("HDF5Array", "vector"),
    function(e1, e2)
        .HDF5Array_delayed_Ops_with_right_vector(.Generic, e1, e2)
)

setMethod("Ops", c("vector", "HDF5Array"),
    function(e1, e2)
        .HDF5Array_delayed_Ops_with_left_vector(.Generic, e1, e2)
)

setMethod("Ops", c("HDF5Array", "HDF5Array"),
    function(e1, e2)
        .HDF5Array_Ops(.Generic, e1, e2)
)

### Support unary operators "+" and "-".
setMethod("+", c("HDF5Array", "missing"),
    function(e1, e2) register_delayed_op(e1, .Generic, Largs=list(0L))
)
setMethod("-", c("HDF5Array", "missing"),
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
    setMethod(.Generic, c("HDF5Array", "vector"),
        function(e1, e2)
            .HDF5Array_delayed_Ops_with_right_vector(.Generic, e1, e2)
    )
    setMethod(.Generic, c("vector", "HDF5Array"),
        function(e1, e2)
            .HDF5Array_delayed_Ops_with_left_vector(.Generic, e1, e2)
    )
    setMethod(.Generic, c("HDF5Array", "HDF5Array"),
        function(e1, e2)
            .HDF5Array_Ops(.Generic, e1, e2)
    )
}


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
    APPLY <- anyNA
    REDUCE <- `||`
    reduced <- FALSE
    BREAKIF <- identity

    x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
    block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, BREAKIF)
}

setMethod("anyNA", "HDF5Array", .HDF5Array_block_anyNA)


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
    is_HDF5Array_or_array <- function(x) is(x, "HDF5Array") || is.array(x)
    if (!all(vapply(objects, is_HDF5Array_or_array, logical(1))))
        stop("the objects to combine must be HDF5Array objects (or NULLs)")
    objects
}

.HDF5Array_block_Summary <- function(.Generic, x, ..., na.rm=FALSE)
{
    objects <- .collect_objects(x, ...)

    GENERIC <- match.fun(.Generic)
    APPLY <- function(subarray) {
        ## We get a warning if 'subarray' is empty (which can't happen, blocks
        ## can't be empty) or if 'na.rm' is TRUE and 'subarray' contains only
        ## NA's or NaN's.
        val <- tryCatch(GENERIC(subarray, na.rm=na.rm), warning=identity)
        if (is(val, "warning"))
            return(NULL)
        val
    }
    REDUCE <- function(reduced, val) {
        if (is.null(reduced) && is.null(val))
            return(NULL)
        GENERIC(reduced, val)
    }
    reduced <- NULL
    BREAKIF <- function(reduced) {
        if (is.null(reduced))
            return(FALSE)
        switch(.Generic,
            max=         is.na(reduced) || reduced == Inf,
            min=         is.na(reduced) || reduced == -Inf,
            range=       is.na(reduced[[1L]]) || all(reduced == c(-Inf, Inf)),
            sum=, prod=  is.na(reduced),
            any=         identical(reduced, TRUE),
            all=         identical(reduced, FALSE),
            FALSE)  # fallback (actually not needed)
    }

    for (x in objects) {
        if (.Generic %in% c("sum", "prod")) {
            x <- .straighten(x, untranspose=TRUE)
        } else {
            x <- .straighten(x, untranspose=TRUE, straighten.index=TRUE)
        }
        reduced <- block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, BREAKIF)
    }
    if (is.null(reduced))
        reduced <- GENERIC()
    reduced
}

setMethod("Summary", "HDF5Array",
    function(x, ..., na.rm=FALSE)
        .HDF5Array_block_Summary(.Generic, x, ..., na.rm=na.rm)
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

    APPLY <- function(subarray) {
        tmp <- as.vector(subarray, mode="numeric")
        subarray_sum <- sum(tmp, na.rm=na.rm)
        subarray_nval <- length(tmp)
        if (na.rm)
            subarray_nval <- subarray_nval - sum(is.na(tmp))
        c(subarray_sum, subarray_nval)
    }
    REDUCE <- `+`
    reduced <- numeric(2)  # sum and nval
    BREAKIF <- function(reduced) is.na(reduced[[1L]])

    x <- .straighten(x, untranspose=TRUE)
    reduced <- block_APPLY_REDUCE(x, APPLY, REDUCE, reduced, BREAKIF)
    reduced[[1L]] / reduced[[2L]]
}

### S3/S4 combo for mean.HDF5Array
mean.HDF5Array <- function(x, trim=0, na.rm=FALSE, ...)
    .HDF5Array_block_mean(x, trim=trim, na.rm=na.rm, ...)
setMethod("mean", "HDF5Array", .HDF5Array_block_mean)


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
.HDF5Array_apply <- function(X, MARGIN, FUN, ...)
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

    ans_names <-  dimnames(X)[[MARGIN]]
    ans <- lapply(setNames(seq_len(X_dim[[MARGIN]]), ans_names),
        function(i) {
            subscript <- rep.int(alist(foo=), length(X_dim))
            subscript[[MARGIN]] <- i
            args <- c(list(X), subscript)
            slice <- do.call("[", args)
            if (length(X_dim) == 3L && is(X, "HDF5Array"))
                slice <- make_HDF5Matrix_from_3D_array(slice, MARGIN)
            FUN(slice, ...)
        })

    ## Try to simplify the answer.
    .simplify_apply_answer(ans)
}

setMethod("apply", "HDF5Array", .HDF5Array_apply)

