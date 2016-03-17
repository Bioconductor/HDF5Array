### =========================================================================
### HDF5Matrix objects
### -------------------------------------------------------------------------


### Extending DataTable gives us a few things for free (head(), tail(),
### etc...)
setClass("HDF5Matrix",
    contains=c("HDF5Array", "DataTable"),
    representation(
        ## x@N1 and x@N2 must be 2 integers such that
        ##     1 <= x@N1 < x@N2 <= length(x@h5index)
        N1="integer",  # single integer
        N2="integer"   # single integer
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.is_valid_N1 <- function(N1, ndim)
{
    isSingleInteger(N1) && N1 >= 1L && N1 <= ndim
}

.validate_HDF5Matrix <- function(x)
{
    if (!.is_valid_N1(x@N1, length(x@h5index))
     || !.is_valid_N1(x@N2, length(x@h5index)))
        return(wmsg("'x@N1' and 'x@N2' must be single integers ",
                    ">= 1 and <= 'length(x@h5index)'"))
    if (x@N1 >= x@N2)
        return("'x@N1' must be < 'x@N2'")
    array_dim <- lengths(x@h5index)
    if (!all(array_dim[-c(x@N1, x@N2)] == 1L))
        return(wmsg("'x@N1' and 'x@N2' are incompatible with the ",
                    "dimensions of the underlying HDF5Array object"))
    TRUE
}

setValidity2("HDF5Matrix", .validate_HDF5Matrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###
### Defining the internal index() getter and setter is enough to make all the
### HDF5Array accessors (length, isEmpty, dim, dimnames, dimnames<-) work on
### an HDF5Matrix object.

setMethod("index", "HDF5Matrix",
    function(x) x@h5index[c(x@N1, x@N2)]
)
setReplaceMethod("index", "HDF5Matrix", 
    function(x, value) { x@h5index[c(x@N1, x@N2)] <- value; x }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.from_HDF5Array_to_HDF5Matrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) < 2L)
        stop(wmsg(class(from), " object with less than 2 dimensions ",
                  "cannot be coerced to an HDF5Matrix object at the moment"))
    if (length(from_dim) == 2L) {
        N1 <- 1L
        N2 <- 2L
    } else {
        idx <- which(from_dim != 1L)
        if (length(idx) > 2L)
            stop(wmsg(class(from), " object with more than 2 effective ",
                      "dimensions cannot be coerced to an HDF5Matrix object. ",
                      slicing_tip))
        if (length(idx) == 2L) {
            N1 <- idx[[1L]]
            N2 <- idx[[2L]]
        } else if (length(idx) == 0L) {
            N1 <- 1L
            N2 <- 2L
        } else {
            ## length(idx) == 1L
            N1 <- idx[[1L]]
            if (N1 == length(from_dim))
                stop(wmsg("A ", class(from), " object where the only ",
                          "effective dimension is its last dimension cannot ",
                          "be coerced to a HDF5Matrix object at the moment"))
            N2 <- N1 + 1L
        }
    }
    new2("HDF5Matrix", from, N1=N1, N2=N2)
}

setAs("HDF5Array", "HDF5Matrix", .from_HDF5Array_to_HDF5Matrix)

### array -> HDF5Array

.from_matrix_to_HDF5Matrix <- function(from)
{
    as(as(from, "HDF5Array"), "HDF5Matrix")
}

setAs("matrix", "HDF5Matrix", .from_matrix_to_HDF5Matrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

HDF5Matrix <- function(file, group, name, type=NA)
{
    hdf5array <- HDF5Array(file, group, name, type=type)
    as(hdf5array, "HDF5Matrix")
}

### 'x' must be an array-like object with 3 dimensions.
### 'MARGIN' is the dimension to drop.
make_HDF5Matrix_from_3D_array <- function(x, MARGIN)
{
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (x_ndim != 3L)
        stop("'x' must have 3 dimensions")
    if (!isSingleNumber(MARGIN))
        stop("'MARGIN' must be a single integer")
    if (!is.integer(MARGIN))
        MARGIN <- as.integer(MARGIN)
    if (MARGIN < 1L || MARGIN > x_ndim)
        stop("'MARGIN' must be >= 1 and <= length(dim(x))")
    if (x_dim[[MARGIN]] != 1L)
        stop("'dim(x)[[MARGIN]]' must be 1")
    if (!is(x, "HDF5Array"))
        x <- as(x, "HDF5Array")

    if (x@is_transposed)
        MARGIN <- x_ndim + 1L - MARGIN
    tmp <- seq_along(x_dim)[-MARGIN]
    N1 <- tmp[[1L]]
    N2 <- tmp[[2L]]
    new2("HDF5Matrix", x, N1=N1, N2=N2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

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

.rsplit_HDF5Matrix <- function(x, m1, m2)
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

.csplit_HDF5Matrix <- function(x, n1, n2)
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

.split_HDF5Matrix <- function(x, m1, m2, n1, n2)
{
    x_ncol <- ncol(x)
    x_colnames <- colnames(x)
    idx1 <- seq_len(n1)
    idx2 <- seq(to=x_ncol, by=1L, length.out=n2)

    x1 <- x[ , idx1, drop=FALSE]
    x2 <- x[ , idx2, drop=FALSE]
    ans1 <- .rsplit_HDF5Matrix(x1, m1, m2)
    ans2 <- .rsplit_HDF5Matrix(x2, m1, m2)
    dots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(dots, ncol=1L), ans2)

    colnames(ans) <- .split_colnames(x_colnames, idx1, idx2)
    ans
}

.prepare_HDF5Matrix_sample <- function(x)
{
    nhead <- get_showHeadLines()
    ntail <- get_showTailLines()
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_nrow <= nhead + ntail + 1L) {
        if (x_ncol <= 9L) {
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
            ans <- .csplit_HDF5Matrix(x, 4L, 4L)
        }
    } else {
        if (x_ncol <= 9L) {
            ans <- .rsplit_HDF5Matrix(x, nhead, ntail)
            ## Only needed because of this bug in base::print.default:
            ##   https://stat.ethz.ch/pipermail/r-devel/2016-March/072479.html
            ## TODO: Remove when the bug is fixed.
            if (is.null(colnames(ans))) {
                idx1 <- seq_len(ncol(ans))
                idx2 <- integer(0)
                colnames(ans) <- .split_colnames(NULL, idx1, idx2)[idx1]
            }
        } else {
            ans <- .split_HDF5Matrix(x, nhead, ntail, 4L, 4L)
        }
    }
    ans
}

setMethod("show", "HDF5Matrix",
    function(object) 
    {
        show_HDF5Array_topline(object)
        #if (isEmpty(object)) {
        #    cat("\n")
        #} else {
            cat(":\n")
            out <- .prepare_HDF5Matrix_sample(object)
            print(out, quote=FALSE, right=TRUE, max=length(out))
        #}
    }
)

