### =========================================================================
### HDF5Matrix objects
### -------------------------------------------------------------------------


### Extending DataTable gives us a few things for free (head(), tail(),
### etc...)
setClass("HDF5Matrix",
    contains=c("HDF5Array", "DataTable"),
    representation(
        ## x@N1 and x@N2 must be 2 integers such that
        ##     1 <= x@N1 < x@N2 <= length(x@index)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

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

### array -> HDF5Array

.from_matrix_to_HDF5Matrix <- function(from)
{
    as(as(from, "HDF5Array"), "HDF5Matrix")
}

setAs("matrix", "HDF5Matrix", .from_matrix_to_HDF5Matrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

HDF5Matrix <- function(file, group, name)
{
    hdf5array <- HDF5Array(file, group, name)
    as(hdf5array, "HDF5Matrix")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.extract_HDF5Matrix_subset <- function(x, i, j, ..., drop=TRUE)
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

    ans1 <- as.matrix(x[idx1, , drop=FALSE])
    ans2 <- as.matrix(x[idx2, , drop=FALSE])
    storage.mode(ans1) <- storage.mode(ans2) <- "character"
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

    ans1 <- as.matrix(x[ , idx1, drop=FALSE])
    ans2 <- as.matrix(x[ , idx2, drop=FALSE])
    storage.mode(ans1) <- storage.mode(ans2) <- "character"
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
            ans <- as.matrix(x)
            storage.mode(ans) <- "character"
        } else {
            ans <- .csplit_HDF5Matrix(x, 4L, 4L)
        }
    } else {
        if (x_ncol <= 9L) {
            ans <- .rsplit_HDF5Matrix(x, nhead, ntail)
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

