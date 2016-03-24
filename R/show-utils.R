### =========================================================================
### Compact display of an array-like object
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 1D array
###

.split_1D_array_names <- function(x_names, idx1, idx2)
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

.prepare_1D_array_sample <- function(x, n1, n2)
{
    x_len <- length(x)
    x_names <- names(x)
    if (x_len <= n1 + n2 + 1L) {
        ans <- format(as.vector(x))
        idx1 <- seq_len(x_len)
        idx2 <- integer(0)
        names(ans) <- .split_1D_array_names(x_names, idx1, idx2)[idx1]
    } else {
        idx1 <- seq_len(n1)
        idx2 <- seq(to=x_len, by=1L, length.out=n2)
        ans1 <- format(as.vector(x[idx1]))
        ans2 <- format(as.vector(x[idx2]))
        ans <- c(ans1, ".", ans2)
        names(ans) <- .split_1D_array_names(x_names, idx1, idx2)
    }
    ans
}

.print_1D_array_data <- function(x, n1, n2)
{
    stopifnot(length(dim(x)) == 1L)
    out <- .prepare_1D_array_sample(x, n1, n2)
    quote <- type(x) == "character"
    print(out, quote=quote, right=TRUE, max=length(out))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2D array
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

.rsplit_2D_array_data <- function(x, m1, m2)
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

.csplit_2D_array_data <- function(x, n1, n2)
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

.split_2D_array_data <- function(x, m1, m2, n1, n2)
{
    x_ncol <- ncol(x)
    x_colnames <- colnames(x)
    idx1 <- seq_len(n1)
    idx2 <- seq(to=x_ncol, by=1L, length.out=n2)

    x1 <- x[ , idx1, drop=FALSE]
    x2 <- x[ , idx2, drop=FALSE]
    ans1 <- .rsplit_2D_array_data(x1, m1, m2)
    ans2 <- .rsplit_2D_array_data(x2, m1, m2)
    dots <- rep.int(".", nrow(ans1))
    ans <- cbind(ans1, matrix(dots, ncol=1L), ans2)

    colnames(ans) <- .split_colnames(x_colnames, idx1, idx2)
    ans
}

.prepare_2D_array_sample <- function(x, m1, m2, n1, n2)
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
            ans <- .csplit_2D_array_data(x, n1, n2)
        }
    } else {
        if (x_ncol <= n1 + n2 + 1L) {
            ans <- .rsplit_2D_array_data(x, m1, m2)
            ## Only needed because of this bug in base::print.default:
            ##   https://stat.ethz.ch/pipermail/r-devel/2016-March/072479.html
            ## TODO: Remove when the bug is fixed.
            if (is.null(colnames(ans))) {
                idx1 <- seq_len(ncol(ans))
                idx2 <- integer(0)
                colnames(ans) <- .split_colnames(NULL, idx1, idx2)[idx1]
            }
        } else {
            ans <- .split_2D_array_data(x, m1, m2, n1, n2)
        }
    }
    ans
}

.print_2D_array_data <- function(x, m1, m2, n1, n2)
{
    stopifnot(length(dim(x)) == 2L)
    out <- .prepare_2D_array_sample(x, m1, m2, n1, n2)
    quote <- type(x) == "character"
    print(out, quote=quote, right=TRUE, max=length(out))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Array of arbitrary dimensions
###

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
        slice <- extract_array_block1(x, subscript)
        dim(slice) <- dim(slice)[1:2]
        .print_2D_array_data(slice, m1, m2, n1, n2)
        cat("\n")
    }
}

.print_nD_array_data <- function(x, n1, n2)
{
    x_dim <- dim(x)
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
    blocks <- new("ArrayBlocks", dim=x_dim, N=3L, by=1L)
    nblock <- length(blocks)
    if (nblock <= z1 + z2 + 1L) {
        idx <- seq_len(nblock)
        .print_2D_slices(x, m1, m2, n1, n2, blocks, idx)
    } else {
        idx1 <- seq_len(z1)
        idx2 <- seq(to=nblock, by=1L, length.out=z2)
        .print_2D_slices(x, m1, m2, n1, n2, blocks, idx1)
        cat("...\n\n")
        .print_2D_slices(x, m1, m2, n1, n2, blocks, idx2)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show_compact_array()
###

.print_array_data <- function(x, n1, n2)
{
    x_dim <- dim(x)
    if (length(x_dim) == 1L)
        return(.print_1D_array_data(x, n1, n2))
    if (length(x_dim) == 2L) {
        nhead <- get_showHeadLines()
        ntail <- get_showTailLines()
        return(.print_2D_array_data(x, nhead, ntail, n1, n2))
    }
    .print_nD_array_data(x, n1, n2)
}

show_compact_array <- function(object)
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
        if (object_type == "integer") {
            n1 <- n2 <- 4L
        } else {
            n1 <- 3L
            n2 <- 2L
        }
        .print_array_data(object, n1, n2)
    }
}

