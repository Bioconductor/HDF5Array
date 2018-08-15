### =========================================================================
### TENxMatrix objects
### -------------------------------------------------------------------------


setClass("TENxMatrixSeed",
    contains="Array",
    representation(
        filepath="character",    # Absolute path to the HDF5 file so the
                                 # object doesn't break when the user
                                 # changes the working directory (e.g. with
                                 # setwd()).
        group="character",       # Name of the group in the HDF5 file
                                 # containing the 10x Genomics data.
        dim="integer",
        dimnames="list",
        col_ranges="data.frame"  # Can't use an IRanges object for this at the
                                 # moment because they don't support Linteger
                                 # start/end values yet.
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### Does NOT access the file.
setMethod("path", "TENxMatrixSeed", function(object) object@filepath)

### Just a placeholder for now. Doesn't actually allow changing the path of
### the object yet.
setReplaceMethod("path", "TENxMatrixSeed",
    function(object, value)
    {
        new_filepath <- normarg_path(value, "the supplied path",
                                            "10x Genomics dataset")
        old_filepath <- path(object)
        if (new_filepath != old_filepath)
            stop(wmsg("changing the path of a TENxMatrixSeed object ",
                      "is not supported yet"))
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() getter
###

### Does NOT access the file.
setMethod("dim", "TENxMatrixSeed", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dimnames() getter
###

### Does NOT access the file.
setMethod("dimnames", "TENxMatrixSeed", function(x) x@dimnames)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level internal disk data extractors
###

### All the 10xGenomics components are monodimensional.
.get_TENx_component <- function(filepath, group, name, idx=NULL)
{
    name <- paste0(group, "/", name)
    if (!is.null(idx))
        idx <- list(idx)
    as.vector(h5read(filepath, name, index=idx))
}

### Return the dimensions of the matrix.
.get_shape <- function(filepath, group)
    .get_TENx_component(filepath, group, "shape")

.get_indptr <- function(filepath, group)
{
    name <- paste0(group, "/indptr")
    as.vector(h5read(filepath, name, bit64conversion="double"))
}

.get_data <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "data", idx=idx)

.linear_get_data <- function(filepath, group, start=NULL, count=NULL)
{
    name <- paste0(group, "/data")
    as.vector(h5read(filepath, name, start=start, count=count))
}

### Return 0-based row indices.
.get_row_indices <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "indices", idx=idx)

.linear_get_row_indices <- function(filepath, group, start=NULL, count=NULL)
{
    name <- paste0(group, "/indices")
    as.vector(h5read(filepath, name, start=start, count=count))
}

### Return the rownames of the matrix.
.get_genes <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "genes", idx=idx)

### Currently unused.
.get_gene_names <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "gene_names", idx=idx)

### Return the colnames of the matrix.
.get_barcodes <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "barcodes", idx=idx)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_data_indices_by_col()
###

### S4Vectors:::fancy_mseq() does not accept 'offset' of type double yet so
### we implement a version that does.
### Will this work if sum(lengths) is > .Machine$integer.max?
.fancy_mseq <- function(lengths, offset=0)
{
    lengths_len <- length(lengths)
    if (lengths_len == 0L)
        return(numeric(0))
    offsets <- offset - cumsum(c(0L, lengths[-lengths_len]))
    seq_len(sum(lengths)) + rep.int(offsets, lengths)
}

### 'j' must be an integer vector containing valid col indices.
### Return data indices in a NumericList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
.get_data_indices_by_col <- function(x, j)
{
    col_ranges <- S4Vectors:::extract_data_frame_rows(x@col_ranges, j)
    start2 <- col_ranges[ , "start"]
    width2 <- col_ranges[ , "width"]
    idx2 <- .fancy_mseq(width2, offset=start2 - 1L)
    ### Will this work if 'idx2' is a long vector?
    relist(idx2, PartitioningByWidth(width2))
}

### 'j' must be NULL or an integer vector containing valid col indices.
### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
.get_data_by_col <- function(x, j)
{
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    data <- .get_data(x@filepath, x@group, idx=idx2)
    relist(data, data_indices)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_sparse_data_from_TENxMatrixSeed()
###

### 'i' must be NULL or an integer vector containing valid row indices.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
### Both 'i' and 'j' can contain duplicates. Duplicates in 'i' have no effect
### on the output but duplicates in 'j' will produce duplicates in the output.
.random_extract_sparse_data_from_TENxMatrixSeed <- function(x, i, j)
{
    stopifnot(is.null(i) || is.numeric(i), is.numeric(j))
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    row_indices <- .get_row_indices(x@filepath, x@group, idx=idx2) + 1L
    col_indices <- rep.int(j, lengths(data_indices))
    if (!is.null(i)) {
        keep_idx <- which(row_indices %in% i)
        idx2 <- idx2[keep_idx]
        row_indices <- row_indices[keep_idx]
        col_indices <- col_indices[keep_idx]
    }
    data <- .get_data(x@filepath, x@group, idx=idx2)
    list(i=row_indices, j=col_indices, data=data)
}

### 'j1' and 'j2' must be 2 single integers representing a valid range of col
### indices.
.extract_sparse_data_from_TENxMatrixSeed_cols <- function(x, j1, j2)
{
    j12 <- j1:j2
    start <- x@col_ranges[j1, "start"]
    count_per_col <- x@col_ranges[j12, "width"]
    count <- sum(count_per_col)
    data <- .linear_get_data(x@filepath, x@group, start=start, count=count)
    row_indices <- .linear_get_row_indices(x@filepath, x@group,
                                           start=start, count=count) + 1L
    col_indices <- rep.int(j12, count_per_col)
    list(i=row_indices, j=col_indices, data=data)
}

### 'i' (or 'j') must be NULL or an integer vector containing valid row (or
### col) indices. 'j' should not be empty.
### The output is not affected by duplicates in 'i' or 'j'.
.linear_extract_sparse_data_from_TENxMatrixSeed <- function(x, i, j)
{
    stopifnot(is.null(i) || is.numeric(i))
    if (is.null(j)) {
        j1 <- 1L
        j2 <- ncol(x)
    } else {
        stopifnot(is.numeric(j), length(j) != 0L)
        j1 <- min(j)
        j2 <- max(j)
    }
    sparse_data <- .extract_sparse_data_from_TENxMatrixSeed_cols(x, j1, j2)
    if (is.null(i) && is.null(j))
        return(sparse_data)
    row_indices <- sparse_data$i
    col_indices <- sparse_data$j
    data <- sparse_data$data
    if (is.null(i)) {
        keep_me <- col_indices %in% j
    } else if (is.null(j)) {
        keep_me <- row_indices %in% i
    } else {
        keep_me <- (row_indices %in% i) & (col_indices %in% j)
    }
    keep_idx <- which(keep_me)
    data <- data[keep_idx]
    row_indices <- row_indices[keep_idx]
    col_indices <- col_indices[keep_idx]
    list(i=row_indices, j=col_indices, data=data)
}

.normarg_method <- function(method, j)
{
    if (method != "auto")
        return(method)
    if (is.null(j))
        return("linear")
    if (length(j) == 0L)
        return("random")
    j1 <- min(j)
    j2 <- max(j)
    ## 'ratio' is > 0 and <= 1. A value close to 1 indicates that the columns
    ## to extract are close from each other (a value of 1 indicating that
    ## they are adjacent e.g. j <- 18:25). A value close to 0 indicates that
    ## they are far apart from each other i.e. that they are separated by many
    ## columns that are not requested. The "linear" method is very efficient
    ## when 'ratio' is close to 1. It is so much more efficient than the
    ## "random" method (typically 10x or 20x faster) that we choose it when
    ## 'ratio' is >= 0.2
    ratio <- length(j) / (j2 - j1 + 1L)
    if (ratio >= 0.2) "linear" else "random"
}

### 'i' (or 'j') must be NULL or an integer vector containing valid row (or
### col) indices.
### Duplicates in 'i' are ok and won't affect the output.
### Duplicates in 'j' are ok but might introduce duplicates in the output
### so should be avoided.
### Return the sparse data in a list (NOT a data.frame) with 3 parallel
### components: i, j, data. We choose list over data.frame to avoid the
### overhead and some annoyances that come with the latter.
.extract_sparse_data_from_TENxMatrixSeed <-
    function(x, i, j, method=c("auto", "random", "linear"))
{
    method <- match.arg(method)
    method <- .normarg_method(method, j)
    if (method == "random") {
        .random_extract_sparse_data_from_TENxMatrixSeed(x, i, j)
    } else {
        .linear_extract_sparse_data_from_TENxMatrixSeed(x, i, j)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.make_matrix_from_sparse_data <- function(dim, sparse_data=NULL)
{
    ans <- array(0L, dim=dim)
    ans[cbind(sparse_data$i, sparse_data$j)] <- sparse_data$data
    ans
}

.extract_array_from_TENxMatrixSeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L))
        return(array(integer(0), dim=ans_dim))  # return an empty matrix
    i <- index[[1L]]
    ui <- if (is.null(i)) NULL else unique(i)
    j <- index[[2L]]
    uj <- if (is.null(j)) NULL else unique(j)
    sparse_data <- .extract_sparse_data_from_TENxMatrixSeed(x, ui, uj)
    if (is.null(ui)) {
        umat_nrow <- nrow(x)
        i2ui <- NULL
    } else {
        sparse_data$i <- match(sparse_data$i, ui)
        umat_nrow <- length(ui)
        i2ui <- match(i, ui)
    }
    if (is.null(uj)) {
        umat_ncol <- ncol(x)
        j2uj <- NULL
    } else {
        sparse_data$j <- match(sparse_data$j, uj)
        umat_ncol <- length(uj)
        j2uj <- match(j, uj)
    }
    umat_dim <- c(umat_nrow, umat_ncol)
    umat <- .make_matrix_from_sparse_data(umat_dim, sparse_data)
    DelayedArray:::subset_by_Nindex(umat, list(i2ui, j2uj))
}

setMethod("extract_array", "TENxMatrixSeed",
    .extract_array_from_TENxMatrixSeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkdim() getter
###

### Does NOT access the file.
setMethod("chunkdim", "TENxMatrixSeed", function(x) c(nrow(x), 1L))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TENxMatrixSeed constructor
###

TENxMatrixSeed <- function(filepath, group="mm10")
{
    filepath <- normarg_path(filepath, "'filepath'", "10x Genomics dataset")
    if (!isSingleString(group))
        stop(wmsg("'group' must be a single string specifying the name ",
                  "of the group in the HDF5 file containing the ",
                  "10x Genomics data"))
    if (group == "")
        stop(wmsg("'group' cannot be the empty string"))

    ## dim
    dim <- .get_shape(filepath, group)
    stopifnot(length(dim) == 2L)

    ## dimnames
    rownames <- .get_genes(filepath, group)
    stopifnot(length(rownames) == dim[[1L]])
    colnames <- .get_barcodes(filepath, group)
    stopifnot(length(colnames) == dim[[2L]])
    dimnames <- list(rownames, colnames)

    ## col_ranges
    ## "/data" and "/indices" are monodimensional arrays of length >= 2^31 so
    ## we need to call h5dim() with 'as.integer=FALSE'.
    data_len <- h5dim(filepath, paste0(group, "/data"), as.integer=FALSE)
    stopifnot(length(data_len) == 1L)
    indices_len <- h5dim(filepath, paste0(group, "/indices"), as.integer=FALSE)
    stopifnot(identical(data_len, indices_len))
    indptr <- .get_indptr(filepath, group)
    stopifnot(length(indptr) == dim[[2L]] + 1L,
              indptr[[1L]] == 0L,
              indptr[[length(indptr)]] == data_len)
    col_ranges <- data.frame(start=indptr[-length(indptr)] + 1,
                             width=as.integer(diff(indptr)))

    new2("TENxMatrixSeed", filepath=filepath,
                           group=group,
                           dim=dim,
                           dimnames=dimnames,
                           col_ranges=col_ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "TENxMatrixSeed",
    function(object)
    {
        cat(DelayedArray:::array_as_one_line_summary(object), ":\n", sep="")
        cat("# dirname: ", dirname(object), "\n", sep="")
        cat("# basename: ", basename(object), "\n", sep="")
        cat("# group: ", object@group, "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TENxMatrix objects
###

setClass("TENxMatrix", contains="DelayedMatrix")

.validate_TENxMatrix <- function(x)
{
    if (!is(x@seed, "TENxMatrixSeed"))
        return(wmsg("'x@seed' must be a TENxMatrixSeed object"))
    TRUE
}

setValidity2("TENxMatrix", .validate_TENxMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "TENxMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="TENxMatrix")
)

### Works directly on a TENxMatrixSeed object, in which case it must be
### called with a single argument.
TENxMatrix <- function(filepath, group="mm10")
{
    if (is(filepath, "TENxMatrixSeed")) {
        if (!missing(group))
            stop(wmsg("TENxMatrix() must be called with a single argument ",
                      "when passed a TENxMatrixSeed object"))
        seed <- filepath
    } else {
        seed <- TENxMatrixSeed(filepath, group)
    }
    DelayedArray(seed)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractNonZeroValuesByCol()
###

### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
setGeneric("extractNonZeroValuesByCol", signature="x",
    function(x, j) standardGeneric("extractNonZeroValuesByCol")
)

setMethod("extractNonZeroValuesByCol", "TENxMatrixSeed",
    function(x, j)
    {
        j <- DelayedArray:::normalizeSingleBracketSubscript2(j, ncol(x),
                                                             colnames(x))
        .get_data_by_col(x, j)
    }
)

setMethod("extractNonZeroValuesByCol", "TENxMatrix",
    function(x, j) extractNonZeroValuesByCol(x@seed, j)
)

