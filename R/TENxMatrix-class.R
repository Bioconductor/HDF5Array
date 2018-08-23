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
setMethod("dimnames", "TENxMatrixSeed",
    function(x) DelayedArray:::simplify_NULL_dimnames(x@dimnames)
)


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
{
    if (!h5exists(filepath, paste0(group, "/genes")))
        return(NULL)
    .get_TENx_component(filepath, group, "genes", idx=idx)
}

### Currently unused.
.get_gene_names <- function(filepath, group, idx=NULL)
{
    if (!h5exists(filepath, paste0(group, "/gene_names")))
        return(NULL)
    .get_TENx_component(filepath, group, "gene_names", idx=idx)
}

### Return the colnames of the matrix.
.get_barcodes <- function(filepath, group, idx=NULL)
{
    if (!h5exists(filepath, paste0(group, "/barcodes")))
        return(NULL)
    .get_TENx_component(filepath, group, "barcodes", idx=idx)
}


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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_data_from_adjacent_cols()
###

### 'j1' and 'j2' must be 2 single integers representing a valid range of
### col indices.
### If 'as.sparse=FALSE', return a NumericList or IntegerList object parallel
### to 'j1:j2' i.e. with one list element per col index in 'j1:j2'.
### If 'as.sparse=TRUE', return a SparseArray object.
.extract_data_from_adjacent_cols <- function(x, j1, j2, as.sparse=FALSE)
{
    j12 <- j1:j2
    start <- x@col_ranges[j1, "start"]
    count_per_col <- x@col_ranges[j12, "width"]
    count <- sum(count_per_col)
    ans_nzdata <- .linear_get_data(x@filepath, x@group,
                                   start=start, count=count)
    if (!as.sparse)
        return(relist(ans_nzdata, PartitioningByWidth(count_per_col)))
    row_indices <- .linear_get_row_indices(x@filepath, x@group,
                                           start=start, count=count) + 1L
    col_indices <- rep.int(j12, count_per_col)
    ans_aind <- cbind(row_indices, col_indices, deparse.level=0L)
    SparseArray(dim(x), ans_aind, ans_nzdata, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_nonzero_data_by_col()
###

### Extract nonzero data using the "random" method. This method is
### based on h5read( , index=idx) which retrieves an arbitrary/random
### subset of the data.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
.random_extract_nonzero_data_by_col <- function(x, j)
{
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    data <- .get_data(x@filepath, x@group, idx=idx2)
    relist(data, data_indices)
}

### Extract nonzero data using the "linear" method. This method is
### based on h5read( , start=start, count=count) which retrieves a
### linear subset of the data and is much faster than doing
### h5read( , index=list(seq(start, length.out=count))).
### 'j' must be NULL or an integer vector containing valid col indices. It
### should not be empty.
.linear_extract_nonzero_data_by_col <- function(x, j)
{
    if (is.null(j)) {
        j1 <- 1L
        j2 <- ncol(x)
    } else {
        stopifnot(is.numeric(j), length(j) != 0L)
        j1 <- min(j)
        j2 <- max(j)
    }
    nonzero_data <- .extract_data_from_adjacent_cols(x, j1, j2)
    if (is.null(j))
        return(nonzero_data)
    nonzero_data[match(j, j1:j2)]
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

### 'j' must be NULL or an integer vector containing valid col indices.
### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
.extract_nonzero_data_by_col <-
    function(x, j, method=c("auto", "random", "linear"))
{
    method <- match.arg(method)
    method <- .normarg_method(method, j)
    if (method == "random") {
        .random_extract_nonzero_data_by_col(x, j)
    } else {
        .linear_extract_nonzero_data_by_col(x, j)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_SparseArray_from_TENxMatrixSeed()
###

### Extract sparse data using the "random" method. This method is
### based on h5read( , index=idx) which retrieves an arbitrary/random
### subset of the data.
### 'i' must be NULL or an integer vector containing valid row indices.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
### Both 'i' and 'j' can contain duplicates. Duplicates in 'i' have no effect
### on the output but duplicates in 'j' will produce duplicates in the output.
### Return a SparseArray object.
.random_extract_SparseArray_from_TENxMatrixSeed <- function(x, i, j)
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
    ans_aind <- cbind(row_indices, col_indices, deparse.level=0L)
    ans_nzdata <- .get_data(x@filepath, x@group, idx=idx2)
    SparseArray(dim(x), ans_aind, ans_nzdata, check=FALSE)
}

### Extract sparse data using the "linear" method. This method is
### based on h5read( , start=start, count=count) which retrieves a
### linear subset of the data and is much faster than doing
### h5read( , index=list(seq(start, length.out=count))).
### 'i' (or 'j') must be NULL or an integer vector containing valid
### row (or col) indices. 'j' should not be empty.
### The output is not affected by duplicates in 'i' or 'j'.
### Return a SparseArray object.
.linear_extract_SparseArray_from_TENxMatrixSeed <- function(x, i, j)
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
    sparse_array <- .extract_data_from_adjacent_cols(x, j1, j2, as.sparse=TRUE)

    ## TODO: Wrap the dance below into some sort of multi-dimensional
    ## subsetting helper for SparseArray objects. Don't use "[" for this
    ## because it wouldn't be clear what to do in the 1-dimension case
    ## (we probably should have a "[" method that supports subsetting along
    ## the **length** of the SparseArray object).
    if (is.null(i) && is.null(j))
        return(sparse_array)
    ans_aind <- aind(sparse_array)
    row_indices <- ans_aind[ , 1L]
    col_indices <- ans_aind[ , 2L]
    if (is.null(i)) {
        keep_me <- col_indices %in% j
    } else if (is.null(j)) {
        keep_me <- row_indices %in% i
    } else {
        keep_me <- (row_indices %in% i) & (col_indices %in% j)
    }
    ## TODO: Define and use a "[" method for SparseArray objects for doing
    ## this.
    keep_idx <- which(keep_me)
    sparse_array@aind <- ans_aind[keep_idx, , drop=FALSE]
    sparse_array@nzdata <- sparse_array@nzdata[keep_idx]
    sparse_array
}

### 'i' (or 'j') must be NULL or an integer vector containing valid
### row (or col) indices.
### Duplicates in 'i' are ok and won't affect the output.
### Duplicates in 'j' are ok but might introduce duplicates in the output
### so should be avoided.
### Return a SparseArray object.
.extract_SparseArray_from_TENxMatrixSeed <-
    function(x, i, j, method=c("auto", "random", "linear"))
{
    method <- match.arg(method)
    method <- .normarg_method(method, j)
    if (method == "random") {
        .random_extract_SparseArray_from_TENxMatrixSeed(x, i, j)
    } else {
        .linear_extract_SparseArray_from_TENxMatrixSeed(x, i, j)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

### 'nrow' and 'ncol' must be single integers.
### 'i'and 'j' must be NULLs or integer vectors. If the latter, they must
### be of length 'nrow' and 'ncol' and contain positive values.
### 'ui' and/or 'uj' must be NULLs (if 'i' and/or 'j' are NULLs) or
### integer vectors equal to 'unique(i)' and 'unique(j)', respectively.
### 'sparse_array' must be a SparseArray object.
### If 'i' and/or 'j' is NULL, the values in 'aind(sparse_array)[ , 1]'
### and/or 'aind(sparse_array)[ , 2]' must be >= 1 and <= 'nrow' and/or
### 'ncol', respectively. Otherwise they must be present in 'i' and/or 'j',
### respectively.
.make_submatrix_from_remapped_SparseArray <- function(nrow, ncol,
                                                      i, j, ui, uj,
                                                      sparse_array)
{
    stopifnot(is(sparse_array, "SparseArray"))
    i2ui <- NULL
    if (is.null(i)) {
        umat_nrow <- nrow
    } else {
        sparse_array@aind[ , 1L] <- match(sparse_array@aind[ , 1L], ui)
        umat_nrow <- length(ui)
        if (!identical(i, ui))
            i2ui <- match(i, ui)
    }
    j2uj <- NULL
    if (is.null(j)) {
        umat_ncol <- ncol
    } else {
        sparse_array@aind[ , 2L] <- match(sparse_array@aind[ , 2L], uj)
        umat_ncol <- length(uj)
        if (!identical(j, ui))
            j2uj <- match(j, uj)
    }
    sparse_array@dim <- c(umat_nrow, umat_ncol)
    umat <- sparse2dense(sparse_array)
    if (is.null(i2ui) && is.null(j2uj))
        return(umat)
    DelayedArray:::subset_by_Nindex(umat, list(i2ui, j2uj))
}

.extract_array_from_TENxMatrixSeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)) {
        ## Return an empty matrix.
        data <- .get_data(x@filepath, x@group, idx=integer(0))
        return(array(data, dim=ans_dim))
    }
    ## Passing 'i' and 'j' thru as.integer() should not be necessary
    ## because the list elements in 'index' are supposed to always be
    ## integer vectors (or NULLs). We do it anyway, because, when called in
    ## some special contexts (e.g. testing), extract_array() could receive
    ## an 'index' where the list elements are numeric vectors or integer
    ## vectors with attributes on them (e.g. names). However, in order to
    ## make sure that a test like 'identical(i, ui)' will behave as expected,
    ## we need to make sure that 'i' (and 'ui') are "naked" integer vectors.
    i <- index[[1L]]
    if (is.null(i)) {
        ui <- NULL
    } else {
        i <- as.integer(i)  # make sure 'i' is a "naked" integer vector
        ui <- unique(i)
    }
    j <- index[[2L]]
    if (is.null(j)) {
        uj <- NULL
    } else {
        j <- as.integer(j)  # make sure 'j' is a "naked" integer vector
        uj <- unique(j)
    }
    sparse_array <- .extract_SparseArray_from_TENxMatrixSeed(x, ui, uj)
    .make_submatrix_from_remapped_SparseArray(ans_dim[[1L]], ans_dim[[2L]],
                                              i, j, ui, uj, sparse_array)
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
    stopifnot(is.null(rownames) || length(rownames) == dim[[1L]])
    colnames <- .get_barcodes(filepath, group)
    stopifnot(is.null(colnames) || length(colnames) == dim[[2L]])
    dimnames <- list(rownames, colnames)

    ## col_ranges
    data_len <- h5length(filepath, paste0(group, "/data"))
    indices_len <- h5length(filepath, paste0(group, "/indices"))
    stopifnot(data_len == indices_len)
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
### Taking advantage of sparsity
###

setGeneric("sparsity", signature="x", function(x) standardGeneric("sparsity"))

setMethod("sparsity", "TENxMatrixSeed",
    function(x)
    {
        data_len <- h5length(x@filepath, paste0(x@group, "/data"))
        1 - data_len / length(x)
    }
)

setMethod("sparsity", "TENxMatrix", function(x) sparsity(x@seed))

.read_sparse_block_from_TENxMatrixSeed <- function(x, viewport)
{
    index <- DelayedArray:::makeNindexFromArrayViewport(
                                      viewport,
                                      expand.RangeNSBS=TRUE)
    i <- index[[1L]]
    j <- index[[2L]]
    sparse_array <- .extract_SparseArray_from_TENxMatrixSeed(x, i, j)
    offsets <- start(viewport) - 1L
    sparse_array@aind[ , 1L] <- sparse_array@aind[ , 1L] - offsets[[1L]]
    sparse_array@aind[ , 2L] <- sparse_array@aind[ , 2L] - offsets[[2L]]
    sparse_array@dim <- dim(viewport)
    sparse_array
}

setMethod("read_sparse_block", "TENxMatrixSeed",
    .read_sparse_block_from_TENxMatrixSeed
)

setMethod("read_sparse_block", "TENxMatrix",
    function(x, viewport) read_sparse_block(x@seed, viewport)
)

### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
### Spelling: "nonzero" preferred over "non-zero". See:
###   https://gcc.gnu.org/ml/gcc/2001-10/msg00610.html
setGeneric("extractNonzeroDataByCol", signature="x",
    function(x, j) standardGeneric("extractNonzeroDataByCol")
)

setMethod("extractNonzeroDataByCol", "TENxMatrixSeed",
    function(x, j)
    {
        j <- DelayedArray:::normalizeSingleBracketSubscript2(j, ncol(x),
                                                             colnames(x))
        .extract_nonzero_data_by_col(x, j)
    }
)

setMethod("extractNonzeroDataByCol", "TENxMatrix",
    function(x, j) extractNonzeroDataByCol(x@seed, j)
)

