### =========================================================================
### TENxMatrixSeed objects
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
.read_tenx_component <- function(filepath, group, name, start=NULL, count=NULL,
                                 as.integer=FALSE)
{
    name <- paste0(group, "/", name)
    if (!is.null(start))
        start <- list(start)
    if (!is.null(count))
        count <- list(count)
    as.vector(h5mread(filepath, name, starts=start, counts=count,
                      as.integer=as.integer))
}

### Return the dimensions of the matrix.
.get_tenx_shape <- function(filepath, group)
    .read_tenx_component(filepath, group, "shape")

.get_tenx_indptr <- function(filepath, group)
    .read_tenx_component(filepath, group, "indptr")

.get_tenx_data <- function(filepath, group, start=NULL, count=NULL)
    .read_tenx_component(filepath, group, "data", start=start, count=count)

### The row indices in the HDF5 file are 0-based but we return them 1-based.
.get_tenx_row_indices <- function(filepath, group, start=NULL, count=NULL)
    .read_tenx_component(filepath, group, "indices", start=start, count=count,
                         as.integer=TRUE) + 1L

### Return the rownames of the matrix.
.get_tenx_genes <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/genes")))
        return(NULL)
    .read_tenx_component(filepath, group, "genes")
}

### Currently unused.
.get_tenx_gene_names <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/gene_names")))
        return(NULL)
    .read_tenx_component(filepath, group, "gene_names")
}

### Return the colnames of the matrix.
.get_tenx_barcodes <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/barcodes")))
        return(NULL)
    .read_tenx_component(filepath, group, "barcodes")
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
### If 'as.sparse=TRUE', return a SparseArraySeed object.
.extract_data_from_adjacent_cols <- function(x, j1, j2, as.sparse=FALSE)
{
    j12 <- j1:j2
    start <- x@col_ranges[j1, "start"]
    count_per_col <- x@col_ranges[j12, "width"]
    count <- sum(count_per_col)
    ans_nzdata <- .get_tenx_data(x@filepath, x@group, start=start, count=count)
    if (!as.sparse)
        return(relist(ans_nzdata, PartitioningByWidth(count_per_col)))
    row_indices <- .get_tenx_row_indices(x@filepath, x@group,
                                         start=start, count=count)
    col_indices <- rep.int(j12, count_per_col)
    ans_nzindex <- cbind(row_indices, col_indices, deparse.level=0L)
    SparseArraySeed(dim(x), ans_nzindex, ans_nzdata, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_nonzero_data_by_col()
###

### Extract nonzero data using the "random" method.
### This method is based on h5mread( , starts=list(start)) which retrieves
### an arbitrary/random subset of the data.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
.random_extract_nonzero_data_by_col <- function(x, j)
{
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    data <- .get_tenx_data(x@filepath, x@group, start=idx2)
    relist(data, data_indices)
}

### Extract nonzero data using the "linear" method.
### This method is based on h5mread( , starts=list(start), counts=list(count))
### which retrieves a linear subset of the data and should be more efficient
### than doing h5mread( , starts=list(seq(start, length.out=count))).
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
### .load_SparseArraySeed_from_TENxMatrixSeed()
###

### Load sparse data using the "random" method.
### This method is based on h5mread( , starts=list(start)) which retrieves
### an arbitrary/random subset of the data.
### 'i' must be NULL or an integer vector containing valid row indices.
### 'j' must be an integer vector containing valid col indices. It cannot
### be NULL.
### Both 'i' and 'j' can contain duplicates. Duplicates in 'i' have no effect
### on the output but duplicates in 'j' will produce duplicates in the output.
### Return a SparseArraySeed object.
.random_load_SparseArraySeed_from_TENxMatrixSeed <- function(x, i, j)
{
    stopifnot(is.null(i) || is.numeric(i), is.numeric(j))
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    row_indices <- .get_tenx_row_indices(x@filepath, x@group, start=idx2)
    col_indices <- rep.int(j, lengths(data_indices))
    if (!is.null(i)) {
        keep_idx <- which(row_indices %in% i)
        idx2 <- idx2[keep_idx]
        row_indices <- row_indices[keep_idx]
        col_indices <- col_indices[keep_idx]
    }
    ans_nzindex <- cbind(row_indices, col_indices, deparse.level=0L)
    ans_nzdata <- .get_tenx_data(x@filepath, x@group, start=idx2)
    SparseArraySeed(dim(x), ans_nzindex, ans_nzdata, check=FALSE)
}

### Load sparse data using the "linear" method.
### This method is based on h5mread( , starts=list(start), counts=list(count))
### which retrieves a linear subset of the data and should be more efficient
### than doing h5mread( , starts=list(seq(start, length.out=count))).
### 'j' must be NULL or a non-empty integer vector containing valid
### col indices. The output is not affected by duplicates in 'j'.
### Return a SparseArraySeed object.
.linear_load_SparseArraySeed_from_TENxMatrixSeed <- function(x, j)
{
    if (is.null(j)) {
        j1 <- 1L
        j2 <- ncol(x)
    } else {
        stopifnot(is.numeric(j), length(j) != 0L)
        j1 <- min(j)
        j2 <- max(j)
    }
    .extract_data_from_adjacent_cols(x, j1, j2, as.sparse=TRUE)
}

### Duplicates in 'index[[1]]' are ok and won't affect the output.
### Duplicates in 'index[[2]]' are ok but might introduce duplicates
### in the output so should be avoided.
### Return a SparseArraySeed object.
.load_SparseArraySeed_from_TENxMatrixSeed <-
    function(x, index, method=c("auto", "random", "linear"))
{
    i <- index[[1L]]
    j <- index[[2L]]
    method <- match.arg(method)
    method <- .normarg_method(method, j)
    if (method == "random") {
        .random_load_SparseArraySeed_from_TENxMatrixSeed(x, i, j)
    } else {
        .linear_load_SparseArraySeed_from_TENxMatrixSeed(x, j)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.extract_array_from_TENxMatrixSeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)) {
        ## Return an empty matrix.
        data0 <- .get_tenx_data(x@filepath, x@group, start=integer(0))
        return(array(data0, dim=ans_dim))
    }
    sas <- .load_SparseArraySeed_from_TENxMatrixSeed(x, index)  # I/O
    extract_array(sas, index)  # in-memory
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
### Constructor
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
    dim <- .get_tenx_shape(filepath, group)
    stopifnot(length(dim) == 2L)

    ## dimnames
    rownames <- .get_tenx_genes(filepath, group)
    stopifnot(is.null(rownames) || length(rownames) == dim[[1L]])
    colnames <- .get_tenx_barcodes(filepath, group)
    stopifnot(is.null(colnames) || length(colnames) == dim[[2L]])
    dimnames <- list(rownames, colnames)

    ## col_ranges
    data_len <- h5length(filepath, paste0(group, "/data"))
    indices_len <- h5length(filepath, paste0(group, "/indices"))
    stopifnot(data_len == indices_len)
    indptr <- .get_tenx_indptr(filepath, group)
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
### Taking advantage of sparsity
###

setMethod("sparsity", "TENxMatrixSeed",
    function(x)
    {
        data_len <- h5length(x@filepath, paste0(x@group, "/data"))
        1 - data_len / length(x)
    }
)

### This is about **structural** sparsity, not about quantitative sparsity
### measured by sparsity().
setMethod("is_sparse", "TENxMatrixSeed", function(x) TRUE)

.extract_sparse_array_from_TENxMatrixSeed <- function(x, index)
{
    sas <- .load_SparseArraySeed_from_TENxMatrixSeed(x, index)  # I/O
    extract_sparse_array(sas, index)  # in-memory
}

setMethod("extract_sparse_array", "TENxMatrixSeed",
    .extract_sparse_array_from_TENxMatrixSeed
)

### The default "read_sparse_block" method defined in DelayedArray would work
### just fine on a TENxMatrixSeed object (thanks to the "extract_sparse_array"
### method for TENxMatrixSeed objects defined above), but we overwrite it with
### the method below which is slightly more efficient. That's because the
### method below calls read_sparse_block() on the SparseArraySeed object
### returned by .load_SparseArraySeed_from_TENxMatrixSeed() and this is
### faster than calling extract_sparse_array() on the same object (which
### is what the "extract_sparse_array" method for TENxMatrixSeed would do
### when called by the default "read_sparse_block" method).
### Not sure the difference is significant enough for this extra method to
### be worth it though, because, time is really dominated by I/O here, that
### is, by the call to .load_SparseArraySeed_from_TENxMatrixSeed().
.read_sparse_block_from_TENxMatrixSeed <- function(x, viewport)
{
    index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
    sas <- .load_SparseArraySeed_from_TENxMatrixSeed(x, index)  # I/O
    ## Unlike the "extract_sparse_array" method for TENxMatrixSeed objects
    ## defined above, we use read_sparse_block() here, which is faster than
    ## using extract_sparse_array().
    read_sparse_block(sas, viewport)  # in-memory
}

setMethod("read_sparse_block", "TENxMatrixSeed",
    .read_sparse_block_from_TENxMatrixSeed
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to dgCMatrix
###

.from_TENxMatrixSeed_to_dgCMatrix <- function(from)
{
    row_indices <- .get_tenx_row_indices(from@filepath, from@group)
    indptr <- .get_tenx_indptr(from@filepath, from@group)
    data <- .get_tenx_data(from@filepath, from@group)
    sparseMatrix(i=row_indices, p=indptr, x=data, dims=dim(from),
                 dimnames=dimnames(from))
}
setAs("TENxMatrixSeed", "dgCMatrix", .from_TENxMatrixSeed_to_dgCMatrix)
setAs("TENxMatrixSeed", "sparseMatrix", .from_TENxMatrixSeed_to_dgCMatrix)

