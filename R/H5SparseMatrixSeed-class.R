### =========================================================================
### H5SparseMatrixSeed objects
### -------------------------------------------------------------------------


setClass("H5SparseMatrixSeed",
    contains="Array",
    representation(
        "VIRTUAL",

    ## --------------------- user supplied slots ---------------------

        ## Absolute path to the HDF5 file so the object won't break when
        ## the user changes the working directory (e.g. with 'setwd()').
        filepath="character",

        ## Name of the group in the HDF5 file where the sparse matrix is
        ## stored.
        group="character",

        ## If 'paste0(group, "/data")' is a group, name of a dataset in
        ## that group. Otherwise, must be set to NULL.
        subdata="character_OR_NULL",

    ## ---------------- automatically populated slots ----------------

        dim="integer",

        ## Can't use an IRanges object for this at the moment because IRanges
        ## objects don't support large integer start/end values yet.
        indptr_ranges="data.frame",

    ## ------------- populated by specialized subclasses -------------

        dimnames="list"
    ),
    prototype(
        dimnames=list(NULL, NULL)
    )
)

.get_data_name <- function(subdata, group=NULL)
{
    name <- "data"
    if (!is.null(subdata))
        name <- paste0(name, "/", subdata)
    if (!is.null(group))
        name <- paste0(group, "/", name)
    name
}

setClass("CSC_H5SparseMatrixSeed", contains="H5SparseMatrixSeed")
setClass("CSR_H5SparseMatrixSeed", contains="H5SparseMatrixSeed")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transposition
###

### S3/S4 combo for t.CSC_H5SparseMatrixSeed
t.CSC_H5SparseMatrixSeed <- function(x)
{
    x@dim <- rev(x@dim)
    class(x) <- "CSR_H5SparseMatrixSeed"
    x
}
setMethod("t", "CSC_H5SparseMatrixSeed", t.CSC_H5SparseMatrixSeed)

### S3/S4 combo for t.CSR_H5SparseMatrixSeed
t.CSR_H5SparseMatrixSeed <- function(x)
{
    x@dim <- rev(x@dim)
    class(x) <- "CSC_H5SparseMatrixSeed"
    x
}
setMethod("t", "CSR_H5SparseMatrixSeed", t.CSR_H5SparseMatrixSeed)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### path() getter/setter
###

### Does NOT access the file.
setMethod("path", "H5SparseMatrixSeed", function(object) object@filepath)

### Just a placeholder for now. Doesn't actually allow changing the path of
### the object yet.
setReplaceMethod("path", "H5SparseMatrixSeed",
    function(object, value)
    {
        new_filepath <- normarg_h5_filepath(value, what1="the supplied path",
                                                   what2="the sparse matrix")
        old_filepath <- path(object)
        if (new_filepath != old_filepath)
            stop(wmsg("changing the path of a ", class(object), " object ",
                      "is not supported yet"))
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() and dimnames() getters
###
### They access the slot, not the file.
###

setMethod("dim", "H5SparseMatrixSeed", function(x) x@dim)

setMethod("dimnames", "H5SparseMatrixSeed",
    function(x) S4Arrays:::simplify_NULL_dimnames(x@dimnames)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level internal h5sparse data readers
###

### All the h5sparse components are monodimensional.
read_h5sparse_component <- function(filepath, group, name,
                                    start=NULL, count=NULL, as.integer=FALSE)
{
    name <- paste0(group, "/", name)
    if (!is.null(start))
        start <- list(start)
    if (!is.null(count))
        count <- list(count)
    as.vector(h5mread(filepath, name, starts=start, counts=count,
                                      as.integer=as.integer))
}

### Returns a numeric vector (integer or double).
.read_h5sparse_dim <- function(filepath, group)
{
    if (h5exists(filepath, paste0(group, "/shape"))) {
        ## 10x format
        return(read_h5sparse_component(filepath, group, "shape"))
    }
    ## h5ad format
    h5attrs <- h5readAttributes(filepath, group)
    shape <- h5attrs$shape
    if (is.null(shape))
        shape <- h5attrs$h5sparse_shape
    if (is.null(shape))
        stop(wmsg("Group \"", group, "\" in HDF5 file \"", filepath,"\" ",
                  "contains no 'shape' dataset and has no 'shape' ",
                  "or 'h5sparse_shape' attribute. As a consequence, the ",
                  "dimensions of the sparse matrix can't be determined."))
    ## We pass 'shape' thru as.vector() to drop its class attribute in case
    ## it's an array.
    rev(as.vector(shape))
}

.read_h5sparse_format <- function(filepath, group)
{
    if (h5exists(filepath, paste0(group, "/shape"))) {
        ## 10x format
        return("csr")
    }
    ## h5ad format
    h5attrs <- h5readAttributes(filepath, group)
    h5sparse_format <- h5attrs[["encoding-type"]]
    if (is.null(h5sparse_format))
        h5sparse_format <- h5attrs[["h5sparse_format"]]
    if (is.null(h5sparse_format))
        return("csr")
    ans <- tolower(substr(h5sparse_format, 1L, 3L))
    if (!(ans %in% c("csr", "csc")))
        stop(wmsg("sparse matrix in group \"", group, "\" in HDF5 ",
                  "file \"", filepath,"\" is stored in unsupported ",
                  "format \"", h5sparse_format, "\""))
    ans
}

.read_h5sparse_indptr <- function(filepath, group)
    read_h5sparse_component(filepath, group, "indptr")

.read_h5sparse_data <-
    function(filepath, group, subdata, start=NULL, count=NULL)
{
    name <- .get_data_name(subdata)
    read_h5sparse_component(filepath, group, name, start=start, count=count)
}

### The row (or column) indices stored in HDF5 dataset "indices" are 0-based
### but we return them 1-based.
.read_h5sparse_indices <- function(filepath, group, start=NULL, count=NULL)
    read_h5sparse_component(filepath, group, "indices",
                            start=start, count=count, as.integer=TRUE) + 1L


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.check_group <- function(filepath, group)
{
    if (!h5exists(filepath, group))
        stop(wmsg("HDF5 group \"", group, "\" does not exist ",
                  "in this HDF5 file"))
    if (h5isdataset(filepath, group)) {
        is_h5ad_X_or_layer <- group == "/X" ||
                              substr(group, 1L, 8L) == "/layers/"
        msg1 <- c("\"", group, "\" is an HDF5 dataset, not an HDF5 group, ",
                  "so it looks like the matrix that you are trying to ",
                  "access is not stored in a sparse format. Please ",
                  "consider using the ")
        if (is_h5ad_X_or_layer) {
            msg2 <- c("H5ADMatrix() constructor if you are trying ",
                      "to access the central matrix of an h5ad file. ",
                      "Otherwise, use the HDF5Array() constructor.")
        } else {
            msg2 <- "HDF5Array() constructor to access this dataset."
        }
        stop(wmsg(msg1, msg2))
    }
    if (!h5isgroup(filepath, group))
        stop(wmsg("HDF5 object \"", group, "\" is not a group"))
}

.check_data_and_subdata <- function(filepath, group, subdata)
{
    data_fullname <- paste0(group, "/data")
    if (!h5exists(filepath, data_fullname))
        stop(wmsg("HDF5 object \"", data_fullname, "\" does not ",
                  "exist in this HDF5 file. Are you sure that HDF5 ",
                  "group \"", group, "\" contains a sparse matrix ",
                  "stored in CSR/CSC/Yale format?"))
    if (is.null(subdata)) {
        if (h5isgroup(filepath, data_fullname))
            stop(wmsg("\"", data_fullname, "\" is an HDF5 group, not an ",
                      "HDF5 dataset. Please use the 'subdata' argument to ",
                      "specify the name of the dataset in this group that ",
                      "contains the matrix data."))
        if (!h5isdataset(filepath, data_fullname))
            stop(wmsg("HDF5 object \"", data_fullname, "\" is not a dataset."))
    } else {
        if (!isSingleString(subdata) || subdata == "")
            stop(wmsg("'subdata' must be NULL or a single non-empty string"))
        if (h5isdataset(filepath, data_fullname))
            stop(wmsg("\"", data_fullname, "\" is an HDF5 dataset, not an ",
                      "HDF5 group. Please note that the 'subdata' argument ",
                      "can be used only when it's a group."))
        if (!h5isgroup(filepath, data_fullname))
            stop(wmsg("HDF5 object \"", data_fullname, "\" is not a group."))
        subdata_fullname <- .get_data_name(subdata, group)
        if (!h5exists(filepath, subdata_fullname))
            stop(wmsg("HDF5 object \"", subdata_fullname, "\" does not ",
                      "exist in this HDF5 file."))
        if (!h5isdataset(filepath, subdata_fullname))
            stop(wmsg("HDF5 object \"", subdata_fullname, "\" is ",
                      "not a dataset."))
    }
}

### Returns an H5SparseMatrixSeed derivative (can be either a
### CSC_H5SparseMatrixSeed or CSR_H5SparseMatrixSeed object).
H5SparseMatrixSeed <- function(filepath, group, subdata=NULL)
{
    ## Check 'filepath', 'group', and 'subdata'.
    filepath <- normarg_h5_filepath(filepath, what2="the sparse matrix")
    group <- normarg_h5_name(group, what1="'group'",
                                    what2="the name of the group",
                                    what3=" that stores the sparse matrix")
    .check_group(filepath, group)
    .check_data_and_subdata(filepath, group, subdata)

    ## Get matrix dimensions.
    dim <- .read_h5sparse_dim(filepath, group)
    stopifnot(length(dim) == 2L)
    dim <- dim_as_integer(dim, filepath, group, what="sparse matrix")

    ## Get sparse format ("csc" or "csr").
    h5sparse_format <- .read_h5sparse_format(filepath, group)
    if (h5sparse_format == "csr") {
        expected_indptr_length <- dim[[2L]] + 1L
        ## Because R has the notions of rows and columns flipped w.r.t.
        ## HDF5, "compressed sparse row" at the HDF5 level translates
        ## into "compressed sparse column" at the R level.
        ans_class <- "CSC_H5SparseMatrixSeed"
    } else {
        expected_indptr_length <- dim[[1L]] + 1L
        ## Because R has the notions of rows and columns flipped w.r.t.
        ## HDF5, "compressed sparse column" at the HDF5 level translates
        ## into "compressed sparse row" at the R level.
        ans_class <- "CSR_H5SparseMatrixSeed"
    }

    ## Get 'indptr_ranges'.
    data_len <- h5length(filepath, .get_data_name(subdata, group))
    indices_len <- h5length(filepath, paste0(group, "/indices"))
    stopifnot(data_len == indices_len)
    indptr <- .read_h5sparse_indptr(filepath, group)
    stopifnot(length(indptr) == expected_indptr_length,
              indptr[[1L]] == 0L,
              indptr[[length(indptr)]] == data_len)
    indptr_ranges <- data.frame(start=indptr[-length(indptr)] + 1,
                                width=as.integer(diff(indptr)))

    new2(ans_class, filepath=filepath, group=group,
                    dim=dim, indptr_ranges=indptr_ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_data_indices_by_col()
###

### base::sequence() does not properly handle a 'from' that is >
### .Machine$integer.max so we implement a variant that does. Note that
### the 2nd argument of this variant is 'offset' instead of 'from' ('offset'
### being the same as 'from' - 1).
### TODO: Does .sequence2() work if sum(lengths) is > .Machine$integer.max?
.sequence2 <- function(lengths, offset=0)
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
    indptr_ranges <- S4Vectors:::extract_data_frame_rows(x@indptr_ranges, j)
    start2 <- indptr_ranges[ , "start"]
    width2 <- indptr_ranges[ , "width"]
    idx2 <- .sequence2(width2, offset=start2 - 1L)
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
    start <- x@indptr_ranges[j1, "start"]
    count_per_col <- x@indptr_ranges[j12, "width"]
    count <- sum(count_per_col)
    ans_nzdata <- .read_h5sparse_data(x@filepath, x@group, x@subdata,
                                      start=start, count=count)
    if (!as.sparse)
        return(relist(ans_nzdata, PartitioningByWidth(count_per_col)))
    row_indices <- .read_h5sparse_indices(x@filepath, x@group,
                                          start=start, count=count)
    col_indices <- rep.int(j12, count_per_col)
    ans_nzindex <- cbind(row_indices, col_indices, deparse.level=0L)
    SparseArraySeed(dim(x), ans_nzindex, ans_nzdata, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .load_sparse_data()
###
### This internal generic is the workhorse behind the extract_array(),
### extract_sparse_array(), and read_sparse_block() methods for
### H5SparseMatrixSeed objects.
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
.load_random_csc_sparse_data <- function(x, i, j)
{
    stopifnot(is.null(i) || is.numeric(i), is.numeric(j))
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    row_indices <- .read_h5sparse_indices(x@filepath, x@group, start=idx2)
    col_indices <- rep.int(j, lengths(data_indices))
    if (!is.null(i)) {
        keep_idx <- which(row_indices %in% i)
        idx2 <- idx2[keep_idx]
        row_indices <- row_indices[keep_idx]
        col_indices <- col_indices[keep_idx]
    }
    ans_nzindex <- cbind(row_indices, col_indices, deparse.level=0L)
    ans_nzdata <- .read_h5sparse_data(x@filepath, x@group, x@subdata,
                                      start=idx2)
    SparseArraySeed(dim(x), ans_nzindex, ans_nzdata, check=FALSE)
}

### Load sparse data using the "linear" method.
### This method is based on h5mread( , starts=list(start), counts=list(count))
### which retrieves a linear subset of the data and should be more efficient
### than doing h5mread( , starts=list(seq(start, length.out=count))).
### 'j' must be NULL or a non-empty integer vector containing valid
### col indices. The output is not affected by duplicates in 'j'.
### Return a SparseArraySeed object.
.load_linear_csc_sparse_data <- function(x, j)
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

### Duplicates in 'index[[1]]' are ok and won't affect the output.
### Duplicates in 'index[[2]]' are ok but might introduce duplicates
### in the output so should be avoided.
### Return a SparseArraySeed object.
.load_csc_sparse_data <- function(x, index, method)
{
    i <- index[[1L]]
    j <- index[[2L]]
    method <- .normarg_method(method, j)
    if (method == "random") {
        ans <- .load_random_csc_sparse_data(x, i, j)
    } else {
        ans <- .load_linear_csc_sparse_data(x, j)
    }
    ans
}

### Return a SparseArraySeed object.
setGeneric(".load_sparse_data", signature="x",
    function(x, index, method=c("auto", "random", "linear"))
        standardGeneric(".load_sparse_data")
)

setMethod(".load_sparse_data", "CSC_H5SparseMatrixSeed",
    function(x, index, method=c("auto", "random", "linear"))
    {
        method <- match.arg(method)
        .load_csc_sparse_data(x, index, method)
    }
)

setMethod(".load_sparse_data", "CSR_H5SparseMatrixSeed",
    function(x, index, method=c("auto", "random", "linear"))
    {
        method <- match.arg(method)
        t(.load_csc_sparse_data(t(x), rev(index), method))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.extract_array_from_H5SparseMatrixSeed <- function(x, index)
{
    ans_dim <- S4Arrays:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)) {
        ## Return an empty matrix.
        data0 <- .read_h5sparse_data(x@filepath, x@group, x@subdata,
                                     start=integer(0))
        return(array(data0, dim=ans_dim))
    }
    sas <- .load_sparse_data(x, index)  # I/O
    extract_array(sas, index)  # in-memory
}

setMethod("extract_array", "H5SparseMatrixSeed",
    .extract_array_from_H5SparseMatrixSeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkdim() getter
###

### Does NOT access the file.
setMethod("chunkdim", "CSC_H5SparseMatrixSeed", function(x) c(nrow(x), 1L))

setMethod("chunkdim", "CSR_H5SparseMatrixSeed", function(x) c(1L, ncol(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Taking advantage of sparsity
###

setMethod("sparsity", "H5SparseMatrixSeed",
    function(x)
    {
        data_len <- h5length(x@filepath, .get_data_name(x@subdata, x@group))
        1 - data_len / length(x)
    }
)

### This is about **structural** sparsity, not about quantitative sparsity
### measured by sparsity().
setMethod("is_sparse", "H5SparseMatrixSeed", function(x) TRUE)

.extract_sparse_array_from_H5SparseMatrixSeed <- function(x, index)
{
    sas <- .load_sparse_data(x, index)  # I/O
    extract_sparse_array(sas, index)  # in-memory
}

setMethod("extract_sparse_array", "H5SparseMatrixSeed",
    .extract_sparse_array_from_H5SparseMatrixSeed
)

### The default read_sparse_block() method defined in DelayedArray would
### work just fine on an H5SparseMatrixSeed derivative (thanks to the
### extract_sparse_array() method for H5SparseMatrixSeed objects defined
### above), but we overwrite it with the method below which should be
### slightly more efficient. That's because the method below calls
### read_sparse_block() on the SparseArraySeed object returned by
### .load_sparse_data(), and this should be faster than calling
### extract_sparse_array() on the same object (which is what the
### extract_sparse_array() method for H5SparseMatrixSeed objects would
### be doing when called by the default read_sparse_block() method).
### Not sure the difference is actually significant enough for this extra
### method to be worth it though, because time is really dominated by I/O
### here, that is, by the call to .load_sparse_data().
.read_sparse_block_from_H5SparseMatrixSeed <- function(x, viewport)
{
    index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
    sas <- .load_sparse_data(x, index)  # I/O
    ## Unlike the extract_sparse_array() method for H5SparseMatrixSeed
    ## objects defined above, we use read_sparse_block() here, which should
    ## be faster than using extract_sparse_array().
    read_sparse_block(sas, viewport)  # in-memory
}

setMethod("read_sparse_block", "H5SparseMatrixSeed",
    .read_sparse_block_from_H5SparseMatrixSeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractNonzeroDataByCol() and extractNonzeroDataByRow()
###
### Spelling: "nonzero" preferred over "non-zero". See:
###   https://gcc.gnu.org/ml/gcc/2001-10/msg00610.html
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
    data <- .read_h5sparse_data(x@filepath, x@group, x@subdata, start=idx2)
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

### 'j' must be NULL or an integer vector containing valid col indices.
### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
.extract_nonzero_csc_sparse_data_by_col <-
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

### Return a NumericList or IntegerList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
setGeneric("extractNonzeroDataByCol", signature="x",
    function(x, j) standardGeneric("extractNonzeroDataByCol")
)

setMethod("extractNonzeroDataByCol", "CSC_H5SparseMatrixSeed",
    function(x, j)
    {
        j <- S4Arrays:::normalizeSingleBracketSubscript2(j, ncol(x),
                                                         colnames(x))
        .extract_nonzero_csc_sparse_data_by_col(x, j)
    }
)

### Return a NumericList or IntegerList object parallel to 'i' i.e. with
### one list element per row index in 'i'.
setGeneric("extractNonzeroDataByRow", signature="x",
    function(x, i) standardGeneric("extractNonzeroDataByRow")
)

setMethod("extractNonzeroDataByRow", "CSR_H5SparseMatrixSeed",
    function(x, i)
    {
        i <- S4Arrays:::normalizeSingleBracketSubscript2(i, nrow(x),
                                                         rownames(x))
        .extract_nonzero_csc_sparse_data_by_col(t(x), i)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to dgCMatrix
###

.from_CSC_H5SparseMatrixSeed_to_dgCMatrix <- function(from)
{
    row_indices <- .read_h5sparse_indices(from@filepath, from@group)
    indptr <- .read_h5sparse_indptr(from@filepath, from@group)
    data <- .read_h5sparse_data(from@filepath, from@group, from@subdata)
    sparseMatrix(i=row_indices, p=indptr, x=data, dims=dim(from),
                 dimnames=dimnames(from))
}

setAs("CSC_H5SparseMatrixSeed", "dgCMatrix",
    .from_CSC_H5SparseMatrixSeed_to_dgCMatrix
)
setAs("CSC_H5SparseMatrixSeed", "sparseMatrix",
    .from_CSC_H5SparseMatrixSeed_to_dgCMatrix
)

.from_CSR_H5SparseMatrixSeed_to_dgCMatrix <- function(from)
{
    col_indices <- .read_h5sparse_indices(from@filepath, from@group)
    indptr <- .read_h5sparse_indptr(from@filepath, from@group)
    data <- .read_h5sparse_data(from@filepath, from@group, from@subdata)
    sparseMatrix(j=col_indices, p=indptr, x=data, dims=dim(from),
                 dimnames=dimnames(from))
}

setAs("CSR_H5SparseMatrixSeed", "dgCMatrix",
    .from_CSR_H5SparseMatrixSeed_to_dgCMatrix
)
setAs("CSR_H5SparseMatrixSeed", "sparseMatrix",
    .from_CSR_H5SparseMatrixSeed_to_dgCMatrix
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod("show", "H5SparseMatrixSeed",
    function(object)
    {
        cat(S4Arrays:::array_as_one_line_summary(object), ":\n", sep="")
        cat("# dirname: ", dirname(object), "\n", sep="")
        cat("# basename: ", basename(object), "\n", sep="")
        cat("# group: ", object@group, "\n", sep="")
    }
)

