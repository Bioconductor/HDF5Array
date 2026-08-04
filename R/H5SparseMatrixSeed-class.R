### =========================================================================
### H5SparseMatrixSeed objects
### -------------------------------------------------------------------------


setClass("H5SparseMatrixSeed",
    contains=c("Array", "OutOfMemoryObject"),
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
    x@dimnames <- rev(x@dimnames)
    class(x) <- class(new("CSR_H5SparseMatrixSeed"))
    x
}
setMethod("t", "CSC_H5SparseMatrixSeed", t.CSC_H5SparseMatrixSeed)

### S3/S4 combo for t.CSR_H5SparseMatrixSeed
t.CSR_H5SparseMatrixSeed <- function(x)
{
    x@dim <- rev(x@dim)
    x@dimnames <- rev(x@dimnames)
    class(x) <- class(new("CSC_H5SparseMatrixSeed"))
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
        new_filepath <- h5mread:::normarg_h5_filepath(value,
                                             what1="the supplied path",
                                             what2="the sparse matrix")
        old_filepath <- path(object)
        if (new_filepath != old_filepath)
            stop(wmsg("changing the path of a ", class(object), " object ",
                      "is not supported yet"))
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### group() getter
###

setMethod("group", "H5SparseMatrixSeed",
    function(object)
    {
        group <- object@group
        if (!startsWith(group, "/"))
            group <- paste0("/", group)
        group
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
### chunkdim() getter
###

### Does NOT access the file.
setMethod("chunkdim", "CSC_H5SparseMatrixSeed",
    function(x) c(nrow(x), min(ncol(x), 1L))
)

setMethod("chunkdim", "CSR_H5SparseMatrixSeed",
    function(x) c(min(nrow(x), 1L), ncol(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse() and nzcount() methods
###

### This is about **structural** sparsity, not about quantitative sparsity
### measured by sparsity().
setMethod("is_sparse", "H5SparseMatrixSeed", function(x) TRUE)

setMethod("nzcount", "H5SparseMatrixSeed",
    function(x) h5length(x@filepath, .get_data_name(x@subdata, x@group))
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
    h5mread(filepath, name, starts=start, counts=count,
            as.vector=TRUE, as.integer=as.integer)
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

.read_h5sparse_layout <- function(filepath, group)
{
    if (h5exists(filepath, paste0(group, "/shape"))) {
        ## 10x format
        return("csr")
    }
    ## h5ad format
    h5attrs <- h5readAttributes(filepath, group)
    h5sparse_layout <- h5attrs[["encoding-type"]]
    if (is.null(h5sparse_layout))
        h5sparse_layout <- h5attrs[["h5sparse_format"]]
    if (is.null(h5sparse_layout))
        return("csr")
    ans <- tolower(substr(h5sparse_layout, 1L, 3L))
    if (!(ans %in% c("csr", "csc")))
        stop(wmsg("sparse matrix in group \"", group, "\" in HDF5 ",
                  "file \"", filepath,"\" is stored in unsupported ",
                  "layout \"", h5sparse_layout, "\""))
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
### and we return them as such.
.read_h5sparse_indices <- function(filepath, group, start=NULL, count=NULL)
    read_h5sparse_component(filepath, group, "indices",
                            start=start, count=count, as.integer=TRUE)


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
                  "stored in CSR/CSC/Yale layout?"))
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

.get_sparse_matrix_dim <- function(filepath, group, dim=NULL)
{
    if (is.null(dim)) {
        dim <- .read_h5sparse_dim(filepath, group)
        stopifnot(length(dim) == 2L)
        return(dim_as_integer(dim, filepath, group, what="sparse matrix"))
    }
    ## Check user-supplied 'dim'.
    if (!is.numeric(dim) || length(dim) != 2L || anyNA(dim))
        stop(wmsg("supplied 'dim' must be an integer vector ",
                  "of length 2 with no NAs"))
    if (!is.integer(dim)) {
        if (any(dim > .Machine$integer.max))
            stop(wmsg("supplied dimensions are too big (all dimensions ",
                      "must be <= '.Machine$integer.max' (= 2^31 - 1))"))
        dim <- as.integer(dim)
    }
    if (any(dim < 0L))
        stop(wmsg("supplied 'dim' cannot contain negative values"))
    dim
}

### Must return "CSC" or "CSR".
.get_sparse_matrix_layout <- function(filepath, group, sparse.layout=NULL)
{
    if (is.null(sparse.layout)) {
        h5sparse_layout <- .read_h5sparse_layout(filepath, group)
        ## Layout in R will be transposed w.r.t. layout used in h5 file.
        ans <- switch(h5sparse_layout, `csr`="CSC", `csc`="CSR",
                      stop(wmsg("unsupported 'h5sparse_layout': ",
                                h5sparse_layout)))
        return(ans)
    }
    ## Check user-supplied 'sparse.layout'.
    if (!isSingleString(sparse.layout))
        stop(wmsg("'sparse.layout' must be a single string"))
    ans <- toupper(sparse.layout)
    if (!(ans %in% c("CSC", "CSR")))
        stop(wmsg("'sparse.layout' must be either \"CSC\" or \"CSR\""))
    ans
}

### Returns an H5SparseMatrixSeed derivative (can be either a
### CSC_H5SparseMatrixSeed or CSR_H5SparseMatrixSeed object).
H5SparseMatrixSeed <- function(filepath, group, subdata=NULL,
                               dim=NULL, sparse.layout=NULL)
{
    ## Check 'filepath', 'group', and 'subdata'.
    filepath <- h5mread:::normarg_h5_filepath(filepath,
                                  what2="the sparse matrix")
    group <- h5mread:::normarg_h5_name(group,
                                  what1="'group'",
                                  what2="the name of the group",
                                  what3=" that stores the sparse matrix")
    .check_group(filepath, group)
    .check_data_and_subdata(filepath, group, subdata)

    ## Get matrix dimensions.
    dim <- .get_sparse_matrix_dim(filepath, group, dim=dim)

    ## Get sparse layout to use ("CSC" or "CSR").
    ## Note that R has the notions of rows and columns flipped w.r.t.
    ## HDF5 so:
    ## - "compressed sparse row" at the HDF5 level translates
    ##   into "compressed sparse column" at the R level,
    ## - "compressed sparse column" at the HDF5 level translates
    ##   into "compressed sparse row" at the R level.
    layout <- .get_sparse_matrix_layout(filepath, group,
                                        sparse.layout=sparse.layout)
    if (layout == "CSC") {
        expected_indptr_len <- dim[[2L]] + 1L
        ans_class <- "CSC_H5SparseMatrixSeed"
    } else {
        expected_indptr_len <- dim[[1L]] + 1L
        ans_class <- "CSR_H5SparseMatrixSeed"
    }

    ## Get 'indptr_ranges'.
    nzcount <- h5length(filepath, .get_data_name(subdata, group))
    indices_len <- h5length(filepath, paste0(group, "/indices"))
    stopifnot(indices_len == nzcount)
    indptr <- .read_h5sparse_indptr(filepath, group)
    stopifnot(length(indptr) == expected_indptr_len,
              indptr[[1L]] == 0L,
              indptr[[length(indptr)]] == nzcount)
    indptr_ranges <- data.frame(start=indptr[-length(indptr)] + 1,
                                width=as.integer(diff(indptr)))

    new2(ans_class, filepath=filepath, group=group,
                    dim=dim, indptr_ranges=indptr_ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .load_CSC_H5SparseMatrixSeed
###
### Loads whole CSC_H5SparseMatrixSeed object 'x' into memory as an
### SVT_SparseMatrix object, or only selected columns if 'j' is specified.
### This is the workhorse behind the extract_sparse_array(), extract_array(),
### and read_block_as_sparse() methods for H5SparseMatrixSeed objects, as
### well as behind coercion from CSC_H5SparseMatrixSeed to SVT_SparseMatrix.
### Does NOT propagate the dimnames.
###
### Notes:
### - SparseArray:::make_SVT_SparseMatrix_from_CSC() will fail if passed
###   long vectors via its 'data' and/or 'row_indices' arguments because R
###   does not support passing long vectors to the .Call interface yet!
###   So we use a block strategy where we load blocks of adjacent columns
###   and convert them to SVT_SparseMatrix objects, then cbind() all the
###   objects together. By default, blocks are made of 125 millions
###   data/indices elements.
### - Supports parallelization via the 'BPPARAM' argument. However some
###   quick testing with 'BiocParallel::MulticoreParam(2)' on a powerful
###   Linux server seemed to indicate that it's not worth it. Execution
###   time remained about the same but memory footprint increased
###   significantly!

.load_CSC_H5SparseMatrixSeed <- function(x, j=NULL,
                                         DATABLOCKLEN=125000000L,
                                         BPPARAM=NULL)
{
    stopifnot(is(x, "CSC_H5SparseMatrixSeed"),
              isSingleInteger(DATABLOCKLEN), DATABLOCKLEN >= 0L)
    if (is.null(j)) {
        ans_ncol <- ncol(x)
        w <- x@indptr_ranges[ , "width"]
    } else {
        stopifnot(is.integer(j))
        ans_ncol <- length(j)
        if (ans_ncol != 0L)
            stopifnot(isStrictlySorted(j),
                      1L <= j[[1L]], j[[ans_ncol]] <= ncol(x))
        w <- x@indptr_ranges[j , "width"]
    }
    ans_dim <- c(nrow(x), ans_ncol)
    ## 'cumsum(as.double(w))' instead of 'cumsum(w)' to avoid integer overflow.
    ans_indptr <- c(0, cumsum(as.double(w)))
    ans_nzcount <- ans_indptr[[length(ans_indptr)]]

    ## DATABLOCKLEN == 0L means no block processing.
    if (DATABLOCKLEN == 0L || ans_nzcount <= DATABLOCKLEN) {
        if (is.null(j)) {
            start <- count <- NULL
        } else {
            start <- x@indptr_ranges[j, "start"]
            count <- x@indptr_ranges[j, "width"]
        }
        ans_data <- .read_h5sparse_data(x@filepath, x@group, x@subdata,
                                        start=start, count=count)
        ans_row_indices <- .read_h5sparse_indices(x@filepath, x@group,
                                        start=start, count=count)
        ans <- SparseArray:::make_SVT_SparseMatrix_from_CSC(ans_dim,
                                        ans_indptr, ans_data, ans_row_indices)
        return(ans)
    }

    ## Compute 'nblock' (will always be >= 2).
    nblock <- ans_nzcount %/% DATABLOCKLEN
    if (ans_nzcount %% DATABLOCKLEN != 0L)
        nblock <- nblock + 1L

    ## Partition column indices in ranges (nb of ranges is guaranteed to be
    ## >= 1 and <= 'min(nblock, ans_ncol)').
    col_ranges <- breakInChunks(ans_ncol, nblock)
    ## There will be zero-width ranges if and only if 'nblock' > 'ans_ncol'.
    ## Drop them.
    col_ranges <- col_ranges[width(col_ranges) != 0L]
    s <- start(col_ranges)
    e <- end(col_ranges)

    ## Load ranges of columns into SVT_SparseMatrix objects.
    objects <- S4Arrays:::bplapply2(seq_along(col_ranges),
        function(b, x, j, s, e) {
            k1 <- s[[b]]
            k2 <- e[[b]]
            jj <- if (is.null(j)) k1:k2 else j[k1:k2]
            ## Set 'DATABLOCKLEN' to 0L to disable block processing.
            .load_CSC_H5SparseMatrixSeed(x, jj, DATABLOCKLEN=0L)
        },
        x, j, s, e,
        BPPARAM=BPPARAM
    )

    ## Combine all objects together.
    do.call(cbind, objects)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sparse_array() and extract_array() methods
###

.extract_sparse_array_from_CSC_H5SparseMatrixSeed <- function(x, index)
{
    j <- index[[2L]]
    if (!is.null(j)) {
        if (!is.integer(j))
            j <- as.integer(j)
        sort_j <- !isStrictlySorted(j)
        if (sort_j) {
            j0 <- j
            j <- unique(sort(j))
        }
    }
    svt <- .load_CSC_H5SparseMatrixSeed(x, j=j)
    index2 <- list(index[[1L]], NULL)
    if (!is.null(j) && sort_j)
        index2[[2L]] <- match(j0, j)
    extract_sparse_array(svt, index2)
}

setMethod("extract_sparse_array", "CSC_H5SparseMatrixSeed",
    function(x, index)
        .extract_sparse_array_from_CSC_H5SparseMatrixSeed(x, index)
)

setMethod("extract_sparse_array", "CSR_H5SparseMatrixSeed",
    function(x, index)
        t(.extract_sparse_array_from_CSC_H5SparseMatrixSeed(t(x), rev(index)))
)

setMethod("extract_array", "H5SparseMatrixSeed",
    function(x, index) as.array(extract_sparse_array(x, index))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractNonzeroDataByCol() and extractNonzeroDataByRow()
###
### TODO: Deprecate these 2 generics and their methods. These 2 generics are
### weird and don't have good/strong use cases. I suspect nobody uses them
### nor is aware of them.

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

### 'j1' and 'j2' must be 2 single integers representing a valid range of
### col indices.
### Returns a NumericList or IntegerList object parallel
### to 'j1:j2' i.e. with one list element per col index in 'j1:j2'.
.extract_data_from_adjacent_cols <- function(x, j1, j2)
{
    j12 <- j1:j2
    start <- x@indptr_ranges[j1, "start"]
    count_per_col <- x@indptr_ranges[j12, "width"]
    count <- sum(count_per_col)
    ans_nzdata <- .read_h5sparse_data(x@filepath, x@group, x@subdata,
                                      start=start, count=count)
    relist(ans_nzdata, PartitioningByWidth(count_per_col))
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

