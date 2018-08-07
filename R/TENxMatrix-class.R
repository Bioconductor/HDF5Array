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

.get_barcodes <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "barcodes", idx=idx)

.get_data <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "data", idx=idx)

### Currently unused.
.get_gene_names <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "gene_names", idx=idx)

.get_genes <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "genes", idx=idx)

### Return 0-based row indices.
.get_row_indices <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "indices", idx=idx)

.get_indptr <- function(filepath, group, idx=NULL)
{
    name <- paste0(group, "/indptr")
    if (!is.null(idx))
        idx <- list(idx)
    as.vector(h5read(filepath, name, index=idx, bit64conversion="double"))
}

.get_shape <- function(filepath, group, idx=NULL)
    .get_TENx_component(filepath, group, "shape", idx=idx)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .get_data_indices_by_col()
###

### S4Vectors:::fancy_mseq() does not accept 'offset' of type double yet so
### we implement a version that does.
.fancy_mseq <- function(lengths, offset=0)
{
    lengths_len <- length(lengths)
    if (lengths_len == 0L)
        return(numeric(0))
    offsets <- offset - cumsum(c(0L, lengths[-lengths_len]))
    seq_len(sum(lengths)) + rep.int(offsets, lengths)
}

### 'j' must be NULL or an integer vector containing valid col indices.
### Return data indices in a NumericList object parallel to 'j' i.e. with
### one list element per col index in 'j'.
.get_data_indices_by_col <- function(x, j)
{
    col_ranges <- x@col_ranges
    if (!is.null(j))
        col_ranges <- S4Vectors:::extract_data_frame_rows(col_ranges, j)
    start2 <- col_ranges[ , "start"]
    width2 <- col_ranges[ , "width"]
    idx2 <- .fancy_mseq(width2, offset=start2 - 1L)
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

### 'i' (or 'j') must be NULL or an integer vector containing valid
### row (or col) indices.
### Return the data in a data.frame with 3 columns: i, j, data.
.extract_sparse_data_from_TENxMatrixSeed <- function(x, i, j)
{
    data_indices <- .get_data_indices_by_col(x, j)
    idx2 <- unlist(data_indices, use.names=FALSE)
    i2 <- .get_row_indices(x@filepath, x@group, idx=idx2) + 1L
    j2 <- rep.int(seq_along(data_indices), lengths(data_indices))
    if (!is.null(i)) {
        m <- findMatches(i2, i)
        idx2 <- idx2[from(m)]
        i2 <- to(m)
        j2 <- j2[from(m)]
    }
    data <- .get_data(x@filepath, x@group, idx=idx2)
    data.frame(i=i2, j=j2, data=data, stringsAsFactors=FALSE)
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
    i <- index[[1L]]
    j <- index[[2L]]
    sparse_data <- .extract_sparse_data_from_TENxMatrixSeed(x, i, j)
    .make_matrix_from_sparse_data(ans_dim, sparse_data)
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

