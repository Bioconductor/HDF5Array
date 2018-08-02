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
                                 # containing the 10xGenomics data.
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
                                            "10xGenomics data")
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

.get_indices <- function(filepath, group, idx=NULL)
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
### extract_array()
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

.extract_array_from_TENxMatrixSeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    ans <- array(0L, dim=ans_dim)

    j <- index[[2L]]
    col_ranges <- x@col_ranges
    if (is.null(j)) {
        j <- seq_len(ncol(ans))
    } else {
        col_ranges <- S4Vectors:::extract_data_frame_rows(col_ranges, j)
    }
    start2 <- col_ranges[ , "start"]
    width2 <- col_ranges[ , "width"]
    idx2 <- .fancy_mseq(width2, offset=start2 - 1)
    i2 <- .get_indices(x@filepath, x@group, idx=idx2) + 1L
    j2 <- rep.int(seq_along(j), width2)

    i <- index[[1L]]
    if (!is.null(i)) {
        m <- findMatches(i2, i)
        idx2 <- idx2[from(m)]
        i2 <- to(m)
        j2 <- j2[from(m)]
    }

    ans[cbind(i2, j2)] <- .get_data(x@filepath, x@group, idx=idx2)
    ans
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
    filepath <- normarg_path(filepath, "'filepath'", "10xGenomics data")
    if (!isSingleString(group))
        stop(wmsg("'group' must be a single string specifying the name ",
                  "of the group in the HDF5 file containing the ",
                  "10xGenomics data"))
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

