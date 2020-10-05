### =========================================================================
### writeHDF5Array()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5RealizationSink objects
###
### The HDF5RealizationSink class is a concrete RealizationSink subclass that
### implements an HDF5Array realization sink.
###

setClass("HDF5RealizationSink",
    contains="RealizationSink",
    representation(
        ## Slots that support the RealizationSink constructor contract.
        dim="integer",              # Naming this slot "dim" makes dim() work
                                    # out of the box.
        dimnames="list",
        type="character",           # Single string.
        as_sparse="logical",        # TRUE or FALSE.

        ## Other slots.
        filepath="character",       # Single string.
        name="character",           # Dataset name.
        chunkdim="integer_OR_NULL"  # An integer vector parallel to the 'dim'
                                    # slot or NULL.
    )
)

setMethod("dimnames", "HDF5RealizationSink",
    function(x)
    {
        ans <- x@dimnames
        if (all(S4Vectors:::sapply_isNULL(ans)))
            return(NULL)
        ans
    }
)

setMethod("type", "HDF5RealizationSink", function(x) x@type)

setMethod("chunkdim", "HDF5RealizationSink", function(x) x@chunkdim)

setMethod("is_sparse", "HDF5RealizationSink", function(x) x@as_sparse)

.normarg_chunkdim <- function(chunkdim, dim)
{
    if (!(is.numeric(chunkdim) || is.logical(chunkdim) && all(is.na(chunkdim))))
        stop(wmsg("'chunkdim' must be NULL or an integer vector"))
    if (!is.integer(chunkdim))
        chunkdim <- as.integer(chunkdim)
    if (length(chunkdim) != length(dim))
        stop(wmsg("'chunkdim' must be an integer vector of length ",
                  "the number of dimensions of the object to write"))
    if (!all(chunkdim <= dim, na.rm=TRUE))
        stop(wmsg("the chunk dimensions specified in 'chunkdim' exceed ",
                  "the dimensions of the object to write"))
    if (any(chunkdim == 0L & dim != 0L, na.rm=TRUE))
        stop(wmsg("'chunkdim' must contain nonzero values unless ",
                  "the zero values correspond to dimensions in the ",
                  "object to write that are also zero"))
    na_idx <- which(is.na(chunkdim))
    chunkdim[na_idx] <- dim[na_idx]
    if (prod(chunkdim) > .Machine$integer.max)
        stop(wmsg("The chunk dimensions in 'chunkdim' are too big. The ",
                  "product of the chunk dimensions should always be <= ",
                  ".Machine$integer.max"))
    chunkdim
}

### Note that the supplied 'as.sparse' value is stored in the 'as_sparse'
### slot of the returned object, and that's all. It doesn't change how the
### data will be laid out to the HDF5 file in anyway (HDF5 doesn't support
### sparse storage at the moment). The only reason we store the supplied
### 'as.sparse' value in the object is so that we can propagate it later
### when we coerce the object to HDF5ArraySeed.
### Unlike with rhdf5::h5createDataset(), if 'chunkdim' is NULL then an
### automatic chunk geometry will be used. To write "unchunked data" (a.k.a.
### contiguous data), 'chunkdim' must be set to 0.
HDF5RealizationSink <- function(dim, dimnames=NULL, type="double",
                                as.sparse=FALSE,
                                filepath=NULL, name=NULL,
                                H5type=NULL, size=NULL,
                                chunkdim=NULL, level=NULL)
{
    if (!isTRUEorFALSE(as.sparse))
        stop(wmsg("'as.sparse' must be TRUE or FALSE"))
    if (is.null(filepath)) {
        filepath <- getHDF5DumpFile(for.use=TRUE)
    } else {
        filepath <- normalize_dump_filepath(filepath)
    }
    if (is.null(name)) {
        name <- getHDF5DumpName(for.use=TRUE)
    } else {
        name <- normalize_dump_name(name)
    }
    if (is.null(chunkdim)) {
        ## TODO: Pass 'x' instead of 'dim' to getHDF5DumpChunkDim() and modify
        ## getHDF5DumpChunkDim() to return 'chunkdim(x)' if it's not NULL.
        ## See TODO comment in dump-management.R
        chunkdim <- getHDF5DumpChunkDim(dim)
    } else if (isSingleNumber(chunkdim) && chunkdim == 0) {
        chunkdim <- NULL  # no chunking
    } else {
        chunkdim <- .normarg_chunkdim(chunkdim, dim)
    }
    if (is.null(level)) {
        if (is.null(chunkdim)) {
            level <- 0L
        } else {
            level <- getHDF5DumpCompressionLevel()
        }
    } else {
        level <- normalize_compression_level(level)
    }
    create_and_log_HDF5_dataset(filepath, name, dim,
                                type=type, H5type=H5type, size=size,
                                chunkdim=chunkdim, level=level)
    if (is.null(dimnames)) {
        dimnames <- vector("list", length(dim))
    } else {
        h5writeDimnames(dimnames, filepath, name)
    }
    new2("HDF5RealizationSink", dim=dim, dimnames=dimnames, type=type,
                                as_sparse=as.sparse,
                                filepath=filepath, name=name,
                                chunkdim=chunkdim)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Writing data to an HDF5RealizationSink object
###

setMethod("write_block", "HDF5RealizationSink",
    function(sink, viewport, block)
    {
        if (!is.array(block))
            block <- as.array(block)
        h5write(block, sink@filepath, sink@name,
                start=start(viewport), count=width(viewport))
        sink
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercing an HDF5RealizationSink object
###

setAs("HDF5RealizationSink", "HDF5ArraySeed",
    function(from) HDF5ArraySeed(from@filepath, from@name,
                                 as.sparse=from@as_sparse)
)

setAs("HDF5RealizationSink", "HDF5Array",
    function(from) DelayedArray(as(from, "HDF5ArraySeed"))
)

setAs("HDF5RealizationSink", "DelayedArray",
    function(from) DelayedArray(as(from, "HDF5ArraySeed"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeHDF5Array()
###

### If 'filepath' and 'name' are NULL (the default), write the dataset to
### the current dump.
### If 'chunkdim' is NULL, an automatic chunk geometry will be used.
### To write "unchunked data" (a.k.a. contiguous data), 'chunkdim' must be
### set to 0.
### Return an HDF5Array object pointing to the newly written HDF5 dataset
### on disk.
writeHDF5Array <- function(x, filepath=NULL, name=NULL,
                              H5type=NULL, chunkdim=NULL, level=NULL,
                              as.sparse=NA,
                              with.dimnames=FALSE, verbose=NA)
{
    if (!(is.logical(as.sparse) && length(as.sparse) == 1L))
        stop(wmsg("'as.sparse' must be NA, TRUE or FALSE"))
    if (!isTRUEorFALSE(with.dimnames))
        stop("'with.dimnames' must be TRUE or FALSE")
    verbose <- DelayedArray:::normarg_verbose(verbose)

    if (is.na(as.sparse))
        as.sparse <- is_sparse(x)
    sink_dimnames <- if (with.dimnames) dimnames(x) else NULL
    ## compute_max_string_size() will trigger block processing if 'x' is a
    ## DelayedArray object of type "character", so it could take a while.
    size <- compute_max_string_size(x)
    sink <- HDF5RealizationSink(dim(x), sink_dimnames, type(x), as.sparse,
                                filepath=filepath, name=name,
                                H5type=H5type, size=size,
                                chunkdim=chunkdim, level=level)
    sink <- BLOCK_write_to_sink(sink, x, verbose=verbose)
    as(sink, "HDF5Array")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to HDF5Array
###
### The methods below write the object to disk. Note that coercion from
### HDF5RealizationSink to HDF5Array is already taken care of by the specific
### method above and doesn't write anything to disk. So coercing to HDF5Array
### in general writes the object to disk *except* when the object to coerce is
### an HDF5RealizationSink object.
###

### Write to current dump.
.as_HDF5Array <- function(from) writeHDF5Array(from, with.dimnames=TRUE)

setAs("ANY", "HDF5Array", .as_HDF5Array)

### Automatic coercion methods from DelayedArray to HDF5Array and from
### DelayedMatrix to HDF5Matrix silently return broken objects (unfortunately
### these dummy automatic coercion methods don't bother to validate the object
### they return). So we overwrite them.
setAs("DelayedArray", "HDF5Array", .as_HDF5Array)
setAs("DelayedMatrix", "HDF5Matrix", .as_HDF5Array)

