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
        dim="integer",              # Naming this slot "dim" makes dim() work
                                    # out of the box.
        dimnames="list",
        type="character",           # Single string.
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

### Unlike with rhdf5::h5createDataset(), if 'chunkdim' is NULL then an
### automatic chunk geometry will be used. To write "unchunked data" (a.k.a.
### contiguous data), 'chunkdim' must be set to 0.
### FIXME: Investigate the possiblity to write the dimnames to the HDF5 file.
HDF5RealizationSink <- function(dim, dimnames=NULL, type="double",
                                filepath=NULL, name=NULL,
                                chunkdim=NULL, level=NULL)
{
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
                                type=type, chunkdim=chunkdim, level=level)
    if (is.null(dimnames)) {
        dimnames <- vector("list", length(dim))
    } else {
        ## TODO: Write the dimnames to the HDF5 file.
    }
    new2("HDF5RealizationSink", dim=dim, dimnames=dimnames, type=type,
                                filepath=filepath, name=name,
                                chunkdim=chunkdim)
}

setMethod("chunkdim", "HDF5RealizationSink", function(x) x@chunkdim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Writing data to an HDF5RealizationSink object
###

setMethod("write_block", "HDF5RealizationSink",
    function(x, viewport, block)
    {
        h5write(block, x@filepath, x@name,
                start=start(viewport), count=width(viewport))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercing an HDF5RealizationSink object
###

### FIXME: This coercion needs to propagate the dimnames *thru* the HDF5 file.
### For more details about this, see FIXME right before definition of
### HDF5RealizationSink() above in this file and right before definition of
### HDF5ArraySeed() in HDF5Array-class.R.
setAs("HDF5RealizationSink", "HDF5ArraySeed",
    function(from) HDF5ArraySeed(from@filepath, from@name, type=from@type)
)

### Note that this coercion currently drops the dimnames but will naturally
### propagate them when coercion from HDF5RealizationSink to HDF5ArraySeed
### propagates them. See FIXME above.
setAs("HDF5RealizationSink", "HDF5Array",
    function(from) DelayedArray(as(from, "HDF5ArraySeed"))
)

setAs("HDF5RealizationSink", "DelayedArray",
    function(from)
    {
        ans <- HDF5Array(as(from, "HDF5ArraySeed"))
        ## Temporarily needed because coercion from HDF5RealizationSink to
        ## HDF5ArraySeed does not propagate the dimnames at the moment. See
        ## FIXME above.
        ## TODO: Remove line below when FIXME above is addressed.
        dimnames(ans) <- dimnames(from)
        ans
    }
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
### FIXME: This needs to write the dimnames to the file. See various FIXMEs
### above in this file about this.
writeHDF5Array <- function(x, filepath=NULL, name=NULL, chunkdim=NULL,
                           level=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    sink <- HDF5RealizationSink(dim(x), dimnames(x), type(x),
                                filepath=filepath, name=name,
                                chunkdim=chunkdim, level=level)
    if (verbose) {
        old_verbose <- DelayedArray:::set_verbose_block_processing(verbose)
        on.exit(DelayedArray:::set_verbose_block_processing(old_verbose))
    }
    BLOCK_write_to_sink(x, sink)
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

.as_HDF5Array <- function(from) writeHDF5Array(from)  # write to current dump

setAs("ANY", "HDF5Array", .as_HDF5Array)

### Automatic coercion methods from DelayedArray to HDF5Array and from
### DelayedMatrix to HDF5Matrix silently return broken objects (unfortunately
### these dummy automatic coercion methods don't bother to validate the object
### they return). So we overwrite them.
setAs("DelayedArray", "HDF5Array", .as_HDF5Array)
setAs("DelayedMatrix", "HDF5Matrix", .as_HDF5Array)

