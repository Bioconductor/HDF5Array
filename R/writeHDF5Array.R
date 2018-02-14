### =========================================================================
### writeHDF5Array()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5RealizationSink objects
###
### The HDF5RealizationSink class is a concrete RealizationSink subclass that
### implements an HDF5 realization sink.
###

setClass("HDF5RealizationSink",
    contains="RealizationSink",
    representation(
        dim="integer",          # Naming this slot "dim" makes dim() work
                                # out of the box.
        dimnames="list",
        type="character",       # Single string.
        filepath="character",   # Single string.
        name="character",       # Dataset name.
        chunkdim="integer"      # Parallel to 'dim' slot.
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
        chunkdim <- getHDF5DumpChunkDim(dim, type)
    } else {
        chunkdim <- as.integer(chunkdim)
    }
    if (is.null(level)) {
        level <- getHDF5DumpCompressionLevel()
    } else {
        level <- normalize_compression_level(level)
    }
    h5createDataset2(filepath, name, dim, type, chunkdim, level)
    appendDatasetCreationToHDF5DumpLog(filepath, name, dim, type,
                                       chunkdim, level)
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

setMethod("write_block_to_sink", "HDF5RealizationSink",
    function(block, sink, viewport)
    {
        stopifnot(identical(dim(sink), refdim(viewport)),
                  identical(dim(block), dim(viewport)))
        index <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        h5write2(block, sink@filepath, sink@name, index=index)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercing a HDF5RealizationSink object.
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
    function(from) HDF5Array(as(from, "HDF5ArraySeed"))
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

### Write the dataset to the current dump if 'filepath' and 'name' are not
### specified.
### Return a HDF5Array object pointing to the newly written HDF5 dataset on
### disk.
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
    write_array_to_sink(x, sink)
    as(sink, "HDF5Array")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to HDF5Array.
###
### The methods below write the object to disk. Note that coercion from
### HDF5RealizationSink to HDF5Array is already taken care of by the specific
### method above and doesn't write anything to disk. So coercing to HDF5Array
### in general writes the object to disk *except* when the object to coerce is
### a HDF5RealizationSink object.
###

.as_HDF5Array <- function(from) writeHDF5Array(from)  # write to current dump

setAs("ANY", "HDF5Array", .as_HDF5Array)

### Automatic coercion method from DelayedArray to HDF5Array silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "HDF5Array", .as_HDF5Array)
setAs("DelayedMatrix", "HDF5Matrix", .as_HDF5Array)

