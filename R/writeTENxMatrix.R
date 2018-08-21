### =========================================================================
### writeTENxMatrix()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers to create a TENxMatrix realization sink and append data
### to it
###

.create_data <- function(filepath, group, dim, type, level=0L)
{
    name <- paste0(group, "/data")
    h5createDataset2(filepath, name,
                     dim=0, maxdim=prod(dim),
                     type=type, chunkdim=16384L, level=level)
}

.create_row_indices <- function(filepath, group, dim, level=0L)
{
    name <- paste0(group, "/indices")
    h5createDataset2(filepath, name,
                     dim=0, maxdim=prod(dim),
                     type="integer", chunkdim=8192L, level=level)
}

.create_indptr <- function(filepath, group, dim, level=0L)
{
    name <- paste0(group, "/indptr")
    ## The values in the "indptr" dataset will be sorted and its last value
    ## (which is also its biggest) should always be the length of the "data"
    ## and "indices" datasets, so it can be >= 2^31. Because we don't know
    ## how to write values that are >= 2^31 to a dataset of type INTEGER
    ## (see https://github.com/grimbough/rhdf5/issues/21), we make the
    ## dataset of type FLOAT.
    h5createDataset2(filepath, name,
                     dim=0, maxdim=dim[[2L]] + 1L,
                     type="double", chunkdim=8192L, level=level)
    h5append(0, filepath, name)
}

### The current length of 'indptr' is the "current 1-based column index"
### i.e. the nb of columns written so far + 1.
.get_current_col_index <- function(filepath, group)
{
    h5length(filepath, paste0(group, "/indptr"))
}

.append_data <- function(filepath, group, data)
{
    name <- paste0(group, "/data")
    h5append(data, filepath, name)
}

.append_row_indices <- function(filepath, group, row_indices)
{
    name <- paste0(group, "/indices")
    h5append(row_indices, filepath, name)
}

### Return the last value in the extended "indptr" dataset.
.append_indptr <- function(filepath, group, col_indices, ncol)
{
    name <- paste0(group, "/indptr")
    old_len <- h5length(filepath, name)
    old_data_len <- as.vector(h5read(filepath, name, start=old_len, count=1L))
    indptr <- end(PartitioningByEnd(col_indices, NG=ncol)) + old_data_len
    new_len <- h5append(indptr, filepath, name)
    as.vector(h5read(filepath, name, start=new_len, count=1L))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TENxRealizationSink objects
###
### The TENxRealizationSink class is a concrete RealizationSink subclass that
### implements a TENxMatrix realization sink.
###

setClass("TENxRealizationSink",
    contains="RealizationSink",
    representation(
        dim="integer",          # Naming this slot "dim" makes dim() work
                                # out of the box.
        dimnames="list",
        type="character",       # Single string.
        filepath="character",   # Single string.
        group="character"       # Name of the group in the HDF5 file
                                # where to write the data.
    )
)

setMethod("dimnames", "TENxRealizationSink",
    function(x)
    {
        ans <- x@dimnames
        if (all(S4Vectors:::sapply_isNULL(ans)))
            return(NULL)
        ans
    }
)

setMethod("type", "TENxRealizationSink", function(x) x@type)

TENxRealizationSink <- function(dim, dimnames=NULL, type="double",
                                filepath=NULL, group=NULL, level=NULL)
{
    if (!(is.integer(dim) && length(dim) == 2L))
        stop(wmsg("'dim' must be an integer vector of length 2"))
    if (S4Vectors:::anyMissingOrOutside(dim, 0L))
        stop(wmsg("'dim' cannot contain NAs or negative values"))
    if (is.null(dimnames)) {
        dimnames <- vector("list", length(dim))
    } else {
        if (!(is.list(dimnames) && length(dimnames) == length(dim)))
            stop(wmsg("'dimnames' must be NULL or a list ",
                      "with 1 list element per dimension"))
    }
    if (is.null(filepath)) {
        filepath <- getHDF5DumpFile(for.use=TRUE)
    } else {
        filepath <- normalize_dump_filepath(filepath)
    }
    if (is.null(group)) {
        group <- getHDF5DumpName(for.use=TRUE)
    } else {
        group <- normalize_dump_name(group)
    }
    ## Let's not compress for now.
    if (is.null(level)) {
        #level <- getHDF5DumpCompressionLevel()
        level <- 0L
    } else {
        #level <- normalize_compression_level(level)
        stop(wmsg("'level' not supported yet"))
    }
    ok <- h5createGroup(filepath, group)
    if (!ok)
        stop(wmsg("failed to create group '", group, "' ",
                  "in file '", filepath, "'"), call.=FALSE)

    h5write(dim, filepath, paste0(group, "/shape"))
    .create_data(filepath, group, dim, type, level=level)
    .create_row_indices(filepath, group, dim, level=level)
    .create_indptr(filepath, group, dim, level=level)

    #appendDatasetCreationToHDF5DumpLog(filepath, name, dim, type,
    #                                   chunkdim, level)
    new2("TENxRealizationSink", dim=dim, dimnames=dimnames, type=type,
                                filepath=filepath, group=group)
}

### Support "append mode" only.
setMethod("write_block", "TENxRealizationSink",
    function(x, viewport, block)
    {
        if (!identical(nrow(block), nrow(x)))
            stop(wmsg("The \"write_block\" method for ", class(x), " objects ",
                      "can only be used to write blocks of full columns ",
                      "i.e. blocks such that 'nrow(block) == nrow(x)'."))

        current_col_idx <- .get_current_col_index(x@filepath, x@group)
        if (!identical(start(viewport)[[2L]], current_col_idx))
            stop(wmsg("The block to write is not adjacent to the last ",
                      "written block. ",
                      "The \"write_block\" method for ", class(x), " objects ",
                      "can only be used to append blocks i.e. the new block ",
                      "to write must be adjacent to the last written block ",
                      "(except for the 1st written block which must be ",
                      "written at the beginning of the sink)."))

        non_zero_data_aind <- which(block != 0L, arr.ind=TRUE)

        ## Append the non-zero data.
        new_data_len1 <- .append_data(x@filepath, x@group,
                                      block[non_zero_data_aind])

        ## Append the 0-based row indices of the non-zero data.
        new_data_len2 <- .append_row_indices(x@filepath, x@group,
                                      non_zero_data_aind[ , "row"] - 1L)
        stopifnot(new_data_len2 == new_data_len1)  # sanity check

        ## Append the "indptr" values.
        new_data_len3 <- .append_indptr(x@filepath, x@group,
                                      non_zero_data_aind[ , "col"],
                                      ncol(block))
        stopifnot(new_data_len3 == new_data_len1)  # sanity check
    }
)

setMethod("close", "TENxRealizationSink",
    function(con)
    {
        current_col_idx <- .get_current_col_index(con@filepath, con@group)
        if (current_col_idx <= ncol(con))
            stop(wmsg("cannot close ", class(con), " object before ",
                      "writing all data to it"))
        stopifnot(current_col_idx == ncol(con) + 1L)  # should never happen
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercing an TENxRealizationSink object.
###

setAs("TENxRealizationSink", "TENxMatrixSeed",
    function(from) TENxMatrixSeed(from@filepath, from@group)
)

setAs("TENxRealizationSink", "TENxMatrix",
    function(from) TENxMatrix(as(from, "TENxMatrixSeed"))
)

setAs("TENxRealizationSink", "DelayedArray",
    function(from) TENxMatrix(as(from, "TENxMatrixSeed"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeTENxMatrix()
###

### Write the dataset to the current dump if 'filepath' and 'group' are not
### specified.
### Return an TENxMatrix object pointing to the newly written HDF5-based
### sparse matrix on disk.
writeTENxMatrix <- function(x, filepath=NULL, group=NULL,
                            level=NULL, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    sink <- TENxRealizationSink(dim(x), dimnames(x), type(x),
                                filepath=filepath, group=group, level=level)
    if (verbose) {
        old_verbose <- DelayedArray:::set_verbose_block_processing(verbose)
        on.exit(DelayedArray:::set_verbose_block_processing(old_verbose))
    }
    write_array_to_sink(x, sink)
    as(sink, "TENxMatrix")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to TENxMatrix.
###
### The methods below write the object to disk. Note that coercion from
### TENxRealizationSink to TENxMatrix is already taken care of by the specific
### method above and doesn't write anything to disk. So coercing to TENxMatrix
### in general writes the object to disk *except* when the object to coerce is
### an TENxRealizationSink object.
###

.as_TENxMatrix <- function(from) writeTENxMatrix(from)  # write to current dump

setAs("ANY", "TENxMatrix", .as_TENxMatrix)

### Automatic coercion method from DelayedArray to TENxMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "TENxMatrix", .as_TENxMatrix)

