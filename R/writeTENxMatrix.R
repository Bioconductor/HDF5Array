### =========================================================================
### writeTENxMatrix()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers to create a TENxMatrix realization sink and append data
### to it
###

.write_TENx_component <- function(filepath, group, name, data,
                                  H5type=NULL, chunk.length=NULL)
{
    name <- paste0(group, "/", name)
    ## Starting with rhdf5 2.57.10, fixed-size string datasets no longer
    ## support NAs. The recommendation now is to use a variable-length string
    ## dataset (by setting 'size' to NULL, the default) if there's a chance
    ## that the character array to write to disk contains NAs.
    #size <- h5mread:::compute_max_string_size(data)
    data_len <- length(data)
    if (is.null(chunk.length) || chunk.length > data_len) {
        chunk_len <- data_len
    } else {
        chunk_len <- chunk.length
    }
    create_and_log_HDF5_dataset(filepath, name, dim=data_len,
                                type=typeof(data), H5type=H5type, # size=size,
                                chunkdim=chunk_len, level=0L)
    h5write(data, filepath, name)
}

.write_shape <- function(filepath, group, shape)
{
    ## Standard HDF5 type H5T_STD_U32LE: unsigned 32-bit integer, little-endian
    .write_TENx_component(filepath, group, "shape", shape,
                          H5type="H5T_STD_U32LE")
}

.write_genes <- function(filepath, group, genes)
{
    .write_TENx_component(filepath, group, "genes", genes,
                          chunk.length=2048L)
}

.write_barcodes <- function(filepath, group, barcodes)
{
    .write_TENx_component(filepath, group, "barcodes", barcodes,
                          chunk.length=4096L)
}

.create_empty_data <- function(filepath, group, maxlen, type, level)
{
    name <- paste0(group, "/data")
    create_and_log_HDF5_dataset(filepath, name, dim=0L, maxdim=maxlen,
                                type=type, chunkdim=16384L, level=level)
}

.create_empty_row_indices <- function(filepath, group, maxlen, level)
{
    name <- paste0(group, "/indices")
    ## Standard HDF5 type H5T_STD_U32LE: unsigned 32-bit integer, little-endian
    create_and_log_HDF5_dataset(filepath, name, dim=0L, maxdim=maxlen,
                                type="integer", H5type="H5T_STD_U32LE",
                                chunkdim=16384L, level=level)
}

.create_empty_indptr <- function(filepath, group, ncol)
{
    name <- paste0(group, "/indptr")
    ## 10x Genomics uses HDF5 standard type H5T_STD_I64LE (64-bit signed
    ## integer, little-endian) for the "indptr" dataset, so we do the same.
    create_and_log_HDF5_dataset(filepath, name, dim=0L, maxdim=ncol+1L,
                                type="integer", H5type="H5T_STD_I64LE",
                                chunkdim=4096L, level=0L)
    h5append(filepath, name, 0)
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
    dim <- S4Arrays:::normarg_dim(dim)
    if (length(dim) != 2L)
        stop(wmsg("TENxMatrix backend only supports ",
                  "realization of matrix-like objects"))
    if (is.null(dimnames)) {
        dimnames <- vector("list", length(dim))
    } else {
        if (!(is.list(dimnames) && length(dimnames) == length(dim)))
            stop(wmsg("'dimnames' must be NULL or a list ",
                      "with 1 list element per dimension"))
    }
    if (is.null(filepath)) {
        filepath <- getHDF5DumpFile()
    } else {
        filepath <- normalize_dump_filepath(filepath)
    }
    if (is.null(group)) {
        group <- getHDF5DumpName(for.use=TRUE)
    } else {
        group <- normalize_dump_name(group)
    }
    if (is.null(level)) {
        level <- getHDF5DumpCompressionLevel()
    } else {
        level <- normalize_compression_level(level)
    }
    ok <- h5createGroup(filepath, group)
    if (!ok)
        stop(wmsg("failed to create group '", group, "' ",
                  "in file '", filepath, "'"), call.=FALSE)

    .write_shape(filepath, group, dim)
    if (!is.null(dimnames)) {
        rownames <- dimnames[[1L]]
        if (!is.null(rownames))
            .write_genes(filepath, group, rownames)
        colnames <- dimnames[[2L]]
        if (!is.null(colnames))
            .write_barcodes(filepath, group, colnames)
    }
    .create_empty_data(filepath, group, prod(dim), type, level)
    .create_empty_row_indices(filepath, group, prod(dim), level)
    .create_empty_indptr(filepath, group, dim[[2L]])
    new2("TENxRealizationSink", dim=dim, dimnames=dimnames, type=type,
                                filepath=filepath, group=group)
}

### Defining this method will force writeTENxMatrix() (thru
### BLOCK_write_to_sink() thru sinkApply()) to write blocks that
### span full columns which is a requirement of the write_block()
### method for TENxRealizationSink objects. See below.
setMethod("chunkdim", "TENxRealizationSink",
    function(x) c(nrow(x), min(ncol(x), 1L))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write_block() and close() methods for TENxRealizationSink objects
###

### The current length of 'indptr' is the "current 1-based column index"
### i.e. the nb of columns written so far + 1.
.get_current_col_index <- function(filepath, group)
{
    h5length(filepath, paste0(group, "/indptr"))
}

.check_viewport <- function(viewport, sink)
{
    if (!identical(nrow(viewport), nrow(sink))) {
        s <- S4Arrays:::make_string_from_ArrayViewport(viewport,
                                                       with.brackets=TRUE)
        stop(wmsg("Invalid viewport passed to write_block(): ", s),
             "\n\n  ",
             wmsg("The write_block() method for ", class(sink), " objects ",
                  "can only be used to write a block to a viewport that ",
                  "spans full columns, that is, to a viewport for which ",
                  "'nrow(viewport)' is equal to 'nrow(sink)'."))
    }
    current_col_idx <- .get_current_col_index(sink@filepath, sink@group)
    if (!identical(start(viewport)[[2L]], current_col_idx)) {
        s <- S4Arrays:::make_string_from_ArrayViewport(viewport,
                                                       with.brackets=TRUE)
        stop(wmsg("Invalid viewport passed to write_block(): ", s),
             "\n\n  ",
             wmsg("The block to write must start at the beginning of the ",
                  "sink or be adjacent to the previously written block."),
             "\n\n  ",
             wmsg("The write_block() method for ", class(sink), " objects ",
                  "can only be used in \"append block mode\", that is, the ",
                  "first block written to the sink must be written at the ",
                  "beginning of the sink and any subsequent block must be ",
                  "written to a viewport that is adjacent to the viewport ",
                  "where the previous block was written."))
    }
}

.append_data <- function(filepath, group, data)
{
    name <- paste0(group, "/data")
    h5append(filepath, name, data)
}

.append_row_indices <- function(filepath, group, row_indices)
{
    name <- paste0(group, "/indices")
    h5append(filepath, name, row_indices)
}

### Returns the last value in the updated "indptr" dataset.
.append_indptr <- function(filepath, group, block_indptr)
{
    name <- paste0(group, "/indptr")
    old_len <- h5length(filepath, name)
    old_data_len <- h5mread(filepath, name, starts=list(old_len),
                            as.vector=TRUE)
    indptr <- as.double(old_data_len) + block_indptr
    h5append(filepath, name, indptr)
    indptr[[length(indptr)]]
}

### Supports "append block mode" only.
setMethod("write_block", "TENxRealizationSink",
    function(sink, viewport, block)
    {
        .check_viewport(viewport, sink)
        nzdata <- nzvals(block)
        Mindex <- nzwhich(block, arr.ind=TRUE)

        ## Append the nonzero data.
        new_data_len1 <- .append_data(sink@filepath, sink@group, nzdata)

        ## Append the 0-based row indices of the nonzero data.
        new_data_len2 <- .append_row_indices(sink@filepath, sink@group,
                                             Mindex[ , 1L] - 1L)
        stopifnot(new_data_len2 == new_data_len1)  # sanity check

        ## Append the "indptr" values.
        partitioning <- PartitioningByEnd(Mindex[ , 2L], NG=ncol(block))
        block_indptr <- end(partitioning)  # sorted integer vector
        new_data_len3 <- .append_indptr(sink@filepath, sink@group, block_indptr)
        stopifnot(new_data_len3 == new_data_len1)  # sanity check
        sink
    }
)

### Only performs some sanity checks (there is actually nothing to close).
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
### Coercing a TENxRealizationSink object
###

setAs("TENxRealizationSink", "TENxMatrixSeed",
    function(from) TENxMatrixSeed(from@filepath, from@group)
)

setAs("TENxRealizationSink", "TENxMatrix",
    function(from) DelayedArray(as(from, "TENxMatrixSeed"))
)

setAs("TENxRealizationSink", "DelayedArray",
    function(from) DelayedArray(as(from, "TENxMatrixSeed"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeTENxMatrix()
###

### Write the dataset to the current dump if 'filepath' and 'group' are not
### specified.
### Returns a TENxMatrix object pointing to the newly written HDF5-based
### sparse matrix on disk.
writeTENxMatrix <- function(x, filepath=NULL, group=NULL,
                               level=NULL, verbose=NA)
{
    verbose <- DelayedArray:::normarg_verbose(verbose)
    sink <- TENxRealizationSink(dim(x), dimnames(x), type(x),
                                filepath=filepath, group=group, level=level)
    sink <- BLOCK_write_to_sink(sink, x, verbose=verbose)
    ans <- as(sink, "TENxMatrix")
    if (verbose)
        message("sparsity: ", round(sparsity(ans), digits=2))
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to TENxMatrix
###
### The methods below write the object to disk. Note that coercion from
### TENxRealizationSink to TENxMatrix is already taken care of by the specific
### method above and doesn't write anything to disk. So coercing to TENxMatrix
### in general writes the object to disk *except* when the object to coerce is
### a TENxRealizationSink object.
###

### Write to current dump.
.as_TENxMatrix <- function(from) writeTENxMatrix(from)

setAs("ANY", "TENxMatrix", .as_TENxMatrix)

### Automatic coercion method from DelayedArray to TENxMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "TENxMatrix", .as_TENxMatrix)
setAs("DelayedMatrix", "TENxMatrix", .as_TENxMatrix)

