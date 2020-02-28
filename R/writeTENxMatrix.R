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
    if (is.character(data)) {
        size <- max(nchar(data)) + 1L
    } else {
        size <- NULL
    }
    data_len <- length(data)
    if (is.null(chunk.length) || chunk.length > data_len) {
        chunk_len <- data_len
    } else {
        chunk_len <- chunk.length
    }
    create_and_log_HDF5_dataset(filepath, name, dim=data_len,
                                type=typeof(data), H5type=H5type, size=size,
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
    ## Standard HDF5 type H5T_STD_U32LE: unsigned 32-bit integer, little-endian
    create_and_log_HDF5_dataset(filepath, name, dim=0L, maxdim=ncol+1L,
                                type="integer", H5type="H5T_STD_U32LE",
                                chunkdim=4096L, level=0L)
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
    old_data_len <- as.vector(h5mread(filepath, name, starts=list(old_len)))
    indptr <- end(PartitioningByEnd(col_indices, NG=ncol)) + old_data_len
    new_len <- h5append(indptr, filepath, name)
    as.vector(h5mread(filepath, name, starts=list(new_len)))
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
    if (!is.integer(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (length(dim) != 2L)
        stop(wmsg("TENxMatrix backend only supports ",
                  "realization of matrix-like objects"))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Writing data to a TENxRealizationSink object
###

.check_viewport <- function(viewport, x)
{
    if (!identical(nrow(viewport), nrow(x)))
        stop(wmsg("The \"write_block\" and \"write_sparse_block\" methods ",
                  "for ", class(x), " objects can only be used to write ",
                  "a block to a viewport that spans full columns i.e. to ",
                  "a viewport such that 'nrow(viewport) == nrow(x)'."))

    current_col_idx <- .get_current_col_index(x@filepath, x@group)
    if (!identical(start(viewport)[[2L]], current_col_idx))
        stop(wmsg("The block to write is not adjacent to the last ",
                  "written block.\n\n",
                  "The \"write_block\" and \"write_sparse_block\" methods ",
                  "for ", class(x), " objects can only be used ",
                  "in \"appending mode\", that is, each block must be ",
                  "written to a viewport that is adjacent to the viewport ",
                  "where the previous block was written (with the exception ",
                  "of the 1st written block which must be written to a ",
                  "viewport that starts at the beginning of the sink)."))
}

### Support "appending mode" only.
setMethod("write_sparse_block", "TENxRealizationSink",
    function(x, viewport, sparse_block)
    {
        .check_viewport(viewport, x)

        ## Append the nonzero data.
        new_data_len1 <- .append_data(x@filepath, x@group,
                                      sparse_block@nzdata)

        ## Append the 0-based row indices of the nonzero data.
        new_data_len2 <- .append_row_indices(x@filepath, x@group,
                                             sparse_block@aind[ , 1L] - 1L)
        stopifnot(new_data_len2 == new_data_len1)  # sanity check

        ## Append the "indptr" values.
        new_data_len3 <- .append_indptr(x@filepath, x@group,
                                        sparse_block@aind[ , 2L],
                                        ncol(viewport))
        stopifnot(new_data_len3 == new_data_len1)  # sanity check
    }
)

setMethod("write_block", "TENxRealizationSink",
    function(x, viewport, block)
    {
        sparse_block <- dense2sparse(block)
        write_sparse_block(x, viewport, sparse_block)
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
### Return a TENxMatrix object pointing to the newly written HDF5-based
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
    BLOCK_write_to_sink(x, sink)
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

.as_TENxMatrix <- function(from) writeTENxMatrix(from)  # write to current dump

setAs("ANY", "TENxMatrix", .as_TENxMatrix)

### Automatic coercion method from DelayedArray to TENxMatrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("DelayedArray", "TENxMatrix", .as_TENxMatrix)
setAs("DelayedMatrix", "TENxMatrix", .as_TENxMatrix)

