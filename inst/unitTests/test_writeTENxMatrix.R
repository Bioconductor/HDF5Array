
test_TENxRealizationSink <- function()
{
    checkTrue(extends("TENxRealizationSink", "RealizationSink"))

    ## --- Some edge cases (0-col or 0-row sink) ---

    for (type in c("integer", "logical", "double", "raw", "character")) {

        ## 0-col sink:

        sink <- TENxRealizationSink(c(20, 0), type=type)
        checkTrue(is(sink, "TENxRealizationSink"))
        checkIdentical(dim(sink), c(20L, 0L))
        checkIdentical(chunkdim(sink), c(20L, 0L))
        checkIdentical(type(sink), type)

        seed <- as(sink, "TENxMatrixSeed")
        checkTrue(is(seed, "TENxMatrixSeed"))
        checkIdentical(dim(seed), dim(sink))
        checkIdentical(chunkdim(seed), chunkdim(sink))
        checkIdentical(type(seed), type(sink))

        ## 0-row sink:

        sink <- TENxRealizationSink(c(0, 15), type=type)
        checkTrue(is(sink, "TENxRealizationSink"))
        checkIdentical(dim(sink), c(0L, 15L))
        checkIdentical(chunkdim(sink), c(0L, 1L))
        checkIdentical(type(sink), type)

        checkException(as(sink, "TENxMatrixSeed"))  # too soon!

        viewport <- ArrayViewport(dim(sink), IRanges(start=c(1, 1),
                                                     end  =c(nrow(sink), 11)))
        block <- matrix(vector(type), nrow=nrow(viewport), ncol=ncol(viewport))
        sink <- write_block(sink, viewport, block)

        checkException(as(sink, "TENxMatrixSeed"))  # still too soon!

        ## Invalid viewport: not adjacent to previously written block!
        viewport <- ArrayViewport(dim(sink), IRanges(start=c(1, 13),
                                                     end  =c(nrow(sink), 15)))
        block <- matrix(vector(type), nrow=nrow(viewport), ncol=ncol(viewport))
        checkException(write_block(sink, viewport, block))

        viewport <- ArrayViewport(dim(sink), IRanges(start=c(1, 12),
                                                     end  =c(nrow(sink), 15)))
        block <- matrix(vector(type), nrow=nrow(viewport), ncol=ncol(viewport))
        sink <- write_block(sink, viewport, block)

        seed <- as(sink, "TENxMatrixSeed")
        checkTrue(is(seed, "TENxMatrixSeed"))
        checkIdentical(dim(seed), dim(sink))
        checkIdentical(chunkdim(seed), chunkdim(sink))
        checkIdentical(type(seed), type(sink))
    }

    ## --- Create a TENxRealizationSink and fill it with zeros ---

    sink <- TENxRealizationSink(c(20, 15))
    checkTrue(is(sink, "TENxRealizationSink"))
    checkIdentical(dim(sink), c(20L, 15L))
    checkIdentical(chunkdim(sink), c(20L, 1L))
    checkIdentical(type(sink), "double")

    ## Invalid viewport: does not span full columns!
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 1), c(19, 10)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    checkException(write_block(sink, viewport, block))

    ## Invalid viewport: does not span full columns!
    viewport <- ArrayViewport(dim(sink), IRanges(c(2, 1), c(20, 10)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    checkException(write_block(sink, viewport, block))

    ## Invalid viewport: does not start at beginning of sink!
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 2), c(20, 10)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    checkException(write_block(sink, viewport, block))

    ## Valid viewport.
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 1), c(20, 10)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    sink <- write_block(sink, viewport, block)

    ## Invalid viewport: not adjacent to previously used viewport!
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 1), c(20, 5)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    checkException(write_block(sink, viewport, block))

    ## Invalid viewport: not adjacent to previously used viewport!
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 10), c(20, 14)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    checkException(write_block(sink, viewport, block))

    ## Invalid viewport: not adjacent to previously used viewport!
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 12), c(20, 14)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    checkException(write_block(sink, viewport, block))

    ## Valid viewport.
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 11), c(20, 14)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    sink <- write_block(sink, viewport, block)

    ## Trying to close too soon.
    checkException(as(sink, "TENxMatrixSeed"))

    ## Valid viewport.
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 15), c(20, 15)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    sink <- write_block(sink, viewport, block)

    seed <- as(sink, "TENxMatrixSeed")
    checkTrue(is(seed, "TENxMatrixSeed"))
    checkIdentical(dim(seed), dim(sink))
    checkIdentical(chunkdim(seed), chunkdim(sink))
    checkIdentical(type(seed), type(sink))

    data <- HDF5Array:::.read_h5sparse_data(path(seed), group(seed),
                                            NULL, type(seed))
    checkIdentical(data, numeric(0))

    indices <- HDF5Array:::.read_h5sparse_indices(path(seed), group(seed))
    checkIdentical(indices, integer(0))

    indptr <- HDF5Array:::.read_h5sparse_indptr(path(seed), group(seed))
    checkIdentical(indptr, numeric(1 + ncol(seed)))

    ## --- Create a TENxRealizationSink and fill it with a sparse matrix ---

    sink <- TENxRealizationSink(c(9, 15), type="integer")
    checkTrue(is(sink, "TENxRealizationSink"))
    checkIdentical(dim(sink), c(9L, 15L))
    checkIdentical(chunkdim(sink), c(9L, 1L))
    checkIdentical(type(sink), "integer")

    set.seed(2009)
    svt <- poissonSparseMatrix(nrow(sink), ncol(sink), density=0.3)

    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 1), c(9, 2)))
    sink <- write_block(sink, viewport, svt[ , 1:2])
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 3), c(9, 5)))
    sink <- write_block(sink, viewport, as.matrix(svt[ , 3:5]))
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 6), c(9, 6)))
    sink <- write_block(sink, viewport, svt[ , 6, drop=FALSE])
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 7), c(9, 10)))
    sink <- write_block(sink, viewport, as(svt[ , 7:10], "dgCMatrix"))

    ## Trying to close too soon.
    checkException(as(sink, "TENxMatrixSeed"))

    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 11), c(9, 15)))
    sink <- write_block(sink, viewport, as(svt[ , 11:15], "COO_SparseMatrix"))

    seed <- as(sink, "TENxMatrixSeed")
    checkTrue(is(seed, "TENxMatrixSeed"))
    checkIdentical(dim(seed), dim(sink))
    checkIdentical(chunkdim(seed), chunkdim(sink))
    checkIdentical(type(seed), type(sink))

    expected_data <- nzvals(svt)
    Mindex <- nzwhich(svt, arr.ind=TRUE)

    data <- HDF5Array:::.read_h5sparse_data(path(seed), group(seed),
                                            NULL, type(seed))
    checkIdentical(data, expected_data)

    indices <- HDF5Array:::.read_h5sparse_indices(path(seed), group(seed))
    checkIdentical(indices, Mindex[ , 1L] - 1L)

    indptr <- HDF5Array:::.read_h5sparse_indptr(path(seed), group(seed))
    checkIdentical(rep.int(seq_len(ncol(seed)), diff(indptr)), Mindex[ , 2L])
}

test_from_TENxRealizationSink_to_TENxMatrixSeed <- function()
{
    for (type in c("integer", "logical", "double", "raw")) {
        for (H5type in c("H5T_NATIVE_UINT8", "H5T_NATIVE_UINT32")) {
            sink <- TENxRealizationSink(c(20, 0), type=type, H5type=H5type)
            checkTrue(is(sink, "TENxRealizationSink"))
            checkIdentical(dim(sink), c(20L, 0L))
            checkIdentical(chunkdim(sink), c(20L, 0L))
            checkIdentical(type(sink), type)

            seed <- as(sink, "TENxMatrixSeed")
            checkTrue(is(seed, "TENxMatrixSeed"))
            checkIdentical(dim(seed), dim(sink))
            checkIdentical(chunkdim(seed), chunkdim(sink))
            checkIdentical(type(seed), type(sink))
        }
    }
}

test_writeTENxMatrix <- function()
{
    set.seed(2009)
    svt <- poissonSparseMatrix(100, 375, density=0.4) * c(0.1, -0.1)
    m <- as.matrix(svt)
    dgc <- as(svt, "dgCMatrix")
    M <- writeHDF5Array(svt)

    objects <- list(svt, m, dgc, M)

    ## By default, blocks are much bigger than 'svt', so, BLOCK_write_to_sink(),
    ## the workhorse behind writeTENxMatrix(), would walk on a grid
    ## made of a single block that covers the full matrix-like object to
    ## write.
    ## So choosing a block size that will induce a grid made of 13 blocks.
    ## This can be seen by calling writeTENxMatrix() with 'verbose=TRUE'.
    setAutoBlockSize(25000)

    filepath <- tempfile(fileext=".h5")

    for (i in seq_along(objects)) {
        x <- objects[[i]]
        group <- paste0("/M", i)
        TENx <- writeTENxMatrix(x, filepath=filepath, group=group)
        checkTrue(is(TENx, "TENxMatrix"))
        checkIdentical(dim(TENx), dim(svt))
        checkIdentical(chunkdim(TENx), c(100L, 1L))
        checkIdentical(path(TENx), filepath)
        checkIdentical(group(TENx), group)
        checkIdentical(extract_sparse_array(TENx, list(NULL, NULL)), svt)
        checkIdentical(as.matrix(TENx), m)
    }
}

