
test_TENxRealizationSink <- function()
{
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

    ## Trying to close prematurely.
    checkException(close(sink))

    ## Valid viewport.
    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 15), c(20, 15)))
    block <- matrix(0, nrow=nrow(viewport), ncol=ncol(viewport))
    sink <- write_block(sink, viewport, block)

    close(sink)

    seed <- as(sink, "TENxMatrixSeed")
    checkTrue(is(seed, "TENxMatrixSeed"))
    checkIdentical(dim(seed), dim(sink))
    checkIdentical(chunkdim(seed), chunkdim(sink))
    checkIdentical(type(seed), type(sink))

    data <- h5mread(path(seed), paste0(group(seed), "/data"))
    checkIdentical(data, numeric(0))
    indices <- h5mread(path(seed), paste0(group(seed), "/indices"))
    checkIdentical(indices, integer(0))
    indptr <- h5mread(path(seed), paste0(group(seed), "/indptr"))
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

    ## Trying to close prematurely.
    checkException(close(sink))

    viewport <- ArrayViewport(dim(sink), IRanges(c(1, 11), c(9, 15)))
    sink <- write_block(sink, viewport, as(svt[ , 11:15], "COO_SparseMatrix"))

    close(sink)

    seed <- as(sink, "TENxMatrixSeed")
    checkTrue(is(seed, "TENxMatrixSeed"))
    checkIdentical(dim(seed), dim(sink))
    checkIdentical(chunkdim(seed), chunkdim(sink))
    checkIdentical(type(seed), type(sink))

    expected_data <- nzvals(svt)
    Mindex <- nzwhich(svt, arr.ind=TRUE)
    data <- h5mread(path(seed), paste0(group(seed), "/data"))
    checkIdentical(data, expected_data)
    indices <- h5mread(path(seed), paste0(group(seed), "/indices"))
    checkIdentical(indices, Mindex[ , 1L] - 1L)
    indptr <- h5mread(path(seed), paste0(group(seed), "/indptr"))
    checkIdentical(rep.int(seq_len(ncol(seed)), diff(indptr)), Mindex[ , 2L])
}

test_writeTENxMatrix <- function()
{
    set.seed(2009)
    svt <- poissonSparseMatrix(100, 375, density=0.4) * 0.2 * c(1, -1)
    m <- as.matrix(svt)
    dgc <- as(svt, "dgCMatrix")

    ## By default, blocks are much bigger than 'svt', so, BLOCK_write_to_sink(),
    ## the workhorse behind writeTENxMatrix(), would walk on a grid of blocks
    ## made of a single block that covers the full matrix-like object to
    ## write.
    ## Choosing a block size that will induce a grid of blocks made
    ## of 13 blocks. This can be seen by calling writeTENxMatrix()
    ## with 'verbose=TRUE'.
    setAutoBlockSize(25000)

    M0 <- writeTENxMatrix(svt)
    checkIdentical(extract_sparse_array(M0, list(NULL, NULL)), svt)
    checkIdentical(as.matrix(M0), m)

    M1 <- writeTENxMatrix(m)
    checkIdentical(extract_sparse_array(M1, list(NULL, NULL)), svt)
    checkIdentical(as.matrix(M1), m)

    M2 <- writeTENxMatrix(dgc)
    checkIdentical(extract_sparse_array(M2, list(NULL, NULL)), svt)
    checkIdentical(as.matrix(M2), m)
}

