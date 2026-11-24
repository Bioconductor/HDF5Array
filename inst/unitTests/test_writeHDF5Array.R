
test_HDF5RealizationSink <- function()
{
    checkTrue(extends("HDF5RealizationSink", "RealizationSink"))

    ## --- Some edge cases (0 columns or 0 rows) ---

    for (type in c("integer", "logical", "double", "raw", "character")) {
        sink <- HDF5RealizationSink(c(20, 0),  # 0 columns
                                    type=type, chunkdim=c(12, 0))
        checkTrue(is(sink, "HDF5RealizationSink"))
        checkIdentical(dim(sink), c(20L, 0L))
        checkIdentical(chunkdim(sink), c(12L, 0L))
        checkIdentical(type(sink), type)

        seed <- as(sink, "HDF5ArraySeed")
        checkTrue(is(seed, "HDF5ArraySeed"))
        checkIdentical(dim(seed), dim(sink))
        checkIdentical(chunkdim(seed), chunkdim(sink))
        checkIdentical(type(seed), type(sink))

        sink <- HDF5RealizationSink(c(0, 15),  # 0 rows
                                    type=type, chunkdim=c(0, 6))
        checkTrue(is(sink, "HDF5RealizationSink"))
        checkIdentical(dim(sink), c(0L, 15L))
        checkIdentical(chunkdim(sink), c(0L, 6L))
        checkIdentical(type(sink), type)

        seed <- as(sink, "HDF5ArraySeed")
        checkTrue(is(seed, "HDF5ArraySeed"))
        checkIdentical(dim(seed), dim(sink))
        checkIdentical(chunkdim(seed), chunkdim(sink))
        checkIdentical(type(seed), type(sink))
    }

    ## --- Create a HDF5RealizationSink and write a few blocks ---
    ## ---                    of data to it                    ---

    sink_dim <- c(9L, 15L, 2L)
    viewports <- list(
        ArrayViewport(sink_dim, IRanges(c(1, 11, 1), c(5, 15, 2))),
        ArrayViewport(sink_dim, IRanges(c(8, 1, 2), c(9, 2, 2))),
        ArrayViewport(sink_dim, IRanges(c(4, 3, 1), c(8, 12, 2))),
        ArrayViewport(sink_dim, IRanges(c(9, 5, 1), c(9, 6, 2))),
        ArrayViewport(sink_dim, IRanges(c(1, 8, 1), c(9, 11, 2))),
        ArrayViewport(sink_dim, IRanges(c(8, 2, 1), c(9, 3, 1))),
        ArrayViewport(sink_dim, IRanges(c(1, 1, 1), c(4, 3, 2)))
    )

    for (type in c("integer", "double")) {
        for (chunkdim in list(c(5L, 6L, 1L), c(2L, 10L, 2L))) {
            sink <- HDF5RealizationSink(sink_dim, type=type, chunkdim=chunkdim)
            checkTrue(is(sink, "HDF5RealizationSink"))
            checkIdentical(dim(sink), sink_dim)
            checkIdentical(chunkdim(sink), chunkdim)
            checkIdentical(type(sink), type)
            expected <- array(vector(type, length=1), sink_dim)
            checkIdentical(as.array(as(sink, "HDF5ArraySeed")), expected)
            ## Write a few blocks of data to 'sink'.
            for (i in seq_along(viewports)) {
                viewport <- viewports[[i]]
                data <- 100L*i + seq_along(viewport) + vector(type, length=1)
                block <- array(data, dim(viewport))
                sink <- write_block(sink, viewport, block)
                seed <- as(sink, "HDF5ArraySeed")
                checkIdentical(read_block(seed, viewport), block)
                expected <- write_block(expected, viewport, block)
                checkIdentical(as.array(seed), expected)
            }
        }
    }
}

test_from_HDF5RealizationSink_to_HDF5ArraySeed <- function()
{
    for (type in c("integer", "logical", "double", "raw")) {
        for (H5type in c("H5T_NATIVE_UINT8", "H5T_NATIVE_UINT32")) {
            sink <- HDF5RealizationSink(c(20, 10), type=type, H5type=H5type)
            checkTrue(is(sink, "HDF5RealizationSink"))
            checkIdentical(dim(sink), c(20L, 10L))
            checkIdentical(type(sink), type)

            seed <- as(sink, "HDF5ArraySeed")
            checkTrue(is(seed, "HDF5ArraySeed"))
            checkIdentical(dim(seed), dim(sink))
            checkIdentical(chunkdim(seed), chunkdim(sink))
            checkIdentical(type(seed), type(sink))
        }
    }
}

## 'a0' should be an ordinary matrix or array, or a dgCMatrix object.
.basic_checks <- function(A, a0)
{
    checkTrue(validObject(A))
    checkTrue(is(A, "HDF5Array"))
    checkIdentical(dim(A), dim(a0))
    checkIdentical(type(A), type(a0))
    ## dimnames() of a dgCMatrix is list(NULL, NULL) but turning 'a0' into
    ## an ordinary array fixes that so we do it **before** comparing
    ## dimnames(A) with dimnames(a).
    a <- as.array(a0)
    checkEquals(dimnames(A), dimnames(a))
    checkIdentical(as.array(A), a)
}

test_writeHDF5Array_2D <- function()
{
    m0 <- matrix(1700:1001, ncol=7, dimnames=list(NULL, letters[1:7]))

    M1 <- writeHDF5Array(m0)
    .basic_checks(M1, m0)
    checkTrue(class(M1) == "HDF5Matrix")
    checkEquals(chunkdim(M1), dim(m0))
    checkTrue(!is_sparse(M1))

    M2 <- writeHDF5Array(M1, chunkdim=c(30, 5))
    .basic_checks(M2, m0)
    checkTrue(class(M2) == "HDF5Matrix")
    checkEquals(chunkdim(M2), c(30, 5))
    checkTrue(!is_sparse(M2))

    M3 <- writeHDF5Array(M2, chunkdim=c(10, 7), as.sparse=TRUE)
    .basic_checks(M3, m0)
    checkTrue(class(M3) == "HDF5Matrix")
    checkEquals(chunkdim(M3), c(10, 7))
    checkTrue(is_sparse(M3))

    M4 <- writeHDF5Array(M3, with.dimnames=FALSE)
    .basic_checks(M4, unname(m0))
    checkTrue(class(M4) == "HDF5Matrix")
    checkEquals(chunkdim(M4), chunkdim(M3))
    checkTrue(is_sparse(M4))

    set.seed(2009L)
    dgc <- rsparsematrix(900, 150, density=0.05)
    setAutoBlockSize(77)
    M5 <- writeHDF5Array(dgc, chunkdim=c(50, 50))
    setAutoBlockSize()
    .basic_checks(M5, dgc)
    checkTrue(class(M5) == "HDF5Matrix")
    checkEquals(chunkdim(M5), c(50, 50))
    checkTrue(is_sparse(M5))
}

test_writeHDF5Array_4D <- function()
{
    a0 <- array(1:3000, c(6:4, 25),
                dimnames=list(letters[1:6], NULL, NULL, LETTERS[1:25]))

    A1 <- writeHDF5Array(a0)
    .basic_checks(A1, a0)
    checkTrue(class(A1) == "HDF5Array")
    checkEquals(chunkdim(A1), dim(a0))
    checkTrue(!is_sparse(A1))

    setAutoBlockSize(77)
    A2 <- writeHDF5Array(A1, chunkdim=c(2, 5, 2, 5))
    setAutoBlockSize()
    .basic_checks(A2, a0)
    checkTrue(class(A2) == "HDF5Array")
    checkEquals(chunkdim(A2), c(2, 5, 2, 5))
    checkTrue(!is_sparse(A2))

    A3 <- writeHDF5Array(A2, as.sparse=TRUE)
    .basic_checks(A3, a0)
    checkTrue(class(A3) == "HDF5Array")
    checkEquals(chunkdim(A3), chunkdim(A2))
    checkTrue(is_sparse(A3))

    A4 <- writeHDF5Array(A3, with.dimnames=FALSE)
    .basic_checks(A4, unname(a0))
    checkTrue(class(A4) == "HDF5Array")
    checkEquals(chunkdim(A4), chunkdim(A3))
    checkTrue(is_sparse(A4))
}

test_writeHDF5Array_on_workers <- function()
{
    library(BiocParallel)
    snow2 <- SnowParam(workers=2)
    make_mat <- function(i) matrix(100*i + 1:12, nrow=2)
    res <- bplapply(1:5, function(i) HDF5Array::writeHDF5Array(make_mat(i)),
                    BPPARAM=snow2)
    for (i in seq_along(res)) {
        M <- res[[i]]
        m <- make_mat(i)
        .basic_checks(M, m)
        checkTrue(class(M) == "HDF5Matrix")
    }
}

