## 'a0' should be an ordinary matrix or array, or a dgCMatrix object.
.basic_checks <- function(A, a0)
{
    checkTrue(validObject(A))
    checkTrue(is(A, "HDF5Array"))
    checkEquals(dim(A), dim(a0))
    checkEquals(type(A), type(a0))
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
    checkEquals(chunkdim(M4), dim(m0))
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

    A3 <- writeHDF5Array(A2, chunkdim=c(3, 5, 4, 2), as.sparse=TRUE)
    .basic_checks(A3, a0)
    checkTrue(class(A3) == "HDF5Array")
    checkEquals(chunkdim(A3), c(3, 5, 4, 2))
    checkTrue(is_sparse(A3))

    A4 <- writeHDF5Array(A3, with.dimnames=FALSE)
    .basic_checks(A4, unname(a0))
    checkTrue(class(A4) == "HDF5Array")
    checkEquals(chunkdim(A4), dim(a0))
    checkTrue(is_sparse(A4))
}

