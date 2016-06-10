test_DelayedArray_constructor <- function()
{
    ## 3-col seed

    DF1 <- DataFrame(aa=1:5, bb=letters[1:5], cc=Rle(c(-0.25, 7.1), 3:2))
    df1 <- as.data.frame(DF1)
    m1 <- as.matrix(df1)

    check_M1 <- function(M1) {
        checkTrue(class(M1) == "DelayedMatrix")
        checkTrue(validObject(M1, complete=TRUE))
        checkIdentical("character", HDF5Array:::type(M1))
        checkIdentical("character", HDF5Array:::type(M1[ , -2]))
        checkIdentical("character", HDF5Array:::type(M1[ , 0]))
        checkIdentical("character", HDF5Array:::type(M1[0 , ]))
        checkIdentical(dim(m1), dim(M1))
        checkIdentical(rownames(m1), rownames(M1))
        checkIdentical(colnames(m1), colnames(M1))
    }
    M1a <- DelayedArray(m1)   # matrix seed
    check_M1(M1a)

    M1b <- DelayedArray(df1)  # data.frame seed
    rownames(M1b) <- NULL
    check_M1(M1b)

    M1c <- DelayedArray(DF1)  # DataFrame seed
    check_M1(M1c)

    ## 2-col seed

    DF2 <- DF1[ , -2]
    df2 <- as.data.frame(DF2)
    m2 <- as.matrix(df2)

    check_M2 <- function(M2) {
        checkTrue(class(M2) == "DelayedMatrix")
        checkTrue(validObject(M2, complete=TRUE))
        checkIdentical("double", HDF5Array:::type(M2))
        checkIdentical("double", HDF5Array:::type(M2[ , -2]))
        checkIdentical("double", HDF5Array:::type(M2[ , 0]))
        checkIdentical("double", HDF5Array:::type(M2[0 , ]))
        checkIdentical(dim(m2), dim(M2))
        checkIdentical(rownames(m2), rownames(M2))
        checkIdentical(colnames(m2), colnames(M2))
    }

    M2a <- DelayedArray(m2)   # matrix seed
    check_M2(M2a)

    M2b <- DelayedArray(df2)  # data.frame seed
    rownames(M2b) <- NULL
    check_M2(M2b)

    M2c <- DelayedArray(DF2)  # DataFrame seed
    check_M2(M2c)
}

