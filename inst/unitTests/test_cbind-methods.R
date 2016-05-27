test_DelayedMatrix_rbind_cbind <- function()
{
    m1 <- matrix(1:15, nrow=3, ncol=5,
                 dimnames=list(NULL, paste0("M1y", 1:5)))
    m2 <- matrix(101:135, nrow=7, ncol=5,
                 dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5)))
    m3 <- matrix(1001:1025, nrow=5, ncol=5,
                 dimnames=list(paste0("M3x", 1:5), NULL))
    M1 <- HDF5Array(m1)
    M2 <- HDF5Array(m2)
    M3 <- HDF5Array(m3)

    target <- rbind(a=m1, b=m2, c=m3)
    current <- rbind(a=M1, b=M2, c=M3)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    current <- cbind(a=t(M1), b=t(M2), c=t(M3))
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(t(target), as.matrix(current))

    ## unary form

    target <- rbind(a=m1)
    current <- rbind(a=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    target <- cbind(a=m1)
    current <- cbind(a=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    ## with empty matrices

    m1 <- matrix(nrow=0, ncol=3, dimnames=list(NULL, letters[1:3]))
    m2 <- matrix(1:15, ncol=3, dimnames=list(NULL, LETTERS[1:3]))
    M1 <- HDF5Array(m1)
    M2 <- HDF5Array(m2)

    target <- rbind(a=m1, a=m2)
    current <- rbind(a=M1, b=M2)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    target <- rbind(a=m2, a=m1)
    current <- rbind(a=M2, b=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))

    target <- rbind(a=m1, a=m1)
    current <- rbind(a=M1, b=M1)
    checkTrue(is(current, "DelayedMatrix"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, as.matrix(current))
}

