a1 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))  # integer array
a2 <- a1 + runif(150) - 0.5                              # numeric array
m2 <- a2[ , , 2]                                         # numeric matrix

block_sizes1 <- c(12L, 20L, 50L, 10000L)
block_sizes2 <- 2L * block_sizes1

test_Arith_and_Math_HDF5Array <- function()
{
    toto1 <- function(a) { 100 / floor(abs((5 * log(a + 0.2) - 1)^3)) }
    toto2 <- function(a) { 100L + (5L * (a - 2L)) %% 7L }

    ## with an integer array
    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- as(a, "HDF5Array")
    ## Not sure what's going on but it seems that this call to checkIdentical()
    ## crashes the RUnit package but only when the tests are run by
    ## 'R CMD check'.
    #checkIdentical(toto2(a), as.array(toto2(A)))

    ## with a numeric array
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- as(a, "HDF5Array")
    checkIdentical(toto1(a), as.array(toto1(A)))
    checkIdentical(toto2(a), as.array(toto2(A)))
    checkIdentical(toto2(toto1(a)), as.array(toto2(toto1(A))))
    checkIdentical(toto1(toto2(a)), as.array(toto1(toto2(A))))

    a <- a[ , 10:4, -2]
    A <- A[ , 10:4, -2]
    checkIdentical(toto1(a), as.array(toto1(A)))
    checkIdentical(toto2(a), as.array(toto2(A)))
    checkIdentical(toto2(toto1(a)), as.array(toto2(toto1(A))))
    checkIdentical(toto1(toto2(a)), as.array(toto1(toto2(A))))

    ## with a numeric matrix
    m <- m2
    m[5, 2] <- NA
    M <- as(m, "HDF5Matrix")
    checkIdentical(toto1(m), as.matrix(toto1(M)))
    checkIdentical(t(toto1(m)), as.matrix(toto1(t(M))))
    checkIdentical(t(toto1(m)), as.matrix(t(toto1(M))))
}

test_anyNA_HDF5Array <- function()
{
    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))

    HDF5Array_block_anyNA <- HDF5Array:::.HDF5Array_block_anyNA

    A1 <- as(a1, "HDF5Array")
    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- as(a, "HDF5Array")

    for (block_size in block_sizes2) {
        options(HDF5Array.block.size=block_size)
        checkIdentical(FALSE, anyNA(A1))
        checkIdentical(FALSE, HDF5Array_block_anyNA(a1))
        checkIdentical(TRUE, anyNA(A))
        checkIdentical(TRUE, HDF5Array_block_anyNA(a))
    }
}

test_Summary_HDF5Array <- function()
{
    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))

    test_Summary <- function(.Generic, a, block_sizes) {
        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(a)
        target2 <- GENERIC(a, na.rm=TRUE)
        A <- as(a, "HDF5Array")
        HDF5Array_block_Summary <- HDF5Array:::.HDF5Array_block_Summary
        for (block_size in block_sizes) {
            options(HDF5Array.block.size=block_size)
            checkIdentical(target1, GENERIC(A))
            checkIdentical(target1, GENERIC(t(A)))
            checkIdentical(target1, HDF5Array_block_Summary(.Generic, a))
            checkIdentical(target2, GENERIC(A, na.rm=TRUE))
            checkIdentical(target2, GENERIC(t(A), na.rm=TRUE))
            checkIdentical(target2, HDF5Array_block_Summary(.Generic, a,
                                                            na.rm=TRUE))
        }
    }

    ## on an integer array
    a <- a1
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    #for (.Generic in c("max", "min", "range", "sum", "prod")) {
    for (.Generic in c("max", "min", "range", "sum"))
        test_Summary(.Generic, a, block_sizes1)

    ## on a numeric array
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    a[2, 10, 2] <- Inf  # same as a[[97]] <- Inf
    for (.Generic in c("max", "min", "range", "sum", "prod"))
        test_Summary(.Generic, a, block_sizes2)

    ## on a logical array
    a <- array(c(rep(NA, 62), rep(TRUE, 87), FALSE), c(5, 10, 3))
    for (.Generic in c("any", "all"))
        test_Summary(.Generic, a, block_sizes1)
}

test_mean_HDF5Array <- function()
{
    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))

    ## on a numeric array
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- as(a, "HDF5Array")

    target1 <- mean(a)
    target2 <- mean(a, na.rm=TRUE)
    target3 <- mean(a[ , 10:4, -2])
    for (block_size in block_sizes2) {
        options(HDF5Array.block.size=block_size)
        checkIdentical(target1, mean(A))
        checkIdentical(target1, mean(t(A)))
        checkIdentical(target2, mean(A, na.rm=TRUE))
        checkIdentical(target2, mean(t(A), na.rm=TRUE))
        checkIdentical(target3, mean(A[ , 10:4, -2]))
        checkIdentical(target3, mean(t(A[ , 10:4, -2])))
    }
}

test_Compare_HDF5Array <- function()
{
    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))

    ## comparing an HDF5Array object with an atomic vector of length 1
    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- as(a, "HDF5Array")

    checkIdentical(a == a[[1]], as.array(A == A[[1]]))
    checkIdentical(a[[1]] == a, as.array(A[[1]] == A))

    checkIdentical(a != a[[1]], as.array(A != A[[1]]))
    checkIdentical(a[[1]] != a, as.array(A[[1]] != A))

    checkIdentical(a <= a[[1]], as.array(A <= A[[1]]))
    checkIdentical(a[[1]] >= a, as.array(A[[1]] >= A))

    checkIdentical(a >= a[[1]], as.array(A >= A[[1]]))
    checkIdentical(a[[1]] <= a, as.array(A[[1]] <= A))

    checkIdentical(a < a[[1]], as.array(A < A[[1]]))
    checkIdentical(a[[1]] > a, as.array(A[[1]] > A))

    checkIdentical(a > a[[1]], as.array(A > A[[1]]))
    checkIdentical(a[[1]] < a, as.array(A[[1]] < A))

    for (block_size in block_sizes2) {
            options(HDF5Array.block.size=block_size)
            checkIdentical(sum(a == a[[1]], na.rm=TRUE),
                           sum(A == A[[1]], na.rm=TRUE))
            checkIdentical(sum(a != a[[1]], na.rm=TRUE),
                           sum(A != A[[1]], na.rm=TRUE))
            checkIdentical(sum(a <= a[[1]], na.rm=TRUE),
                           sum(A <= A[[1]], na.rm=TRUE))
            checkIdentical(sum(a >= a[[1]], na.rm=TRUE),
                           sum(A >= A[[1]], na.rm=TRUE))
            checkIdentical(sum(a < a[[1]], na.rm=TRUE),
                           sum(A < A[[1]], na.rm=TRUE))
            checkIdentical(sum(a > a[[1]], na.rm=TRUE),
                           sum(A > A[[1]], na.rm=TRUE))
    }

    ## comparing 2 HDF5Array objects
    A1 <- as(a1, "HDF5Array")
    A2 <- as(a2, "HDF5Array")
    a3 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))
    a3[2, 9, 2] <- NA  # same as a3[[92]] <- NA
    A3 <- as(a3, "HDF5Array")

    for (.Generic in c("==", "!=", "<=", ">=", "<", ">")) {
        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(a1, a2)
        target2 <- GENERIC(a2, a1)
        target3 <- GENERIC(a1, a3)
        ## The smallest block size is very slow -> skip it.
        for (block_size in block_sizes2[-1L]) {
            options(HDF5Array.block.size=block_size)
            checkIdentical(target1, GENERIC(A1, A2))
            checkIdentical(target2, GENERIC(A2, A1))
            checkIdentical(target3, GENERIC(A1, A3))
        }
    }
}

