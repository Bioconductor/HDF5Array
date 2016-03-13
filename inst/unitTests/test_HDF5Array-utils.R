
test_anyNA_HDF5Array <- function()
{
    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))

    HDF5Array_block_anyNA <- HDF5Array:::.HDF5Array_block_anyNA

    a1 <- a2 <- array(1:300, c(5, 10, 2, 3))
    a2[2, 9, 2, 1] <- NA  ## same as a2[[92]] <- NA
    A1 <- as(a1, "HDF5Array")
    A2 <- as(a2, "HDF5Array")

    target1 <- anyNA(a1)
    target2 <- anyNA(a2)
    for (block_size in c(12L, 20L, 50L, 10000L)) {
        options(HDF5Array.block.size=block_size)

        current <- anyNA(A1)
        checkIdentical(target1, current)
        current <- HDF5Array_block_anyNA(a1)
        checkIdentical(target1, current)

        current <- anyNA(A2)
        checkIdentical(target2, current)
        current <- HDF5Array_block_anyNA(a2)
        checkIdentical(target2, current)
    }
}

test_Summary_HDF5Array <- function()
{
    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))

    HDF5Array_block_Summary <- HDF5Array:::.HDF5Array_block_Summary

    a1 <- array(1:300, c(5, 10, 2, 3))
    a1[2, 9, 2, 1] <- NA  ## same as a1[[92]] <- NA
    a2 <- a1
    storage.mode(a2) <- "double"
    a2[2, 10, 2, 1] <- Inf  ## same as a2[[97]] <- Inf
    A1 <- as(a1, "HDF5Array")
    A2 <- as(a2, "HDF5Array")

    for (.Generic in c("max", "min", "range", "sum", "prod")) {
        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(a1)
        target1narm <- GENERIC(a1, na.rm=TRUE)
        target2 <- GENERIC(a2)
        target2narm <- GENERIC(a2, na.rm=TRUE)
        for (block_size in c(12L, 20L, 50L, 10000L)) {
            options(HDF5Array.block.size=block_size)

            current <- GENERIC(A1)
            checkIdentical(target1, current)
            current <- HDF5Array_block_Summary(.Generic, a1)
            checkIdentical(target1, current)

            current <- GENERIC(A1, na.rm=TRUE)
            checkIdentical(target1narm, current)
            current <- HDF5Array_block_Summary(.Generic, a1, na.rm=TRUE)
            checkIdentical(target1narm, current)

            current <- GENERIC(A2)
            checkIdentical(target2, current)
            current <- HDF5Array_block_Summary(.Generic, a2)
            checkIdentical(target2, current)

            current <- GENERIC(A2, na.rm=TRUE)
            checkIdentical(target2narm, current)
            current <- HDF5Array_block_Summary(.Generic, a2, na.rm=TRUE)
            checkIdentical(target2narm, current)
        }
    }

    a1 <- array(c(rep(NA, 62), rep(TRUE, 237), FALSE), c(5, 10, 2, 3))
    A1 <- as(a1, "HDF5Array")

    for (.Generic in c("any", "all")) {
        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(a1)
        target1narm <- GENERIC(a1, na.rm=TRUE)
        for (block_size in c(12L, 20L, 50L, 10000L)) {
            options(HDF5Array.block.size=block_size)

            current <- GENERIC(A1)
            checkIdentical(target1, current)
            current <- HDF5Array_block_Summary(.Generic, a1)
            checkIdentical(target1, current)

            current <- GENERIC(A1, na.rm=TRUE)
            checkIdentical(target1narm, current)
            current <- HDF5Array_block_Summary(.Generic, a1, na.rm=TRUE)
            checkIdentical(target1narm, current)
        }
    }
}

