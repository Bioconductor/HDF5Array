m2 <- matrix(runif(60), ncol=6)                        # numeric matrix

block_sizes1 <- c(12L, 20L, 50L, 10000L)
block_sizes2 <- 2L * block_sizes1

test_HDF5Matrix_row_col_summary <- function()
{
    test_row_col_summary <- function(.Generic, m, block_sizes) {
        GENERIC <- match.fun(.Generic)
        target1 <- GENERIC(m)
        target2 <- GENERIC(m, na.rm=TRUE)
        target3 <- GENERIC(t(m))
        target4 <- GENERIC(t(m), na.rm=TRUE)
        M <- as(m, "HDF5Matrix")
        for (block_size in block_sizes) {
            options(HDF5Array.block.size=block_size)
            checkIdentical(target1, GENERIC(M))
            checkIdentical(target2, GENERIC(M, na.rm=TRUE))
            checkIdentical(target3, GENERIC(t(M)))
            checkIdentical(target4, GENERIC(t(M), na.rm=TRUE))
        }
    }

    ## on a numeric matrix
    m <- m2
    m[2, 4] <- NA
    m[5, 4] <- Inf
    m[6, 3] <- -Inf
    for (.Generic in c("rowSums", "colSums", "rowMeans", "colMeans"))
        test_row_col_summary(.Generic, m, block_sizes2)
}

