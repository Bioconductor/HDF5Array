m2 <- matrix(runif(60), ncol=6)                        # numeric matrix

block_sizes1 <- c(12L, 20L, 50L, 10000L)
block_sizes2 <- 2L * block_sizes1

test_HDF5Matrix_row_col_summary <- function()
{
    test_row_col_summary <- function(FUN, m, block_sizes) {
        FUN <- match.fun(FUN)
        target1 <- FUN(m)
        target2 <- FUN(m, na.rm=TRUE)
        target3 <- FUN(t(m))
        target4 <- FUN(t(m), na.rm=TRUE)
        M <- as(m, "HDF5Matrix")
        for (block_size in block_sizes) {
            options(HDF5Array.block.size=block_size)
            checkIdentical(target1, FUN(M))
            checkIdentical(target2, FUN(M, na.rm=TRUE))
            checkIdentical(target3, FUN(t(M)))
            checkIdentical(target4, FUN(t(M), na.rm=TRUE))
        }
    }

    ## on a numeric matrix
    m <- m2
    m[2, 4] <- NA
    m[5, 4] <- Inf
    m[6, 3] <- -Inf
    for (FUN in c("rowSums", "colSums", "rowMeans", "colMeans"))
        test_row_col_summary(FUN, m, block_sizes2)

    library(genefilter)
    ## Note that the matrixStats package defines another rowVars() function.
    test_row_col_summary(genefilter::rowVars, m, block_sizes2)
}

