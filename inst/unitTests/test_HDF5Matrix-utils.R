Arith_members <- c("+", "-", "*", "/", "^", "%%", "%/%")
Compare_members <- c("==", "!=", "<=", ">=", "<", ">")
Logic_members <- c("&", "|")  # currently untested

a1 <- array(sample(5L, 150, replace=TRUE), c(5, 10, 3))  # integer array
a2 <- a1 + runif(150) - 0.5                              # numeric array
m2 <- matrix(runif(60), ncol=6)                          # numeric matrix

block_sizes1 <- c(12L, 20L, 50L, 10000L)
block_sizes2 <- 2L * block_sizes1

test_HDF5Matrix_delayed_Ops <- function()
{
    test_delayed_Ops_on_matrix <- function(.Generic, m, M) {
        GENERIC <- match.fun(.Generic)

        target_current <- list(
            list(GENERIC(m, m[ , 1]), GENERIC(M, M[ , 1])),
            list(GENERIC(m[ , 2], m), GENERIC(M[ , 2], M))
        )
        for (i in seq_along(target_current)) {
            target <- target_current[[i]][[1L]]
            current <- target_current[[i]][[2L]]
            checkIdentical(target, as.matrix(current))
            checkIdentical(t(target), as.matrix(t(current)))
            checkIdentical(target[-2, 8:5], as.matrix(current[-2, 8:5]))
            checkIdentical(t(target[-2, 8:5]), as.matrix(t(current[-2, 8:5])))
            checkIdentical(target[-2, 0], as.matrix(current[-2, 0]))
            checkIdentical(t(target[-2, 0]), as.matrix(t(current[-2, 0])))
            checkIdentical(target[0, ], as.matrix(current[0, ]))
            checkIdentical(t(target[0, ]), as.matrix(t(current[0, ])))
        }

        target_current <- list(
            list(GENERIC(t(m), 8:-1), GENERIC(t(M), 8:-1)),
            list(GENERIC(8:-1, t(m)), GENERIC(8:-1, t(M))),

            list(GENERIC(t(m), m[1 , ]), GENERIC(t(M), M[1 , ])),
            list(GENERIC(m[2 , ], t(m)), GENERIC(M[2 , ], t(M))),

            list(GENERIC(t(m), m[1 , 6:10]), GENERIC(t(M), M[1 , 6:10])),
            list(GENERIC(m[2 , 8:7], t(m)), GENERIC(M[2 , 8:7], t(M)))
        )
        for (i in seq_along(target_current)) {
            target <- target_current[[i]][[1L]]
            current <- target_current[[i]][[2L]]
            checkIdentical(target, as.matrix(current))
            checkIdentical(target[1:3 , ], as.matrix(current[1:3 , ]))
            checkIdentical(target[ , 1:3], as.matrix(current[ , 1:3]))
            checkIdentical(t(target), as.matrix(t(current)))
            checkIdentical(t(target)[1:3 , ], as.matrix(t(current)[1:3 , ]))
            checkIdentical(t(target)[ , 1:3], as.matrix(t(current)[ , 1:3]))
            checkIdentical(target[8:5, -2], as.matrix(current[8:5, -2]))
            checkIdentical(t(target[8:5, -2]), as.matrix(t(current[8:5, -2])))
            checkIdentical(target[0, -2], as.matrix(current[0, -2]))
            checkIdentical(t(target[0, -2]), as.matrix(t(current[0, -2])))
            checkIdentical(target[ , 0], as.matrix(current[ , 0]))
            checkIdentical(t(target[ , 0]), as.matrix(t(current[ , 0])))
        }
    }

    a <- a2
    a[2, 9, 2] <- NA  # same as a[[92]] <- NA
    A <- as(a, "HDF5Array")
    m <- a[ , , 2]
    M <- as(m, "HDF5Matrix")

    ## "Logic" members currently untested.
    for (.Generic in c(Arith_members, Compare_members))
        test_delayed_Ops_on_matrix(.Generic, m, M)

    M <- as(A[ , , 2], "HDF5Matrix")
    for (.Generic in c(Arith_members, Compare_members))
        test_delayed_Ops_on_matrix(.Generic, m, M)
}

test_HDF5Matrix_row_col_summary <- function()
{
    test_row_col_summary <- function(FUN, m, block_sizes) {
        on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))
        FUN <- match.fun(FUN)

        target1 <- FUN(m)
        target2 <- FUN(m, na.rm=TRUE)
        target3 <- FUN(t(m))
        target4 <- FUN(t(m), na.rm=TRUE)
        M <- as(m, "HDF5Matrix")
        for (block_size in block_sizes) {
            options(HDF5Array.block.size=block_size)
            checkEquals(target1, FUN(M))
            checkEquals(target2, FUN(M, na.rm=TRUE))
            checkEquals(target3, FUN(t(M)))
            checkEquals(target4, FUN(t(M), na.rm=TRUE))
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
    ## Note that the matrixStats package also defines a rowVars() function.
    test_row_col_summary(genefilter::rowVars, m, block_sizes2)
}

test_HDF5Matrix_mult <- function()
{
    m <- m2
    m[2, 4] <- NA
    m[5, 4] <- Inf
    m[6, 3] <- -Inf
    M <- as(m, "HDF5Matrix")

    Lm <- rbind(rep(1L, 10), rep(c(1L, 0L), 5), rep(-100L, 10))
    Rm <- rbind(Lm + 7.05, 0.1 * Lm)

    on.exit(options(HDF5Array.block.size=HDF5Array:::DEFAULT_BLOCK_SIZE))
    for (block_size in block_sizes2) {
        options(HDF5Array.block.size=block_size)
        P <- Lm %*% M
        checkEquals(Lm %*% m, as.matrix(P))
        P <- M %*% Rm
        checkEquals(m %*% Rm, as.matrix(P))
    }
}

