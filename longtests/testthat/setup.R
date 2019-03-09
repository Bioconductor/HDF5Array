# Setting up the test suite.

testpkg <- system.file("testpkg", package="beachmat")
devtools::install(testpkg, quick=TRUE)
library(beachtest)

# Setting up functions for generating the matrices to use for testing.

library(HDF5Array)
integer_hFUN <- function(nr=15, nc=10) {
    as(matrix(rpois(nr*nc, lambda=5), nr, nc), "HDF5Array")
}

logical_hFUN <- function(nr=15, nc=10) {
    as(matrix(rbinom(nr*nc, 1, 0.5)==1L, nr, nc), "HDF5Array")
}

numeric_hFUN <- function(nr=15, nc=10) {
    as(matrix(rnorm(nr*nc), nr, nc), "HDF5Array")
}

genwords <- function(n = 5000) {
    all.choices <- c(rep("", 4), LETTERS) # to get variable length strings.
    a <- do.call(paste0, replicate(5, sample(all.choices, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

character_hFUN <- function(nr=15, nc=10) {
    as(matrix(genwords(nr*nc), nr, nc), "HDF5Array")
}

check_write_HDF5 <- function(FUN, mode) 
# Checking that the output function in '.Call' does not overwrite the 
# underlying HDF5 file and change the values of other HDF5Matrix objects. 
{
    mats <- list(FUN(), FUN(5, 30), FUN(30, 5))
    ref.mats <- lapply(mats, as.matrix)
    all.paths <- sapply(mats, path)
    new.mats <- vector("list", length(mats)) 

    for (i in seq_along(mats)) { 
        out <- .Call(paste0("set_col_all_", mode), mats[[i]], seq_len(ncol(mats[[i]])), package="beachtest")[[1]]
        new.mats[[i]] <- out
        expect_s4_class(out, "HDF5Matrix")
        expect_identical(ref.mats[[i]], as.matrix(out))

        # Checking that all the other matrices are still what they should be.
        for (j in seq_along(mats)) { 
            ref <- ref.mats[[j]]
            expect_identical(ref, as.matrix(mats[[j]]))
        }

        all.paths <- c(all.paths, BiocGenerics::path(out))
    }
    expect_false(any(duplicated(all.paths)))

    # Checking that the old and realized files are in the log.
    expect_message(log <- showHDF5DumpLog())
    for (is.new in c(TRUE, FALSE)) {
        for (i in seq_along(mats)) {
            if (is.new) {
                current <- new.mats[[i]]
            } else {
                current <- mats[[i]]
            }   

            j <- which(log$name==seed(current)@name & log$file==path(current))
            expect_true(length(j)==1L)
            expect_identical(if (mode=="numeric") "double" else mode, log$type[j])
            expect_identical(sprintf("%ix%i", nrow(current), ncol(current)), log$dims[j])
        }
    }
    return(invisible(NULL))
}
