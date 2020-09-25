test_h5mread_2D <- function()
{
    do_2D_tests <- function(m, M, noreduce=FALSE, as.integer=FALSE, method=0L) {
        read <- function(starts=NULL, counts=NULL)
            h5mread(M@seed@filepath, M@seed@name,
                    starts=starts, counts=counts,
                    noreduce=noreduce, as.integer=as.integer, method=method)

        current <- read()
        checkIdentical(m, current)

        current <- read(list(NULL, NULL))
        checkIdentical(m, current)

        current <- read(list(integer(0), integer(0)))
        checkIdentical(m[NULL, NULL, drop=FALSE], current)

        ## With indices strictly sorted

        current <- read(list(c(2:5, 7:10), NULL))
        checkIdentical(m[c(2:5, 7:10), , drop=FALSE], current)

        current <- read(list(NULL, 1:2))
        checkIdentical(m[ , 1:2, drop=FALSE], current)

        current <- read(list(7:10, c(1:2, 5)))
        checkIdentical(m[7:10, c(1:2, 5), drop=FALSE], current)

        ## With indices in any order and with duplicates

        i <- c(2:6, 6:3, 1, 1, 1, 9:8)

        current <- read(list(i, NULL))
        checkIdentical(m[i, , drop=FALSE], current)

        current <- read(list(i, c(6:5, 5)))
        checkIdentical(m[i, c(6:5, 5), drop=FALSE], current)

        ## Only methods 1 and 3 support 'counts'.
        if (!(method %in% c(1, 3)))
            return()

        starts <- list(integer(0), 4L)
        counts <- list(integer(0), 2L)
        current <- read(starts, counts)
        checkIdentical(m[NULL, 4:5, drop=FALSE], current)

        starts <- list(5L, integer(0))
        counts <- list(4L, integer(0))
        current <- read(starts, counts)
        checkIdentical(m[5:8, NULL, drop=FALSE], current)

        starts <- list(c(2L, 5L), 4L)
        counts <- list(c(3L, 4L), 2L)
        current <- read(starts, counts)
        checkIdentical(m[2:8, 4:5, drop=FALSE], current)
    }

    do_2D_sparse_tests <- function(M, as.integer=FALSE) {
        test_with <- function(starts=NULL) {
            sas <- h5mread(M@seed@filepath, M@seed@name,
                           starts=starts,
                           as.integer=as.integer, as.sparse=TRUE)
            checkTrue(is(sas, "SparseArraySeed"))
            current <- sparse2dense(sas)
            target <- h5mread(M@seed@filepath, M@seed@name,
                              starts=starts, as.integer=as.integer)
            checkIdentical(target, current)
        }
        test_with()
        test_with(list(NULL, NULL))
	test_with(list(integer(0), NULL))
	test_with(list(NULL, integer(0)))
	test_with(list(integer(0), integer(0)))
        test_with(list(c(2:5, 7:10), NULL))
        test_with(list(NULL, 1:2))
        test_with(list(7:10, c(1:2, 5)))
    }

    chunkdims <- list(0,         # no chunking (i.e. contiguous data)
                      c(10, 6),
                      c(10, 1),
                      c( 1, 6),
                      c( 4, 5),
                      c( 4, 3),
                      c( 2, 2))

    ## with an integer matrix

    m0 <- matrix(1:60, ncol=6)
    for (chunkdim in chunkdims) {
        M0 <- writeHDF5Array(m0, chunkdim=chunkdim)
        do_2D_tests(m0, M0, method=1L)
        do_2D_tests(m0, M0, noreduce=TRUE, method=1L)
        do_2D_tests(m0, M0, method=2L)
        do_2D_tests(m0, M0, method=3L)
        do_2D_tests(m0, M0, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m0, M0, method=4L)
            do_2D_tests(m0, M0, method=6L)
            do_2D_tests(m0, M0, method=7L)
            do_2D_sparse_tests(M0)
        }
        do_2D_tests(m0, M0)
    }

    ## with a logical matrix

    m1 <- m0 %% 3L == 0L
    for (chunkdim in chunkdims) {
        M1 <- writeHDF5Array(m1, chunkdim=chunkdim)
        do_2D_tests(m1, M1, method=1L)
        do_2D_tests(m1, M1, noreduce=TRUE, method=1L)
        do_2D_tests(m1, M1, method=2L)
        do_2D_tests(m1, M1, method=3L)
        do_2D_tests(m1, M1, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m1, M1, method=4L)
            do_2D_tests(m1, M1, method=6L)
            do_2D_tests(m1, M1, method=7L)
            do_2D_sparse_tests(M1)
        }
        do_2D_tests(m1, M1)
        storage.mode(m1) <- "integer"
        do_2D_tests(m1, M1, as.integer=TRUE)
        do_2D_tests(m1, M1, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m1, M1, as.integer=TRUE, method=4L)
            do_2D_tests(m1, M1, as.integer=TRUE, method=6L)
            do_2D_tests(m1, M1, as.integer=TRUE, method=7L)
            do_2D_sparse_tests(M1, as.integer=TRUE)
        }
    }

    ## with a numeric matrix

    m2 <- matrix(10 * runif(60), ncol=6)
    for (chunkdim in chunkdims) {
        M2 <- writeHDF5Array(m2, chunkdim=chunkdim)
        do_2D_tests(m2, M2, method=1L)
        do_2D_tests(m2, M2, noreduce=TRUE, method=1L)
        do_2D_tests(m2, M2, method=2L)
        do_2D_tests(m2, M2, method=3L)
        do_2D_tests(m2, M2, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m2, M2, method=4L)
            do_2D_tests(m2, M2, method=6L)
            do_2D_tests(m2, M2, method=7L)
            do_2D_sparse_tests(M2)
        }
        do_2D_tests(m2, M2)
        storage.mode(m2) <- "integer"
        do_2D_tests(m2, M2, as.integer=TRUE)
        do_2D_tests(m2, M2, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m2, M2, as.integer=TRUE, method=4L)
            do_2D_tests(m2, M2, as.integer=TRUE, method=6L)
            do_2D_tests(m2, M2, as.integer=TRUE, method=7L)
            do_2D_sparse_tests(M2, as.integer=TRUE)
        }
    }

    ## with a character matrix

    m3 <- matrix(as.character(1:60), ncol=6)
    for (chunkdim in chunkdims) {
        M3 <- writeHDF5Array(m3, chunkdim=chunkdim)
        do_2D_tests(m3, M3, method=4L)
        if (!identical(chunkdim, 0)) {
            do_2D_sparse_tests(M3)
        }
    }

    m3[cbind(5:10, 6:1)] <- NA_character_
    for (chunkdim in chunkdims) {
        M3 <- writeHDF5Array(m3, chunkdim=chunkdim)
        do_2D_tests(m3, M3, method=4L)
        if (!identical(chunkdim, 0)) {
            do_2D_sparse_tests(M3)
        }
    }

    ## with a raw matrix

    m4 <- m0
    storage.mode(m4) <- "raw"
    for (chunkdim in chunkdims) {
        M4 <- writeHDF5Array(m4, chunkdim=chunkdim)
        do_2D_tests(m4, M4, method=1L)
        do_2D_tests(m4, M4, noreduce=TRUE, method=1L)
        do_2D_tests(m4, M4, method=2L)
        do_2D_tests(m4, M4, method=3L)
        do_2D_tests(m4, M4, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m4, M4, method=4L)
            do_2D_tests(m4, M4, method=6L)
            do_2D_tests(m4, M4, method=7L)
            do_2D_sparse_tests(M4)
        }
        do_2D_tests(m4, M4)
        do_2D_tests(m0, M4, as.integer=TRUE)
        do_2D_tests(m0, M4, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_2D_tests(m0, M4, as.integer=TRUE, method=4L)
            do_2D_tests(m0, M4, as.integer=TRUE, method=6L)
            do_2D_tests(m0, M4, as.integer=TRUE, method=7L)
            do_2D_sparse_tests(M4, as.integer=TRUE)
        }
    }
}

test_h5mread_3D <- function()
{
    DIM <- c(10, 15, 6)

    do_3D_tests <- function(a, A, noreduce=FALSE, as.integer=FALSE, method=0L) {
        read <- function(starts=NULL, counts=NULL)
            h5mread(A@seed@filepath, A@seed@name,
                    starts=starts, counts=counts,
                    noreduce=noreduce, as.integer=as.integer, method=method)

        current <- read()
        checkIdentical(a, current)

        current <- read(list(NULL, NULL, NULL))
        checkIdentical(a, current)

        current <- read(list(integer(0), integer(0), integer(0)))
        checkIdentical(a[NULL, NULL, NULL, drop=FALSE], current)

        ## With indices strictly sorted

        current <- read(list(c(2:5, 7:10), NULL, NULL))
        checkIdentical(a[c(2:5, 7:10), , , drop=FALSE], current)

        current <- read(list(NULL, c(11:12, 14), NULL))
        checkIdentical(a[ , c(11:12, 14), , drop=FALSE], current)

        current <- read(list(NULL, NULL, 1:2))
        checkIdentical(a[ , , 1:2, drop=FALSE], current)

        current <- read(list(7:10, c(11:12, 14), c(1:2, 5)))
        checkIdentical(a[7:10, c(11:12, 14), c(1:2, 5), drop=FALSE], current)

        ## With indices in any order and with duplicates

        i <- c(2:6, 6:3, 1, 1, 1, 9:8)

        current <- read(list(i, NULL, NULL))
        checkIdentical(a[i, , , drop=FALSE], current)

        current <- read(list(i, NULL, c(6:5, 5)))
        checkIdentical(a[i, , c(6:5, 5), drop=FALSE], current)

        ## Only methods 1 and 3 support 'counts'.
        if (!(method %in% c(1, 3)))
            return()

        starts <- list(integer(0), NULL, 4L)
        counts <- list(integer(0), NULL, 2L)
        current <- read(starts, counts)
        checkIdentical(a[NULL, , 4:5, drop=FALSE], current)

        starts <- list(5L, integer(0), NULL)
        counts <- list(4L, integer(0), NULL)
        current <- read(starts, counts)
        checkIdentical(a[5:8, NULL, , drop=FALSE], current)

        starts <- list(c(2L, 5L), NULL, 4L)
        counts <- list(c(3L, 4L), NULL, 2L)
        current <- read(starts, counts)
        checkIdentical(a[2:8, , 4:5, drop=FALSE], current)

        starts <- list(c(2L, 5L), 11L, 4L)
        counts <- list(c(3L, 4L),  5L, 2L)
        current <- read(starts, counts)
        checkIdentical(a[2:8, 11:15, 4:5, drop=FALSE], current)
    }

    do_3D_sparse_tests <- function(A, as.integer=FALSE) {
        test_with <- function(starts=NULL) {
            sas <- h5mread(A@seed@filepath, A@seed@name,
                           starts=starts,
                           as.integer=as.integer, as.sparse=TRUE)
            checkTrue(is(sas, "SparseArraySeed"))
            current <- sparse2dense(sas)
            target <- h5mread(A@seed@filepath, A@seed@name,
                              starts=starts, as.integer=as.integer)
            checkIdentical(target, current)
        }
        test_with()
        test_with(list(NULL, NULL, NULL))
	test_with(list(integer(0), NULL, NULL))
	test_with(list(NULL, integer(0), NULL))
	test_with(list(NULL, NULL, integer(0)))
	test_with(list(NULL, integer(0), integer(0)))
	test_with(list(integer(0), integer(0), integer(0)))
        test_with(list(c(2:5, 7:10), NULL, NULL))
        test_with(list(NULL, c(11:12, 14), NULL))
        test_with(list(NULL, NULL, 1:2))
        test_with(list(7:10, c(11:12, 14), c(1:2, 5)))
    }

    chunkdims <- list(0,             # no chunking (i.e. contiguous data)
                      c(10, 15, 6),
                      c(10, 15, 1),
                      c(10,  1, 6),
                      c( 1, 15, 6),
                      c(10,  1, 1),
                      c( 1, 15, 1),
                      c( 1,  1, 6),
                      c( 4,  5, 5),
                      c( 2,  2, 2))

    ## with an integer array

    a0 <- array(1:900, DIM)
    for (chunkdim in chunkdims) {
        A0 <- writeHDF5Array(a0, chunkdim=chunkdim)
        do_3D_tests(a0, A0, method=1L)
        do_3D_tests(a0, A0, noreduce=TRUE, method=1L)
        do_3D_tests(a0, A0, method=2L)
        do_3D_tests(a0, A0, method=3L)
        do_3D_tests(a0, A0, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a0, A0, method=4L)
            do_3D_tests(a0, A0, method=6L)
            do_3D_tests(a0, A0, method=7L)
            do_3D_sparse_tests(A0)
        }
        do_3D_tests(a0, A0)
    }

    ## with a logical array

    a1 <- a0 %% 3L == 0L
    for (chunkdim in chunkdims) {
        A1 <- writeHDF5Array(a1, chunkdim=chunkdim)
        do_3D_tests(a1, A1, method=1L)
        do_3D_tests(a1, A1, noreduce=TRUE, method=1L)
        do_3D_tests(a1, A1, method=2L)
        do_3D_tests(a1, A1, method=3L)
        do_3D_tests(a1, A1, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a1, A1, method=4L)
            do_3D_tests(a1, A1, method=6L)
            do_3D_tests(a1, A1, method=7L)
            do_3D_sparse_tests(A1)
        }
        do_3D_tests(a1, A1)
        storage.mode(a1) <- "integer"
        do_3D_tests(a1, A1, as.integer=TRUE)
        do_3D_tests(a1, A1, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a1, A1, as.integer=TRUE, method=4L)
            do_3D_tests(a1, A1, as.integer=TRUE, method=6L)
            do_3D_tests(a1, A1, as.integer=TRUE, method=7L)
            do_3D_sparse_tests(A1, as.integer=TRUE)
        }
    }

    ## with a numeric array

    a2 <- array(10 * runif(900), DIM)
    for (chunkdim in chunkdims) {
        A2 <- writeHDF5Array(a2, chunkdim=chunkdim)
        do_3D_tests(a2, A2, method=1L)
        do_3D_tests(a2, A2, noreduce=TRUE, method=1L)
        do_3D_tests(a2, A2, method=2L)
        do_3D_tests(a2, A2, method=3L)
        do_3D_tests(a2, A2, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a2, A2, method=4L)
            do_3D_tests(a2, A2, method=6L)
            do_3D_tests(a2, A2, method=7L)
            do_3D_sparse_tests(A2)
        }
        do_3D_tests(a2, A2)
        storage.mode(a2) <- "integer"
        do_3D_tests(a2, A2, as.integer=TRUE)
        do_3D_tests(a2, A2, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a2, A2, as.integer=TRUE, method=4L)
            do_3D_tests(a2, A2, as.integer=TRUE, method=6L)
            do_3D_tests(a2, A2, as.integer=TRUE, method=7L)
            do_3D_sparse_tests(A2, as.integer=TRUE)
        }
    }

    ## with a character array

    a3 <- array(as.character(1:900), DIM)
    for (chunkdim in chunkdims) {
        A3 <- writeHDF5Array(a3, chunkdim=chunkdim)
        do_3D_tests(a3, A3, method=4L)
        if (!identical(chunkdim, 0)) {
            do_3D_sparse_tests(A3)
        }
    }

    ## with a raw array

    a0 <- a0 %% 256L
    a4 <- a0
    storage.mode(a4) <- "raw"
    for (chunkdim in chunkdims) {
        A4 <- writeHDF5Array(a4, chunkdim=chunkdim)
        do_3D_tests(a4, A4, method=1L)
        do_3D_tests(a4, A4, noreduce=TRUE, method=1L)
        do_3D_tests(a4, A4, method=2L)
        do_3D_tests(a4, A4, method=3L)
        do_3D_tests(a4, A4, noreduce=TRUE, method=3L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a4, A4, method=4L)
            do_3D_tests(a4, A4, method=6L)
            do_3D_tests(a4, A4, method=7L)
            do_3D_sparse_tests(A4)
        }
        do_3D_tests(a4, A4)
        do_3D_tests(a0, A4, as.integer=TRUE)
        do_3D_tests(a0, A4, as.integer=TRUE, method=1L)
        if (!identical(chunkdim, 0)) {
            do_3D_tests(a0, A4, as.integer=TRUE, method=4L)
            do_3D_tests(a0, A4, as.integer=TRUE, method=6L)
            do_3D_tests(a0, A4, as.integer=TRUE, method=7L)
            do_3D_sparse_tests(A4)
        }
    }
}

