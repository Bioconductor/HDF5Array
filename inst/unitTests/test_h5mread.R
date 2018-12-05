test_check_selection <- function()
{
    check_selection <- HDF5Array:::check_selection

    ## specifying 'starts' only (no 'counts')

    checkIdentical(integer(0), check_selection(list()))
    checkIdentical(integer(0), check_selection(list(), dim=integer(0)))
    checkIdentical(15L, check_selection(list(NULL), dim=15))
    checkIdentical(5L, check_selection(list(c(6, 5, 2, 2, 10)), dim=15))
    checkIdentical(30L, check_selection(list(c(15:1, 1:15)), dim=15))
    checkIdentical(5L, check_selection(list(c(6, 5, 2, 2, 10)), dim=3e9))
    checkIdentical(1L, check_selection(list(1e18)))
    checkIdentical(1L, check_selection(list(1e18), dim=1e18))

    checkException(check_selection(list(), dim=15))
    checkException(check_selection(list(NULL), dim=3e9))
    checkException(check_selection(list(0)))
    checkException(check_selection(list(NA)))
    checkException(check_selection(list(NA_integer_)))
    checkException(check_selection(list(NA_real_)))
    checkException(check_selection(list(NaN)))
    checkException(check_selection(list(Inf)))
    checkException(check_selection(list(-Inf)))
    checkException(check_selection(list(1e19)))
    checkException(check_selection(list(18), NULL, 15))

    ## specifying 'starts' and 'counts'

    checkIdentical(integer(0), check_selection(list(), list()))
    checkIdentical(15L, check_selection(list(NULL), list(NULL), dim=15))
    checkIdentical(5L, check_selection(list(c(6, 5, 2, 2, 10)),
                                       list(NULL), dim=15))
    checkIdentical(5L, check_selection(list(c(14, 5, 8)),
                                       list(c( 2, 0, 3)), dim=15))

    checkException(check_selection(list(NULL), list()))
    checkException(check_selection(list(), list(NULL)))
    checkException(check_selection(list(NULL), list(3)))
    checkException(check_selection(list(6), list(-1)))
    checkException(check_selection(list(11), list(6), dim=15))
    checkException(check_selection(list(1), list(3e9)))
    checkException(check_selection(list(c(3, 5, 2)), list(1e9, 1e9, 1e9)))
}

test_check_ordered_selection <- function()
{
    check_ordered_selection <- HDF5Array:::check_ordered_selection

    # TODO!
}

test_reduce_selection <- function()
{
    reduce_selection <- HDF5Array:::reduce_selection

    ## specifying 'starts' only (no 'counts')

    current <- reduce_selection(list(c(2:5, 7:10)))
    checkIdentical(list(c(2L, 7L)), current[[1L]])
    checkIdentical(list(c(4L, 4L)), current[[2L]])

    current <- reduce_selection(list(7:10, c(1:2, 5)))
    checkIdentical(list(7L, c(1L, 5L)), current[[1L]])
    checkIdentical(list(4L, c(2L, 1L)), current[[2L]])

    starts <- list(NULL, c(2, 6))
    checkIdentical(NULL, reduce_selection(starts))  # no reduction
    current <- reduce_selection(starts, dim=5:6)
    checkIdentical(list(NULL, c(2L, 6L)), current[[1L]])
    checkIdentical(list(NULL, NULL), current[[2L]])

    starts <- list(NULL, c(4, 6:7), c(2, 5), integer(0), numeric(0))
    current <- reduce_selection(starts)
    checkIdentical(list(NULL, c(4L, 6L), c(2L, 5L), integer(0), integer(0)),
                   current[[1L]])
    checkIdentical(list(NULL, c(1L, 2L), NULL, NULL, NULL),
                   current[[2L]])

    ## specifying 'starts' and 'counts'

    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, NULL)
    checkIdentical(NULL, reduce_selection(starts, counts))  # no reduction

    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, c(3, 4))
    checkIdentical(NULL, reduce_selection(starts, counts))  # no reduction

    starts <- list(NULL, 5)
    counts <- list(NULL, 0)
    current <- reduce_selection(starts, counts)
    checkIdentical(list(NULL, integer(0)), current[[1L]])
    checkIdentical(list(NULL, integer(0)), current[[2L]])

    starts <- list(NULL, c(6, 9))
    counts <- list(NULL, c(2, 0))
    current <- reduce_selection(starts, counts)
    checkIdentical(list(NULL, 6L), current[[1L]])
    checkIdentical(list(NULL, 2L), current[[2L]])

    starts <- list(NULL, c(2, 5, 5, 6, 11, 11))
    counts <- list(NULL, c(3, 0, 1, 4,  0,  5))
    current <- reduce_selection(starts, counts)
    checkIdentical(list(NULL, c(2L, 11L)), current[[1L]])
    checkIdentical(list(NULL, c(8L,  5L)), current[[2L]])

    starts <- list(NULL, c(2, 6, 10), c(5:10, 15, 3e9 + 1:10), c(99, 2e9))
    counts <- list(NULL, c(3, 4, 1), NULL, c(6e8, 5e8))
    current <- reduce_selection(starts, counts)
    checkIdentical(list(NULL, c(2L, 6L), c(5, 15, 3e9 + 1), c(99L, 2e9L)),
                   current[[1L]])
    checkIdentical(list(NULL, c(3L, 5L), c(6L, 1L, 10L), c(6e8L, 5e8L)),
                   current[[2L]])
}

test_map_starts_to_chunks <- function()
{
    map_starts_to_chunks <- HDF5Array:::map_starts_to_chunks

    ## 1 dimension

    current <- map_starts_to_chunks(list(13:22), 85, 1)
    target <- list(list(1:10), list(as.double(12:21)))
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(13:22), 85, 10)
    target <- list(list(c(8L, 10L)), list(c(1, 2)))
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(c(1:15, 49:51)), 85, 10)
    target <- list(list(c(10L, 15L, 17L, 18L)), list(c(0, 1, 4, 5)))
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(2*(10:35)), 85, 10)
    target <- list(list(1L + 5L*(0:5)), list(as.double(1:6)))
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(c(6e9, 6e9 + 1)), 9e9, 3)
    target <- list(list(1:2), list(c(1999999999, 2000000000)))
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(8e9), 9e9, 3)
    target <- list(list(1L), list(2666666666))
    checkIdentical(target, current)

    ## more dimensions

    current <- map_starts_to_chunks(list(NULL, 13:22, NULL),
                                    c(0, 85, 999), c(0, 10, 1))
    target <- list(list(NULL, c(8L, 10L), NULL), list(NULL, c(1, 2), NULL))
    checkIdentical(target, current)

    ## edge cases

    current <- map_starts_to_chunks(list(), integer(0), integer(0))
    target <- list(list(), list())
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(NULL), 0, 0)
    target <- list(list(NULL), list(NULL))
    checkIdentical(target, current)

    checkException(map_starts_to_chunks(list(NULL), 1, 0))

    current <- map_starts_to_chunks(list(NULL), 0, 1)
    target <- list(list(NULL), list(NULL))
    checkIdentical(target, current)

    current <- map_starts_to_chunks(list(integer(0)), 0, 1)
    target <- list(list(integer(0)), list(numeric(0)))
    checkIdentical(target, current)
}

test_h5mread_2D <- function()
{
    do_2D_tests <- function(m, M, noreduce=FALSE, as.integer=FALSE, method=0L) {
        read <- function(starts, counts=NULL)
            h5mread(M@seed@filepath, M@seed@name,
                    starts=starts, counts=counts,
                    noreduce=noreduce, as.integer=as.integer, method=method)

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

    chunkdims <- list(c(10, 6),
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
        do_2D_tests(m0, M0, method=4L)
        do_2D_tests(m0, M0, method=6L)
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
        do_2D_tests(m1, M1, method=4L)
        do_2D_tests(m1, M1, method=6L)
        do_2D_tests(m1, M1)
        storage.mode(m1) <- "integer"
        do_2D_tests(m1, M1, as.integer=TRUE)
        do_2D_tests(m1, M1, as.integer=TRUE, method=1L)
        do_2D_tests(m1, M1, as.integer=TRUE, method=6L)
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
        do_2D_tests(m2, M2, method=4L)
        do_2D_tests(m2, M2, method=6L)
        do_2D_tests(m2, M2)
        storage.mode(m2) <- "integer"
        do_2D_tests(m2, M2, as.integer=TRUE)
        do_2D_tests(m2, M2, as.integer=TRUE, method=1L)
        do_2D_tests(m2, M2, as.integer=TRUE, method=6L)
    }

    ## with a character matrix

    m3 <- matrix(as.character(1:60), ncol=6)
    for (chunkdim in chunkdims) {
        M3 <- writeHDF5Array(m3, chunkdim=chunkdim)
        do_2D_tests(m3, M3, method=4L)
    }
}

test_h5mread_3D <- function()
{
    DIM <- c(10, 15, 6)

    do_3D_tests <- function(a, A, noreduce=FALSE, as.integer=FALSE, method=0L) {
        read <- function(starts, counts=NULL)
            h5mread(A@seed@filepath, A@seed@name,
                    starts=starts, counts=counts,
                    noreduce=noreduce, as.integer=as.integer, method=method)

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

    chunkdims <- list(c(10, 15, 6),
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
        do_3D_tests(a0, A0, method=4L)
        do_3D_tests(a0, A0, method=6L)
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
        do_3D_tests(a1, A1, method=4L)
        do_3D_tests(a1, A1, method=6L)
        do_3D_tests(a1, A1)
        storage.mode(a1) <- "integer"
        do_3D_tests(a1, A1, as.integer=TRUE)
        do_3D_tests(a1, A1, as.integer=TRUE, method=1L)
        do_3D_tests(a1, A1, as.integer=TRUE, method=6L)
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
        do_3D_tests(a2, A2, method=4L)
        do_3D_tests(a2, A2, method=6L)
        do_3D_tests(a2, A2)
        storage.mode(a2) <- "integer"
        do_3D_tests(a2, A2, as.integer=TRUE)
        do_3D_tests(a2, A2, as.integer=TRUE, method=1L)
        do_3D_tests(a2, A2, as.integer=TRUE, method=6L)
    }

    ## with a character array

    a3 <- array(as.character(1:900), DIM)
    for (chunkdim in chunkdims) {
        A3 <- writeHDF5Array(a3, chunkdim=chunkdim)
        do_3D_tests(a3, A3, method=4L)
    }
}

