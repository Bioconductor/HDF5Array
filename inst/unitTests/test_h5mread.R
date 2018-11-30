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

test_h5mread <- function()
{
    do_tests <- function(m, filepath, name, as.integer=FALSE, method=0L) {
        read <- function(starts, counts=NULL)
            h5mread(filepath, name, starts=starts, counts=counts,
                    as.integer=as.integer, method=method)

        current <- read(list(NULL, NULL))
        checkIdentical(m, current)

        current <- read(list(c(2:5, 7:10), NULL))
        checkIdentical(m[c(2:5, 7:10), , drop=FALSE], current)

        current <- read(list(NULL, 1:2))
        checkIdentical(m[ , 1:2, drop=FALSE], current)

        current <- read(list(7:10, c(1:2, 5)))
        checkIdentical(m[7:10, c(1:2, 5), drop=FALSE], current)

        current <- read(list(integer(0), integer(0)))
        checkIdentical(m[integer(0), integer(0), drop=FALSE], current)

        if (method > 1L)
            return()

        starts <- list(integer(0), 4L)
        counts <- list(integer(0), 2L)
        current <- read(starts, counts)
        checkIdentical(m[integer(0), 4:5, drop=FALSE], current)

        starts <- list(5L, integer(0))
        counts <- list(4L, integer(0))
        current <- read(starts, counts)
        checkIdentical(m[5:8, integer(0), drop=FALSE], current)

        starts <- list(c(2L, 5L), 4L)
        counts <- list(c(3L, 4L), 2L)
        current <- read(starts, counts)
        checkIdentical(m[2:8, 4:5, drop=FALSE], current)
    }

    ## with an array of integers

    m0 <- matrix(1:60, ncol=6)
    M0 <- as(m0, "HDF5Array")
    do_tests(m0, M0@seed@filepath, M0@seed@name)
    do_tests(m0, M0@seed@filepath, M0@seed@name, method=1L)
    do_tests(m0, M0@seed@filepath, M0@seed@name, method=6L)

    ## with an array of doubles

    m1 <- matrix(10 * runif(60), ncol=6)
    M1 <- as(m1, "HDF5Array")
    do_tests(m1, M1@seed@filepath, M1@seed@name)
    do_tests(m1, M1@seed@filepath, M1@seed@name, method=1L)
    do_tests(m1, M1@seed@filepath, M1@seed@name, method=6L)

    storage.mode(m1) <- "integer"
    do_tests(m1, M1@seed@filepath, M1@seed@name, as.integer=TRUE)
    do_tests(m1, M1@seed@filepath, M1@seed@name, as.integer=TRUE, method=1L)
    do_tests(m1, M1@seed@filepath, M1@seed@name, as.integer=TRUE, method=6L)
}

