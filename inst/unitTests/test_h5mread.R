test_reduce_selection <- function()
{
    ## specifying 'starts' only (no 'counts')

    starts <- list(NULL, c(2, 6))

    checkIdentical(NULL, reduce_selection(starts))

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
    checkIdentical(NULL, reduce_selection(starts, counts))

    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, c(3, 4))
    checkIdentical(NULL, reduce_selection(starts, counts))

    starts <- list(NULL, c(2, 6, 10), c(5:10, 15, 3e9 + 1:10), c(99, 2e9))
    counts <- list(NULL, c(3, 4, 1), NULL, c(6e8, 5e8))
    current <- reduce_selection(starts, counts)
    checkIdentical(list(NULL, c(2L, 6L), c(5, 15, 3e9 + 1), c(99L, 2e9L)),
                   current[[1L]])
    checkIdentical(list(NULL, c(3L, 5L), c(6L, 1L, 10L), c(6e8L, 5e8L)),
                   current[[2L]])
}

test_h5mread <- function()
{
    do_tests <- function(m, filepath, name, as.integer=FALSE) {
        read <- function(starts, counts=NULL)
            h5mread(filepath, name, starts=starts, counts=counts,
                    as.integer=as.integer)

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

    ## with an array of doubles

    m1 <- matrix(10 * runif(60), ncol=6)
    M1 <- as(m1, "HDF5Array")
    do_tests(m1, M1@seed@filepath, M1@seed@name)

    storage.mode(m1) <- "integer"
    do_tests(m1, M1@seed@filepath, M1@seed@name, as.integer=TRUE)
}

