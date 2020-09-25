test_check_uaselection <- function()
{
    check_uaselection <- HDF5Array:::check_uaselection

    ## specifying 'starts' only (no 'counts')

    checkIdentical(6:4, check_uaselection(6:4))
    checkIdentical(c(6L, 1L, 2L),
                   check_uaselection(6:4, starts=list(NULL, 5, 3:4)))

    checkIdentical(integer(0), check_uaselection(integer(0)))
    checkIdentical(integer(0), check_uaselection(integer(0), starts=NULL))
    checkIdentical(integer(0), check_uaselection(integer(0), starts=list()))

    checkIdentical(15L, check_uaselection(15))
    checkIdentical(15L, check_uaselection(15, starts=NULL))
    checkIdentical(15L, check_uaselection(15, starts=list(NULL)))

    checkIdentical(5L, check_uaselection(15, starts=list(c(6, 5, 2, 2, 10))))
    checkIdentical(30L, check_uaselection(15, starts=list(c(15:1, 1:15))))
    checkIdentical(5L, check_uaselection(3e9, starts=list(c(6, 5, 2, 2, 10))))
    checkIdentical(1L, check_uaselection(-1, starts=list(1e18)))
    checkIdentical(1L, check_uaselection(1e18, starts=list(1e18)))

    checkException(check_uaselection(15, starts=list()))
    checkException(check_uaselection(3e9))
    checkException(check_uaselection(-1, starts=list(0)))
    checkException(check_uaselection(-1, starts=list(NA)))
    checkException(check_uaselection(-1, starts=list(NA_integer_)))
    checkException(check_uaselection(-1, starts=list(NA_real_)))
    checkException(check_uaselection(-1, starts=list(NaN)))
    checkException(check_uaselection(-1, starts=list(Inf)))
    checkException(check_uaselection(-1, starts=list(-Inf)))
    checkException(check_uaselection(-1, starts=list(1e19)))
    checkException(check_uaselection(15, starts=list(18)))

    ## specifying 'starts' and 'counts'

    checkIdentical(integer(0), check_uaselection(integer(0), list(), list()))
    checkIdentical(15L, check_uaselection(15, starts=list(NULL),
                                              counts=list(NULL)))
    checkIdentical(5L, check_uaselection(15, starts=list(c(6, 5, 2, 2, 10)),
                                             counts=list(NULL)))
    checkIdentical(5L, check_uaselection(15, starts=list(c(14, 5, 8)),
                                             counts=list(c( 2, 0, 3))))

    checkException(check_uaselection(-1, starts=list(NULL),
                                         counts=list()))
    checkException(check_uaselection(-1, starts=list(),
                                         counts=list(NULL)))
    checkException(check_uaselection(-1, starts=list(NULL),
                                         counts=list(3)))
    checkException(check_uaselection(-1, starts=list(6),
                                         counts=list(-1)))
    checkException(check_uaselection(15, starts=list(11),
                                         counts=list(6)))
    checkException(check_uaselection(-1, starts=list(1),
                                         counts=list(3e9)))
    checkException(check_uaselection(-1, starts=list(c(3, 5, 2)),
                                         counts=list(1e9, 1e9, 1e9)))
}

test_check_ordered_uaselection <- function()
{
    check_ordered_uaselection <- HDF5Array:::check_ordered_uaselection

    # TODO!
}

test_reduce_uaselection <- function()
{
    reduce_uaselection <- HDF5Array:::reduce_uaselection

    ## specifying 'starts' only (no 'counts')

    current <- reduce_uaselection(-1, starts=list(c(2:5, 7:10)))
    checkIdentical(list(c(2L, 7L)), current[[1L]])
    checkIdentical(list(c(4L, 4L)), current[[2L]])

    current <- reduce_uaselection(c(-1, -1), starts=list(7:10, c(1:2, 5)))
    checkIdentical(list(7L, c(1L, 5L)), current[[1L]])
    checkIdentical(list(4L, c(2L, 1L)), current[[2L]])

    starts <- list(NULL, c(2, 6))
    checkIdentical(NULL, reduce_uaselection(c(-1, -1), starts))  # no reduction
    current <- reduce_uaselection(5:6, starts)
    checkIdentical(list(NULL, c(2L, 6L)), current[[1L]])
    checkIdentical(list(NULL, NULL), current[[2L]])

    starts <- list(NULL, c(4, 6:7), c(2, 5), integer(0), numeric(0))
    current <- reduce_uaselection(rep(-1, 5), starts)
    checkIdentical(list(NULL, c(4L, 6L), c(2L, 5L), integer(0), integer(0)),
                   current[[1L]])
    checkIdentical(list(NULL, c(1L, 2L), NULL, NULL, NULL),
                   current[[2L]])

    ## specifying 'starts' and 'counts'

    dim <- c(-1, -1)  # unspecified dimensions
    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, NULL)
    checkIdentical(NULL, reduce_uaselection(dim, starts, counts))  # no reduction

    starts <- list(NULL, c(2, 6))
    counts <- list(NULL, c(3, 4))
    checkIdentical(NULL, reduce_uaselection(dim, starts, counts))  # no reduction

    starts <- list(NULL, 5)
    counts <- list(NULL, 0)
    current <- reduce_uaselection(dim, starts, counts)
    checkIdentical(list(NULL, integer(0)), current[[1L]])
    checkIdentical(list(NULL, integer(0)), current[[2L]])

    starts <- list(NULL, c(6, 9))
    counts <- list(NULL, c(2, 0))
    current <- reduce_uaselection(dim, starts, counts)
    checkIdentical(list(NULL, 6L), current[[1L]])
    checkIdentical(list(NULL, 2L), current[[2L]])

    starts <- list(NULL, c(2, 5, 5, 6, 11, 11))
    counts <- list(NULL, c(3, 0, 1, 4,  0,  5))
    current <- reduce_uaselection(dim, starts, counts)
    checkIdentical(list(NULL, c(2L, 11L)), current[[1L]])
    checkIdentical(list(NULL, c(8L,  5L)), current[[2L]])

    starts <- list(NULL, c(2, 6, 10), c(5:10, 15, 3e9 + 1:10), c(99, 2e9))
    counts <- list(NULL, c(3, 4, 1), NULL, c(6e8, 5e8))
    current <- reduce_uaselection(rep(-1, 4), starts, counts)
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

