test_h5mread_from_reshaped_2D <- function()
{
    m0 <- matrix(1:60, ncol=6)
    M0 <- writeHDF5Array(m0, name="M0")

    ## Invalid/unsupported reshapings.
    checkException(h5mread_from_reshaped(path(M0), "M0", dim=c(11, 6),
                                         starts=list(NULL, NULL)))
    checkException(h5mread_from_reshaped(path(M0), "M0", dim=c(5, 2, 6),
                                         starts=list(NULL, NULL, NULL)))
    checkException(h5mread_from_reshaped(path(M0), "M0", dim=c(6, 10),
                                         starts=list(NULL, NULL)))

    ## No reshaping.
    dim <- dim(m0)
    checkException(h5mread_from_reshaped(path(M0), "M0", dim,
                                         starts=list(NULL)))
    checkException(h5mread_from_reshaped(path(M0), "M0", dim,
                                         starts=list(NULL, NULL, NULL)))
    starts <- list(NULL, NULL)
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m0, current)
    starts <- list(c(3:7, 10:6), NULL)
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m0[starts[[1]], , drop=FALSE], current)

    ## Collapse the 2 dimensions.
    dim <- 60
    checkException(h5mread_from_reshaped(path(M0), "M0", dim,
                                         starts=list(NULL, NULL)))
    checkException(h5mread_from_reshaped(path(M0), "M0", dim,
                                         starts=list("xx")))
    m1 <- `dim<-`(m0, dim)  # reshape 'm0'
    starts <- list(NULL)
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m1, current)
    starts <- list(1:60)
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m1, current)
    starts <- list(c(60:1, 1:60))
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m1[starts[[1]], drop=FALSE], current)
    starts <- list(c(31:29, 8:12))
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m1[starts[[1]], drop=FALSE], current)
    starts <- list(integer(0))
    current <- h5mread_from_reshaped(path(M0), "M0", dim, starts=starts)
    checkIdentical(m1[starts[[1]], drop=FALSE], current)
}

test_h5mread_from_reshaped_3D <- function()
{
    a0 <- array(1:350, c(10, 5, 7))
    A0 <- writeHDF5Array(a0, name="A0")

    ## Invalid/unsupported reshapings.
    checkException(h5mread_from_reshaped(path(A0), "A0", dim=c(10, 4, 7),
                                         starts=list(NULL, NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim=c(10, 5, 7, 1),
                                         starts=list(NULL, NULL, NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim=c(2, 25, 7),
                                         starts=list(NULL, NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim=c(5, 10, 7),
                                         starts=list(NULL, NULL, NULL)))

    ## No reshaping.
    dim <- dim(a0)
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, NULL, NULL, NULL)))
    starts <- list(NULL, NULL, NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a0, current)
    starts <- list(c(3:7, 10:6), NULL, 7:6)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a0[starts[[1]], , starts[[3]], drop=FALSE], current)

    ## Collapse the first 2 dimensions.
    dim <- c(50, 7)
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list("xx", NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, "xx")))
    a1 <- `dim<-`(a0, dim)  # reshape 'a0'
    starts <- list(NULL, NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(1:50, NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(c(50:1, 1:50), NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], , drop=FALSE], current)
    starts <- list(c(31:29, 8:12), 7:5)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], drop=FALSE], current)
    starts <- list(integer(0), NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], , drop=FALSE], current)
    starts <- list(integer(0), 7:5)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], drop=FALSE], current)
    starts <- list(integer(0), integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], drop=FALSE], current)

    ## Collapse the last 2 dimensions.
    dim <- c(10, 35)
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list("xx", NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, "xx")))
    a1 <- `dim<-`(a0, dim)  # reshape 'a0'
    starts <- list(NULL, NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(NULL, 1:35)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(NULL, c(35:1, 1:35))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[ , starts[[2]], drop=FALSE], current)
    starts <- list(7:5, c(31:29, 8:16))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], drop=FALSE], current)
    starts <- list(NULL, integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[ , starts[[1]], drop=FALSE], current)
    starts <- list(7:5, integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], drop=FALSE], current)
    starts <- list(integer(0), integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], drop=FALSE], current)

    ## Collapse the 3 dimensions.
    dim <- 350
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list(NULL, NULL)))
    checkException(h5mread_from_reshaped(path(A0), "A0", dim,
                                         starts=list("xx")))
    a1 <- `dim<-`(a0, dim)  # reshape 'a0'
    starts <- list(NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(1:350)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(c(350:1, 1:350))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], drop=FALSE], current)
    starts <- list(c(71:29, 68:133))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], drop=FALSE], current)
    starts <- list(integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], drop=FALSE], current)
}

test_h5mread_from_reshaped_4D <- function()
{
    a0 <- array(runif(4200), c(15, 4, 7, 10))
    A0 <- writeHDF5Array(a0, name="A0")

    ## Collapse the 2nd & 3rd dimensions.
    dim <- c(15, 28, 10)
    a1 <- `dim<-`(a0, dim)  # reshape 'a0'
    starts <- list(NULL, NULL, NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(NULL, 1:28, NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1, current)
    starts <- list(NULL, c(28:1, 1:28), NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[ , starts[[2]], , drop=FALSE], current)
    starts <- list(NULL, c(21:9, 7:15), 7:5)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[ , starts[[2]], starts[[3]], drop=FALSE], current)
    starts <- list(NULL, integer(0), NULL)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[ , starts[[2]], , drop=FALSE], current)
    starts <- list(8:14, integer(0), 7:5)
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], starts[[3]], drop=FALSE],
                   current)
    starts <- list(14:8, integer(0), integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], starts[[3]], drop=FALSE],
                   current)
    starts <- list(integer(0), integer(0), integer(0))
    current <- h5mread_from_reshaped(path(A0), "A0", dim, starts=starts)
    checkIdentical(a1[starts[[1]], starts[[2]], starts[[3]], drop=FALSE],
                   current)
}

