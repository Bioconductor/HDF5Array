\name{h5mread}

\alias{class:DSetHandle}
\alias{DSetHandle-class}
\alias{DSetHandle}
\alias{destroy_DSetHandle}
\alias{show,DSetHandle-method}
\alias{get_h5mread_returned_type}

\alias{h5mread}

\title{An alternative to \code{rhdf5::h5read}}

\description{
  \code{h5mread} is the result of experimenting with alternative
  \code{rhdf5::h5read} implementations.

  It should still be considered experimental!
}

\usage{
h5mread(filepath, name, starts=NULL, counts=NULL, noreduce=FALSE,
        as.integer=FALSE, method=0L)

get_h5mread_returned_type(filepath, name, as.integer=FALSE)
}

\arguments{
  \item{filepath}{
    The path (as a single character string) to the HDF5 file where the
    dataset to read from is located.
  }
  \item{name}{
    The name of the dataset in the HDF5 file.
  }
  \item{starts, counts}{
    \code{starts} and \code{counts} are used to specify the \emph{array
    selection}. Each argument can be either \code{NULL} or a list with
    one list element per dimension in the dataset.

    If \code{starts} and \code{counts} are both \code{NULL}, then the entire
    dataset is read.

    If \code{starts} is a list, each list element in it must be a vector of
    valid positive indices along the corresponding dimension in the dataset.
    An empty vector (\code{integer(0)}) is accepted and indicates an empty
    selection along that dimension. A \code{NULL} is accepted and indicates
    a \emph{full} selection along the dimension so has the same meaning
    as a missing subscript when subsetting an array-like object with \code{[}.
    (Note that for \code{[} a \code{NULL} subscript indicates an empty
    selection.)

    Each list element in \code{counts} must be \code{NULL} or a vector
    of non-negative integers of the same length as the corresponding
    list element in \code{starts}. Each value in the vector indicates how
    many positions to select starting from the associated start value.
    A \code{NULL} indicates that a single position is selected for each
    value along the corresponding dimension.

    If \code{counts} is \code{NULL}, then each index in each \code{starts}
    list element indicates a single position selection along the corresponding
    dimension. Note that in this case the \code{starts} argument is
    equivalent to the \code{index} argument of \code{\link[rhdf5]{h5read}}
    and \code{\link[DelayedArray]{extract_array}} (with the caveat that
    \code{\link[rhdf5]{h5read}} doesn't accept empty selections).

    Finally note that when \code{counts} is not \code{NULL} then the
    selection described by \code{starts} and \code{counts} must be
    \emph{strictly ascending} along each dimension.
  }
  \item{noreduce}{
    TODO
  }
  \item{as.integer}{
    TODO
  }
  \item{method}{
    TODO
  }
}

\details{
  COMING SOON...
}

\value{
  An array for \code{h5mread}.

  The type of the array that will be returned by \code{h5mread} for
  \code{get_h5mread_returned_type}.
  Equivalent to:
  \preformatted{  typeof(h5mread(filepath, name, rep(list(integer(0)), ndim)))
  } where \code{ndim} is the number of dimensions (a.k.a. the \emph{rank}
  in HDF5 jargon) of the dataset. \code{get_h5mread_returned_type} is
  provided for convenience.
}

\seealso{
  \itemize{
    \item \code{\link[rhdf5]{h5read}} in the \pkg{rhdf5} package.

    \item \code{\link[DelayedArray]{type}} in the \pkg{DelayedArray}
          package.

    \item \code{\link[DelayedArray]{extract_array}} in the \pkg{DelayedArray}
          package.

    \item The \code{\link[TENxBrainData]{TENxBrainData}} dataset (in the
          \pkg{TENxBrainData} package).

    \item \code{\link{h5mread_from_reshaped}} to read data from a virtually
          reshaped HDF5 dataset.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC USAGE
## ---------------------------------------------------------------------
m0 <- matrix((runif(600) - 0.5) * 10, ncol=12)
M0 <- writeHDF5Array(m0, name="M0")

m <- h5mread(path(M0), "M0")
stopifnot(identical(m0, m))

m <- h5mread(path(M0), "M0", starts=list(NULL, c(3, 12:8)))
stopifnot(identical(m0[ , c(3, 12:8)], m))

m <- h5mread(path(M0), "M0", starts=list(integer(0), c(3, 12:8)))
stopifnot(identical(m0[NULL , c(3, 12:8)], m))

m <- h5mread(path(M0), "M0", starts=list(1:5, NULL), as.integer=TRUE)
storage.mode(m0) <- "integer"
stopifnot(identical(m0[1:5, ], m))

a0 <- array(1:350, c(10, 5, 7))
A0 <- writeHDF5Array(a0, filepath=path(M0), name="A0")
h5ls(path(A0))

a <- h5mread(path(A0), "A0", starts=list(c(2, 7), NULL, 6),
                             counts=list(c(4, 2), NULL, NULL))
stopifnot(identical(a0[c(2:5, 7:8), , 6, drop=FALSE], a))

## ---------------------------------------------------------------------
## PERFORMANCE
## ---------------------------------------------------------------------
library(ExperimentHub)
hub <- ExperimentHub()

## With the "sparse" TENxBrainData dataset
## ---------------------------------------
fname0 <- hub[["EH1039"]]
h5ls(fname0)  # all datasets are 1D datasets

index <- list(77 * sample(34088679, 5000, replace=TRUE))
## h5mread() about 3x faster than h5read():
system.time(a <- h5mread(fname0, "mm10/data", index))
system.time(b <- h5read(fname0, "mm10/data", index=index))
stopifnot(identical(a, b))

index <- list(sample(1306127, 7500, replace=TRUE))
## h5mread() about 20x faster than h5read():
system.time(a <- h5mread(fname0, "mm10/barcodes", index))
system.time(b <- h5read(fname0, "mm10/barcodes", index=index))
stopifnot(identical(a, b))

## With the "dense" TENxBrainData dataset
## ---------------------------------------
fname1 <- hub[["EH1040"]]
h5ls(fname1)  # "counts" is a 2D dataset

index <- list(sample(  27998, 250, replace=TRUE),
              sample(1306127, 250, replace=TRUE))
## h5mread() about 2x faster than h5read():
system.time(a <- h5mread(fname1, "counts", index))
system.time(b <- h5read(fname1, "counts", index=index))
stopifnot(identical(a, b))

## The bigger the selection, the greater the speedup between
## h5read() and h5mread():
\dontrun{
  index <- list(sample(  27998, 1000, replace=TRUE),
                sample(1306127, 1000, replace=TRUE))
  ## h5mread() about 8x faster than h5read() (22s vs 3min):
  system.time(a <- h5mread(fname1, "counts", index))
  system.time(b <- h5read(fname1, "counts", index=index))
  stopifnot(identical(a, b))
}
}
\keyword{methods}
