\name{cbind-methods}

\alias{cbind-methods}

\alias{dim,RowBinder-method}
\alias{dim,ColBinder-method}
\alias{dimnames,RowBinder-method}
\alias{dimnames,ColBinder-method}

\alias{rbind}
\alias{rbind,DelayedMatrix-method}
\alias{rbind,DelayedArray-method}

\alias{cbind}
\alias{cbind,DelayedMatrix-method}
\alias{cbind,DelayedArray-method}

\title{Bind DelayedArray objects along their rows or columns}

\description{
  Methods for binding DelayedArray objects along their rows or columns.
}

\details{
  \code{rbind} and \code{cbind} methods are defined for \link{DelayedArray}
  objects. They perform delayed binding along the rows (\code{rbind}) or
  columns (\code{cbind}) of the objects passed to them.
}

\seealso{
  \itemize{
    \item \code{\link[base]{cbind}} in the \pkg{base} package for the
          corresponding operations.

    \item \link{DelayedArray-utils} for common operations on
          \link{DelayedArray} objects.

    \item \link{DelayedArray} objects.

    \item \link{HDF5Array} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
library(rhdf5)
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

M1 <- HDF5Array(toy_h5, "M1")
M2 <- HDF5Array(toy_h5, "M2")

M <- rbind(M1, t(M2))
M
colMeans(M)
}
\keyword{methods}
