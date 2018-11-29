\name{h5mread}

\alias{h5mread}
\alias{reduce_selection}

\title{An alternative to \code{rhdf5::h5read}}

\description{
  \code{h5mread} is the result of experimenting with alternative
  \code{rhdf5::h5read} implementations.

  It should still be considered experimental!
}

\usage{
h5mread(filepath, name, starts, counts=NULL, noreduce=FALSE,
        as.integer=FALSE, method=0L)
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
    2 lists of the same length specifying the array selection.
    If \code{counts} is NULL, then all the counts are considered to be 1.
  }
  \item{noreduce}{
  }
  \item{as.integer}{
  }
  \item{method}{
  }
}

\value{
  An array.
}

\seealso{
  \itemize{
    \item \code{\link[rhdf5]{h5read}} in the \pkg{rhdf5}.
  }
}

\examples{
}
\keyword{methods}
