\name{HDF5-dump-management}

\alias{HDF5-dump-management}
\alias{dump-management}

\alias{setHDF5DumpFile}
\alias{getHDF5DumpFile}
\alias{lsHDF5DumpFile}

\alias{setHDF5DumpName}
\alias{getHDF5DumpName}

\alias{showHDF5DumpLog}

\title{HDF5 dump management}

\description{
  A set of utilities to control the location of automatically created HDF5
  datasets.
}

\usage{
setHDF5DumpFile(file)
getHDF5DumpFile()
lsHDF5DumpFile()

setHDF5DumpName(name)
getHDF5DumpName()

showHDF5DumpLog()
}

\arguments{
  \item{file}{
    The path (as a single string) to the \emph{HDF5 dump file}, that is,
    to the (new or existing) HDF5 file where to write automatically created
    HDF5 datasets. If \code{file} is missing, then a new file with an
    automatic name is used for each automatically created HDF5 dataset.
  }
  \item{name}{
    The name of the \emph{next automatic HDF5 dataset} to be written to the
    current \emph{HDF5 dump file}. The supplied name is for a one-time use
    only.
  }
}

\value{
  \code{getHDF5DumpFile} returns the path to the HDF5 file where to write the
  \emph{next automatic HDF5 dataset}.

  \code{getHDF5DumpName} returns the name of the \emph{next automatic HDF5
  dataset}.
}

\note{
  \code{lsHDF5DumpFile()} is a just convenience wrapper for
  \code{rhdf5::\link[rhdf5]{h5ls}(getHDF5DumpFile())}.
}

\seealso{
  \itemize{
    \item \code{\link{writeHDF5Array}} for writting an array-like object
          to an HDF5 file.

    \item \link{HDF5Array} objects.

    \item The \code{\link[rhdf5]{h5ls}} function in the \pkg{rhdf5} package,
          on which \code{lsHDF5DumpFile} is based.

  }
}

\examples{
getHDF5DumpFile()

## Use setHDF5DumpFile() to change the current HDF5 dump file.
## If the specified file exists, then it must be in HDF5 format or
## an error will be raised. If it doesn't exist, then it will be
## created.
#setHDF5DumpFile("path/to/some/HDF5/file")

lsHDF5DumpFile()

a <- array(1:600, c(150, 4))
A <- as(a, "HDF5Array")
lsHDF5DumpFile()
A

b <- array(runif(6000), c(4, 2, 150))
B <- as(b, "HDF5Array")
lsHDF5DumpFile()
B

C <- (log(2 * A + 0.88) - 5)^3 * t(drop(B[ , 1, ]))
as(C, "HDF5Array")  # realize C on disk
lsHDF5DumpFile()

## Matrix multiplication is not delayed: the output matrix is realized
## block by block. The current "realization backend" controls where
## realization happens e.g. in memory if set to NULL or in an HDF5 file
## if set to "HDF5Array". See '?realize' in the DelayedArray package for
## more information about "realization backends".
setRealizationBackend("HDF5Array")
m <- matrix(runif(20), nrow=4)
P <- C \%*\% m
lsHDF5DumpFile()

## See all the HDF5 datasets created in the current session so far:
showHDF5DumpLog()
}
\keyword{methods}
