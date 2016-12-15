\name{setHDF5DumpFile}

\alias{setHDF5DumpFile}
\alias{getHDF5DumpFile}
\alias{lsHDF5DumpFile}

\alias{setHDF5DumpName}
\alias{getHDF5DumpName}

\alias{class:HDF5DatasetDump}
\alias{HDF5DatasetDump-class}
\alias{HDF5DatasetDump}

\alias{write_to_dump,array,HDF5DatasetDump-method}
\alias{close,HDF5DatasetDump-method}

\title{Control the location of automatically created HDF5 datasets}

\description{
  Utility functions to control the location of automatically created HDF5
  datasets.
}

\usage{
setHDF5DumpFile(file=paste0(tempfile(), ".h5"))
getHDF5DumpFile()
lsHDF5DumpFile()

setHDF5DumpName(name)
getHDF5DumpName()
}

\arguments{
  \item{file}{
    The path to the \emph{current output HDF5 file}, that is, to the HDF5
    file where all newly created datasets shall be written.
  }
  \item{name}{
    The name of the \emph{next} dataset to be written to the current output
    HDF5 file. This is for a one-time use only.
  }
}

\note{
  \code{lsHDF5DumpFile()} is a just convenience wrapper for
  \code{rhdf5::\link[rhdf5]{h5ls}(getHDF5DumpFile())}.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \link{HDF5Array} objects.

    \item The \code{\link[rhdf5]{h5ls}} function in the \pkg{rhdf5} package,
          on which \code{lsHDF5DumpFile} is based.

  }
}

\examples{
getHDF5DumpFile()

## Use setHDF5DumpFile() to change the current output HDF5 file.
## If the specified file exists, then it must be in HDF5 format or
## an error will be raised. If it doesn't exist, then it will be
## created.
#setHDF5DumpFile("path/to/some/HDF5/file")

lsHDF5DumpFile()

a <- array(1:600, c(150, 4))
h5a <- HDF5Dataset(a)
lsHDF5DumpFile()
A <- HDF5Array(h5a)   # DelayedArray object
A

b <- array(runif(6000), c(4, 2, 150))
h5b <- HDF5Dataset(b)
lsHDF5DumpFile()
B <- HDF5Array(h5b)   # DelayedArray object
B

C <- (log(2 * A + 0.88) - 5)^3 * t(drop(B[ , 1, ]))
C
HDF5Dataset(C)        # realize C on disk
lsHDF5DumpFile()

## Matrix multiplication is not delayed:
P <- C \%*\% matrix(runif(20), nrow=4)
lsHDF5DumpFile()
}
\keyword{methods}
