\name{setHDF5OutputFile}

\alias{setHDF5OutputFile}
\alias{getHDF5OutputFile}
\alias{lsHDF5OutputFile}

\alias{setHDF5OutputName}
\alias{getHDF5OutputName}

\title{Control the location of automatically created HDF5 datasets}

\description{
  Utility functions to control the location of automatically created HDF5
  datasets.
}

\usage{
setHDF5OutputFile(file=paste0(tempfile(), ".h5"))
getHDF5OutputFile()
lsHDF5OutputFile()

setHDF5OutputName(name)
getHDF5OutputName()
}

\arguments{
  \item{file}{
    The path to the \emph{current output HDF5 file}, that is, the HDF5 file
    where to write all newly created datasets.
  }
  \item{name}{
    The name of the \emph{next} dataset to be written to the current output
    HDF5 file. This is for a one-time use only.
  }
}

\note{
  \code{lsHDF5OutputFile()} is a just convenience wrapper for
  \code{rhdf5::\link[rhdf5]{h5ls}(getHDF5OutputFile())}.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \link{HDF5Array} objects.

    \item The \code{\link[rhdf5]{h5ls}} function in the \pkg{rhdf5} package,
          on which \code{lsHDF5OutputFile} is based.

  }
}

\examples{
getHDF5OutputFile()

## Use setHDF5OutputFile() to change the current output HDF5 file.
## If the specified file exists, then it must be in HDF5 format or
## an error will be raised. If it doesn't exist, then it will be
## created.
#setHDF5OutputFile("path/to/some/HDF5/file")

lsHDF5OutputFile()

a <- array(1:600, c(150, 4))
h5a <- HDF5Dataset(a)
lsHDF5OutputFile()
A <- HDF5Array(h5a)   # DelayedArray object
A

b <- array(runif(6000), c(4, 2, 150))
h5b <- HDF5Dataset(b)
lsHDF5OutputFile()
B <- HDF5Array(h5b)   # DelayedArray object
B

C <- (log(2 * A + 0.88) - 5)^3 * t(drop(B[ , 1, ]))
C
HDF5Dataset(C)        # realize C on disk
lsHDF5OutputFile()

## Matrix multiplication is not delayed:
P <- C \%*\% matrix(runif(20), nrow=4)
lsHDF5OutputFile()
}
\keyword{methods}
