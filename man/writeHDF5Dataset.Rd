\name{setHDF5DumpFile}

\alias{setHDF5DumpFile}
\alias{getHDF5DumpFile}
\alias{lsHDF5DumpFile}

\alias{setHDF5DumpName}
\alias{getHDF5DumpName}

\alias{class:HDF5DatasetDump}
\alias{HDF5DatasetDump-class}
\alias{HDF5DatasetDump}

\alias{dimnames,HDF5DatasetDump-method}
\alias{write_to_dump,array,HDF5DatasetDump-method}
\alias{close,HDF5DatasetDump-method}

\alias{coerce,HDF5DatasetDump,HDF5Dataset-method}
\alias{writeHDF5Dataset}
\alias{coerce,ANY,HDF5Dataset-method}
\alias{coerce,HDF5DatasetDump,DelayedArray-method}
\alias{coerce,ANY,HDF5Array-method}

\title{Write an array-like object to an HDF5 file}

\description{
  A function for writting an array-like object to an HDF5 file and other
  low-level utilities to control the location of automatically created
  HDF5 datasets.
}

\usage{
## Write an array-like object to an HDF5 file:
writeHDF5Dataset(x, file, name)

## Control the location of automatically created HDF5 datasets (a.k.a. HDF5
## dump management):
setHDF5DumpFile(file=paste0(tempfile(), ".h5"))
getHDF5DumpFile()
lsHDF5DumpFile()

setHDF5DumpName(name)
getHDF5DumpName()
}

\arguments{
  \item{x}{
    The array-like object to write to an HDF5 file e.g. an ordinary array,
    a \link{DelayedArray} object, a \link{HDF5Dataset} object, etc...

    If \code{x} is a \link{DelayedArray} object, \code{writeHDF5Dataset}
    \emph{realizes} it on disk, that is, all the delayed operations carried
    by the object are executed while the object is written to disk.
    See "On-disk realization of a DelayedArray object as an HDF5 dataset"
    section below for more information.

    If \code{x} is an \link{HDF5Dataset} object, \code{writeHDF5Dataset}
    first wraps it in a \link{DelayedArray} object in order to trigger the
    block-processing strategy used by the on-disk realization mechanism.
  }
  \item{file}{
    For \code{writeHDF5Dataset}: The path (as a single character string)
    to the (new or existing) HDF5 file where to write the dataset.

    For \code{setHDF5DumpFile}: The path to the HDF5 file where to write
    automatically created HDF5 datasets.
  }
  \item{name}{
    For \code{writeHDF5Dataset}: The name of the HDF5 dataset to write.

    For \code{setHDF5DumpName}: The name of the \emph{next automatic dataset}
    to be written to the current output HDF5 file. The supplied name is for
    a one-time use only.
  }
}

\details{
  Use \code{\link{setHDF5DumpFile}} and \code{\link{setHDF5DumpName}} to
  control the location of automatically created HDF5 datasets.
}

\value{
  An invisible HDF5Dataset object for \code{writeHDF5Dataset()}.
}

\section{On-disk realization of a DelayedArray object as an HDF5 dataset}{
  When passed an \link{DelayedArray} object, \code{writeHDF5Dataset}
  \emph{realizes} it on disk, that is, all the delayed operations carried
  by the object are executed on-the-fly while the object is written to disk.
  This uses a block-processing strategy so that the full object is not
  realized at once in memory. Instead the object is processed block by block
  i.e. the blocks are realized in memory and written to disk one at a time.

  In other words, \code{writeHDF5Dataset(x, ...)} is semantically equivalent
  to \code{writeHDF5Dataset(as.array(x), ...)}, except that \code{as.array(x)}
  is not called because this would realize the full object at once in memory.

  See \code{?\link{DelayedArray}} for general information about
  \link{DelayedArray} objects.
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
## ---------------------------------------------------------------------
## writeHDF5Dataset()
## ---------------------------------------------------------------------
library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

cov0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")

cov1 <- drop(cov0[ , , 29000001:29000007])

out_file <- tempfile()
writeHDF5Dataset(cov1, out_file, "cov1")
h5ls(out_file)

## ---------------------------------------------------------------------
## AUTOMATIC HDF5 DATASETS AND DUMP MANAGEMENT
## ---------------------------------------------------------------------
getHDF5DumpFile()

## Use setHDF5DumpFile() to change the current output HDF5 file.
## If the specified file exists, then it must be in HDF5 format or
## an error will be raised. If it doesn't exist, then it will be
## created.
#setHDF5DumpFile("path/to/some/HDF5/file")

lsHDF5DumpFile()

a <- array(1:600, c(150, 4))
h5a <- as(a, "HDF5Dataset")
lsHDF5DumpFile()
A <- as(h5a, "HDF5Array")   # DelayedArray object
A

b <- array(runif(6000), c(4, 2, 150))
h5b <- as(b, "HDF5Dataset")
lsHDF5DumpFile()
B <- as(h5b, "HDF5Array")   # DelayedArray object
B

C <- (log(2 * A + 0.88) - 5)^3 * t(drop(B[ , 1, ]))
C
as(C, "HDF5Dataset")        # realize C on disk
lsHDF5DumpFile()

## Matrix multiplication is not delayed:
P <- C \%*\% matrix(runif(20), nrow=4)
lsHDF5DumpFile()
}
\keyword{methods}
