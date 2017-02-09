\name{writeHDF5Array}

\alias{setHDF5DumpFile}
\alias{getHDF5DumpFile}
\alias{lsHDF5DumpFile}

\alias{setHDF5DumpName}
\alias{getHDF5DumpName}

\alias{class:HDF5RealizationSink}
\alias{HDF5RealizationSink-class}
\alias{HDF5RealizationSink}

\alias{dimnames,HDF5RealizationSink-method}
\alias{write_to_sink,array,HDF5RealizationSink-method}

\alias{coerce,HDF5RealizationSink,HDF5ArraySeed-method}
\alias{coerce,HDF5RealizationSink,HDF5Array-method}
\alias{coerce,HDF5RealizationSink,DelayedArray-method}

\alias{writeHDF5Array}
\alias{writeHDF5Dataset}

\alias{coerce,ANY,HDF5Array-method}
\alias{coerce,DelayedArray,HDF5Array-method}

\title{Write an array-like object to an HDF5 file}

\description{
  A function for writting an array-like object to an HDF5 file and other
  low-level utilities to control the location of automatically created
  HDF5 datasets.
}

\usage{
## Write an array-like object to an HDF5 file:
writeHDF5Array(x, file=NULL, name=NULL, verbose=FALSE)

## Control the location of automatically created HDF5 datasets (a.k.a.
## HDF5 dump management):
setHDF5DumpFile(file=paste0(tempfile(), ".h5"))
getHDF5DumpFile()
lsHDF5DumpFile()

setHDF5DumpName(name)
getHDF5DumpName()
}

\arguments{
  \item{x}{
    The array-like object to write to an HDF5 file.

    If \code{x} is a \link{DelayedArray} object, \code{writeHDF5Array}
    \emph{realizes} it on disk, that is, all the delayed operations carried
    by the object are executed while the object is written to disk.
    See "On-disk realization of a DelayedArray object as an HDF5 dataset"
    section below for more information.
  }
  \item{file}{
    For \code{writeHDF5Array}: \code{NULL} or the path (as a single string)
    to the (new or existing) HDF5 file where to write the dataset.
    If \code{NULL} then the dataset is written to the \emph{HDF5 dump
    file} (see below) i.e. the path returned by \code{getHDF5DumpFile()}
    is used.

    For \code{setHDF5DumpFile}: The path (as a single string) to the
    \emph{HDF5 dump file}, that is, to the (new or existing) HDF5 file
    where to write automatically created HDF5 datasets.
  }
  \item{name}{
    For \code{writeHDF5Array}: \code{NULL} or the name of the HDF5 dataset
    to write. If \code{NULL} then the name returned by \code{getHDF5DumpName()}
    is used.

    For \code{setHDF5DumpName}: The name of the \emph{next automatic dataset}
    to be written to the current \emph{HDF5 dump file}. The supplied name is
    for a one-time use only.
  }
  \item{verbose}{
    Set to \code{TRUE} to make the function display progress.
  }
}

\details{
  Please note that, depending on the size of the data to write to disk
  and the performance of the disk, \code{writeHDF5Array} can take a long
  time to complete. Use \code{verbose=FALSE} to see its progress.

  Use \code{\link{setHDF5DumpFile}} and \code{\link{setHDF5DumpName}} to
  control the location of automatically created HDF5 datasets.
}

\value{
  A HDF5Array object pointing to the newly written HDF5 dataset on disk
  for \code{writeHDF5Array()}.
}

\section{On-disk realization of a DelayedArray object as an HDF5 dataset}{
  When passed an \link{DelayedArray} object, \code{writeHDF5Array}
  \emph{realizes} it on disk, that is, all the delayed operations carried
  by the object are executed on-the-fly while the object is written to disk.
  This uses a block-processing strategy so that the full object is not
  realized at once in memory. Instead the object is processed block by block
  i.e. the blocks are realized in memory and written to disk one at a time.

  In other words, \code{writeHDF5Array(x, ...)} is semantically equivalent
  to \code{writeHDF5Array(as.array(x), ...)}, except that \code{as.array(x)}
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

    \item \code{\link[SummarizedExperiment]{saveHDF5SummarizedExperiment}} and
          \code{\link[SummarizedExperiment]{loadHDF5SummarizedExperiment}} in
          the \pkg{SummarizedExperiment} package for saving/loading a
          HDF5-based SummarizedExperiment object to/from disk.

    \item The \code{\link[rhdf5]{h5ls}} function in the \pkg{rhdf5} package,
          on which \code{lsHDF5DumpFile} is based.

  }
}

\examples{
## ---------------------------------------------------------------------
## writeHDF5Array()
## ---------------------------------------------------------------------
library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

cov0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")

cov1 <- drop(cov0[ , , 29000001:29000007])

out_file <- tempfile()
writeHDF5Array(cov1, out_file, "cov1")
h5ls(out_file)

## ---------------------------------------------------------------------
## AUTOMATIC HDF5 DATASETS AND DUMP MANAGEMENT
## ---------------------------------------------------------------------
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

## Matrix multiplication is not delayed. It will realize the result and
## thus one needs to set a "realize() backend" first if none is set yet.
setRealizeBackend("HDF5Array")

P <- C \%*\% matrix(runif(20), nrow=4)
lsHDF5DumpFile()
}
\keyword{methods}
