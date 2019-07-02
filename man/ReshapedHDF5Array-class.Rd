\name{ReshapedHDF5Array-class}
\docType{class}

\alias{class:ReshapedHDF5Array}
\alias{ReshapedHDF5Array-class}
\alias{ReshapedHDF5Array}

\alias{DelayedArray,ReshapedHDF5ArraySeed-method}

\alias{class:ReshapedHDF5Matrix}
\alias{ReshapedHDF5Matrix-class}
\alias{ReshapedHDF5Matrix}

\alias{matrixClass,ReshapedHDF5Array-method}
\alias{coerce,ReshapedHDF5Array,ReshapedHDF5Matrix-method}
\alias{coerce,ReshapedHDF5Matrix,ReshapedHDF5Array-method}
\alias{coerce,ANY,ReshapedHDF5Matrix-method}

\title{Virtually reshaped HDF5 datasets as DelayedArray objects}

\description{
  The ReshapedHDF5Array class is a \link[DelayedArray]{DelayedArray}
  subclass for representing an HDF5 dataset with a user-supplied upfront
  virtual reshaping.

  All the operations available for \link[DelayedArray]{DelayedArray}
  objects work on ReshapedHDF5Array objects.
}

\usage{
## Constructor function:
ReshapedHDF5Array(filepath, name, dim, type=NA)
}

\arguments{
  \item{filepath, name, type}{
    See \code{?\link{HDF5Array}} for a description of these arguments.
  }
  \item{dim}{
    A vector of dimensions that describes the virtual reshaping i.e.
    the reshaping that is virtually applied upfront to the HDF5 dataset
    when the ReshapedHDF5Array object gets constructed.

    Note that the HDF5 dataset is treated as read-only so is not
    \emph{effectively} reshaped, that is, the dataset dimensions encoded
    in the HDF5 file are not mmodified.

    Also please note that arbitrary reshapings are not supported.
    Only reshapings that reduce the number of dimensions by collapsing
    a group of consecutive dimensions into a single dimension are supported.
    For example, reshaping a 10 x 3 x 5 x 1000 array as a 10 x 15 x 1000
    array or as a 150 x 1000 matrix is supported.
  }
}

\value{
  A ReshapedHDF5Array object.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing HDF5 datasets as
          \link[DelayedArray]{DelayedArray} objects without upfront
          virtual reshaping.

    \item \link[DelayedArray]{DelayedArray} objects in the \pkg{DelayedArray}
          package.

    \item \code{\link{writeHDF5Array}} for writing an array-like object
          to an HDF5 file.

    \item \code{\link{saveHDF5SummarizedExperiment}} and
          \code{\link{loadHDF5SummarizedExperiment}} in this
          package (the \pkg{HDF5Array} package) for saving/loading
          an HDF5-based \link[SummarizedExperiment]{SummarizedExperiment}
          object to/from disk.

    \item The \link{ReshapedHDF5ArraySeed} helper class.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## CONSTRUCTION
## ---------------------------------------------------------------------
library(h5vcData)
tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")

library(rhdf5)  # for h5ls()
h5ls(tally_file)

## Pick up "Coverages" dataset for Human chromosome 16 and collapse its
## first 2 dimensions:
cvg <- ReshapedHDF5Array(tally_file, "/ExampleStudy/16/Coverages",
                         dim=c(12, 90354753))
cvg

is(cvg, "DelayedArray")  # TRUE
seed(cvg)
path(cvg)
dim(cvg)
chunkdim(cvg)
}
\keyword{classes}
\keyword{methods}
