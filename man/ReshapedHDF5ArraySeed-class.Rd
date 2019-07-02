\name{ReshapedHDF5ArraySeed-class}
\docType{class}

\alias{class:ReshapedHDF5ArraySeed}
\alias{ReshapedHDF5ArraySeed-class}
\alias{ReshapedHDF5ArraySeed}

\alias{dim,ReshapedHDF5ArraySeed-method}
\alias{extract_array,ReshapedHDF5ArraySeed-method}
\alias{chunkdim,ReshapedHDF5ArraySeed-method}

\title{ReshapedHDF5ArraySeed objects}

\description{
  ReshapedHDF5ArraySeed is a low-level helper class for representing a
  pointer to a virtually reshaped HDF5 dataset.

  ReshapedHDF5ArraySeed objects are not intended to be used directly.
  Most end users should create and manipulate \link{ReshapedHDF5Array}
  objects instead. See \code{?\link{ReshapedHDF5Array}} for more information.
}

\usage{
## Constructor function:
ReshapedHDF5ArraySeed(filepath, name, dim, type=NA)
}

\arguments{
  \item{filepath, name, dim, type}{
    See \code{?\link{ReshapedHDF5Array}} for a description of these arguments.
  }
}

\details{
  No operation can be performed directly on a ReshapedHDF5ArraySeed object.
  It first needs to be wrapped in a \link[DelayedArray]{DelayedArray}
  object. The result of this wrapping is a \link{ReshapedHDF5Array} object
  (a \link{ReshapedHDF5Array} object is just a ReshapedHDF5ArraySeed object
  wrapped in a \link[DelayedArray]{DelayedArray} object).
}

\value{
  A ReshapedHDF5ArraySeed object.
}

\seealso{
  \itemize{
    \item \link{ReshapedHDF5Array} objects.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
library(h5vcData)
tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")

library(rhdf5)  # for h5ls()
h5ls(tally_file)

## Collapse the first 2 dimensions:
seed <- ReshapedHDF5ArraySeed(tally_file, "/ExampleStudy/16/Coverages",
                              dim=c(12, 90354753))
seed
path(seed)
dim(seed)
chunkdim(seed)
}
\keyword{classes}
\keyword{methods}
