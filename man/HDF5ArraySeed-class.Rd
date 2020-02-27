\name{HDF5ArraySeed-class}
\docType{class}

\alias{class:HDF5ArraySeed}
\alias{HDF5ArraySeed-class}
\alias{HDF5ArraySeed}

\alias{path,HDF5ArraySeed-method}
\alias{path<-,HDF5ArraySeed-method}
\alias{type,HDF5ArraySeed-method}
\alias{dim,HDF5ArraySeed-method}
\alias{extract_array,HDF5ArraySeed-method}
\alias{chunkdim,HDF5ArraySeed-method}

\title{HDF5ArraySeed objects}

\description{
  HDF5ArraySeed is a low-level helper class for representing a pointer
  to an HDF5 dataset.
  HDF5ArraySeed objects are not intended to be used directly. Most
  end users should create and manipulate \link{HDF5Array} objects
  instead. See \code{?\link{HDF5Array}} for more information.
}

\usage{
## Constructor function:
HDF5ArraySeed(filepath, name, type=NA)
}

\arguments{
  \item{filepath, name, type}{
    See \code{?\link{HDF5Array}} for a description of these arguments.
  }
}

\details{
  No operation can be performed directly on an HDF5ArraySeed object.
  It first needs to be wrapped in a \link[DelayedArray]{DelayedArray}
  object. The result of this wrapping is an \link{HDF5Array} object
  (an \link{HDF5Array} object is just an HDF5ArraySeed object wrapped
  in a \link[DelayedArray]{DelayedArray} object).
}

\value{
  An HDF5ArraySeed object.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
library(h5vcData)
tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")

library(rhdf5)  # for h5ls()
h5ls(tally_file)

seed <- HDF5ArraySeed(tally_file, "/ExampleStudy/16/Coverages")
seed
path(seed)
dim(seed)
chunkdim(seed)
}
\keyword{classes}
\keyword{methods}
