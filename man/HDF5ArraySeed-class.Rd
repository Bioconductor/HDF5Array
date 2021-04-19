\name{HDF5ArraySeed-class}
\docType{class}

\alias{class:HDF5ArraySeed}
\alias{HDF5ArraySeed-class}
\alias{HDF5ArraySeed}

\alias{path,HDF5ArraySeed-method}
\alias{path<-,HDF5ArraySeed-method}
\alias{type,HDF5ArraySeed-method}
\alias{dim,HDF5ArraySeed-method}
\alias{dimnames,HDF5ArraySeed-method}
\alias{extract_array,HDF5ArraySeed-method}
\alias{is_sparse,HDF5ArraySeed-method}
\alias{is_sparse<-,HDF5ArraySeed-method}
\alias{extract_sparse_array,HDF5ArraySeed-method}
\alias{chunkdim,HDF5ArraySeed-method}

\alias{updateObject,HDF5ArraySeed-method}

\title{HDF5ArraySeed objects}

\description{
  HDF5ArraySeed is a low-level helper class for representing a
  pointer to an HDF5 dataset.

  Note that an HDF5ArraySeed object is not intended to be used directly.
  Most end users will typically create and manipulate a higher-level
  \link{HDF5Array} object instead. See \code{?\link{HDF5Array}} for
  more information.
}

\usage{
## --- Constructor function ---

HDF5ArraySeed(filepath, name, as.sparse=FALSE, type=NA)

## --- Accessors --------------

\S4method{path}{HDF5ArraySeed}(object)

\S4method{path}{HDF5ArraySeed}(object) <- value

\S4method{dim}{HDF5ArraySeed}(x)

\S4method{dimnames}{HDF5ArraySeed}(x)

\S4method{type}{HDF5ArraySeed}(x)

\S4method{is_sparse}{HDF5ArraySeed}(x)

\S4method{is_sparse}{HDF5ArraySeed}(x) <- value

\S4method{chunkdim}{HDF5ArraySeed}(x)

## --- Data extraction --------

\S4method{extract_array}{HDF5ArraySeed}(x, index)

\S4method{extract_sparse_array}{HDF5ArraySeed}(x, index)
}

\arguments{
  \item{filepath, name, as.sparse, type}{
    See \code{?\link{HDF5Array}} for a description of these arguments.
  }
  \item{object, x}{
    An HDF5ArraySeed object or derivative.
  }
  \item{value}{
    For the \code{path()} setter: The new path (as a single string) to the
    HDF5 file where the dataset is located.

    For the \code{is_sparse()} setter: \code{TRUE} or \code{FALSE}.
  }
  \item{index}{
    See \code{?\link[DelayedArray]{extract_array}} in the \pkg{DelayedArray}
    package.
  }
}

\details{
  The HDF5ArraySeed class has one direct subclass: \link{Dense_H5ADMatrixSeed}.
  See \code{?\link{Dense_H5ADMatrixSeed}} for more information.

  Note that the implementation of HDF5ArraySeed objects follows the widely
  adopted convention of transposing HDF5 matrices when they get loaded into R.

  Finally note that an HDF5ArraySeed object supports a very limited set
  of methods:
  \itemize{
    \item \code{path()}: Returns the path to the HDF5 file where the dataset
          is located.
    \item \code{dim()}, \code{dimnames()}.
    \item \code{type()}, \code{is_sparse()}, \code{chunkdim()},
          \code{extract_array()}, \code{extract_sparse_array()}: These generics
          are defined and documented in the \pkg{DelayedArrray} package.
  }
}

\value{
  \code{HDF5ArraySeed()} returns an HDF5ArraySeed object.
}

\section{HDF5ArraySeed vs HDF5Array objects}{
  In order to have access to the full set of operations that are available
  for \link[DelayedArray]{DelayedArray} objects, an HDF5ArraySeed object
  first needs to be wrapped in a \link[DelayedArray]{DelayedArray} object,
  typically by calling the \code{\link[DelayedArray]{DelayedArray}()}
  constructor on it.

  This is what the \code{\link{HDF5Array}()} constructor function does.

  Note that the result of this wrapping is an \link{HDF5Array} object,
  which is just an HDF5ArraySeed object wrapped in a
  \link[DelayedArray]{DelayedArray} object.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects.

    \item \code{\link{h5ls}} to list the content of an HDF5 file.
  }
}

\examples{
library(h5vcData)
tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

name <- "/ExampleStudy/16/Coverages"  # name of the dataset of interest
seed1 <- HDF5ArraySeed(tally_file, name)
seed1
path(seed1)
dim(seed1)
chunkdim(seed1)

seed2 <- HDF5ArraySeed(tally_file, name, as.sparse=TRUE)
seed2

## Alternatively:
is_sparse(seed1) <- TRUE
seed1  # same as 'seed2'

DelayedArray(seed1)
stopifnot(class(DelayedArray(seed1)) == "HDF5Array")
}
\keyword{classes}
\keyword{methods}
