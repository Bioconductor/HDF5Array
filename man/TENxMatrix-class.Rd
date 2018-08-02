\name{TENxMatrix-class}
\docType{class}

% TENxMatrixSeed

\alias{class:TENxMatrixSeed}
\alias{TENxMatrixSeed-class}
\alias{TENxMatrixSeed}

\alias{path,TENxMatrixSeed-method}
\alias{path<-,TENxMatrixSeed-method}
\alias{dim,TENxMatrixSeed-method}
\alias{dimnames,TENxMatrixSeed-method}
\alias{extract_array,TENxMatrixSeed-method}
\alias{chunkdim,TENxMatrixSeed-method}

% TENxMatrix

\alias{class:TENxMatrix}
\alias{TENxMatrix-class}

\alias{DelayedArray,TENxMatrixSeed-method}
\alias{TENxMatrix}

\title{10x Genomics datasets as DelayedArray objects}

\description{
  The 1.3 Million Brain Cell Dataset and other datasets published by
  10x Genomics use an HDF5-based sparse matrix representation instead of
  the conventional (i.e. dense) HDF5 representation.

  We provide 2 classes for representing a 10x Genomics dataset as a
  matrix-like object in R:
  \itemize{
    \item TENxMatrix: A high-level class that extends \link{DelayedArray}.
          All the operations available for \link{DelayedArray} objects work
          on TENxMatrix objects.

    \item TENxMatrixSeed: A low-level class for pointing to a 10x Genomics
          dataset.
          No operation can be performed directly on a TENxMatrixSeed object.
          It needs to be wrapped in a \link{DelayedArray} object first.
          A TENxMatrix object is just a TENxMatrixSeed object wrapped in
          a \link{DelayedArray} object.
  }
}

\usage{
## Constructor functions:
TENxMatrix(filepath, group="mm10")
TENxMatrixSeed(filepath, group="mm10")
}

\arguments{
  \item{filepath}{
    The path (as a single character string) to the HDF5 file where the
    10x Genomics dataset is located.
  }
  \item{group}{
    The name of the group in the HDF5 file containing the 10x Genomics data.
  }
}

\value{
  A TENxMatrix object for \code{TENxMatrix()}.

  A TENxMatrixSeed object for \code{TENxMatrixSeed()}.
}

\note{
  If your dataset uses the HDF5-based sparse matrix representation from
  10x Genomics, use the \code{TENxMatrix()} constructor.

  If your dataset uses the conventional (i.e. dense) HDF5 representation,
  use the \code{\link{HDF5Array}()} constructor.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{HDF5Array} objects for representing conventional (i.e.
          dense) HDF5 datasets as \link{DelayedArray} objects.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which TENxMatrix
          and TENxMatrixSeed objects are implemented.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
}
\keyword{methods}
\keyword{classes}
