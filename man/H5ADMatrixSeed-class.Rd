\name{H5ADMatrixSeed-class}
\docType{class}

\alias{class:Dense_H5ADMatrixSeed}
\alias{Dense_H5ADMatrixSeed-class}
\alias{Dense_H5ADMatrixSeed}

\alias{class:Sparse_H5ADMatrixSeed}
\alias{Sparse_H5ADMatrixSeed-class}
\alias{Sparse_H5ADMatrixSeed}

\alias{class:H5ADMatrixSeed}
\alias{H5ADMatrixSeed-class}
\alias{H5ADMatrixSeed}

\title{H5ADMatrixSeed objects}

\description{
  *** WORK-IN-PROGRESS! ***

  H5ADMatrixSeed is a low-level helper class used to represent a pointer
  to a matrix stored in a \code{.h5ad} file.

  It is a union class with two members, Dense_H5ADMatrixSeed and
  Sparse_H5ADMatrixSeed:
  \itemize{
    \item The Dense_H5ADMatrixSeed class is used when the matrix is stored
          as a conventional HDF5 dataset in the \code{.h5ad} file. It is
          a direct entension of the \link{HDF5ArraySeed} class.

    \item The Sparse_H5ADMatrixSeed class is used when the matrix is sparse
          and stored in the CSR/CSC/Yale format ("Compressed Sparse Row")
          in the \code{.h5ad} file. It is a direct entension of the
          \link{H5SparseMatrixSeed} class.
  }

  Note that an H5ADMatrixSeed object is not intended to be used directly.
  Most end users will typically create and manipulate a higher-level
  \link{H5ADMatrix} object instead. See \code{?\link{H5ADMatrix}} for
  more information.
}

\usage{
## Constructor function:
H5ADMatrixSeed(filepath)
}

\arguments{
  \item{filepath}{
    See \code{?\link{H5ADMatrix}} for a description of this argument.
  }
}

\details{
  A Dense_H5ADMatrixSeed or Sparse_H5ADMatrixSeed object supports
  the same limited set of methods as an \link{HDF5ArraySeed} or
  \link{H5SparseMatrixSeed} object, respectively.
  See \code{?\link{HDF5ArraySeed}} and \code{?\link{H5SparseMatrixSeed}}
  for the details.
}

\value{
  \code{H5ADMatrixSeed()} returns an H5ADMatrixSeed object
  (Dense_H5ADMatrixSeed or Sparse_H5ADMatrixSeed).
}

\section{H5ADMatrixSeed vs H5ADMatrix objects}{
  In order to have access to the full set of operations that are available
  for \link[DelayedArray]{DelayedMatrix} objects, a Dense_H5ADMatrixSeed
  or Sparse_H5ADMatrixSeed object first needs to be wrapped in a
  \link[DelayedArray]{DelayedMatrix} object, typically by calling
  the \code{\link[DelayedArray]{DelayedArray}()} constructor on it.

  This is what the \code{\link{H5ADMatrix}()} constructor function does.

  Note that the result of this wrapping is an \link{H5ADMatrix} object,
  which is just an H5ADMatrixSeed object (Dense_H5ADMatrixSeed or
  Sparse_H5ADMatrixSeed) wrapped in a \link[DelayedArray]{DelayedMatrix}
  object.
}

\seealso{
  \itemize{
    \item \link{H5ADMatrix} objects.

    \item \link{HDF5ArraySeed} and \link{H5SparseMatrixSeed} objects.
  }
}

\examples{
library(zellkonverter)
h5ad_file <- system.file("extdata", "krumsiek11.h5ad",
                         package="zellkonverter")
seed <- H5ADMatrixSeed(h5ad_file)
seed
path(seed)
dim(seed)
is_sparse(seed)

DelayedArray(seed)
stopifnot(class(DelayedArray(seed)) == "H5ADMatrix")
}
\keyword{classes}
\keyword{methods}
