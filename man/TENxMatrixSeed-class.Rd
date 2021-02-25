\name{TENxMatrixSeed-class}
\docType{class}

\alias{class:TENxMatrixSeed}
\alias{TENxMatrixSeed-class}
\alias{TENxMatrixSeed}

\title{TENxMatrixSeed objects}

\description{
  TENxMatrixSeed is a low-level helper class that is a direct extension
  of the \link{H5SparseMatrixSeed} class. It is used to represent a
  pointer to an HDF5 sparse matrix that is stored in the CSR/CSC/Yale
  format ("Compressed Sparse Row") and follows the 10x Genomics convention
  for storing the dimensions of the matrix.

  Note that a TENxMatrixSeed object is not intended to be used directly.
  Most end users will typically create and manipulate a higher-level
  \link{TENxMatrix} object instead. See \code{?\link{TENxMatrix}} for
  more information.
}

\usage{
## Constructor function:
TENxMatrixSeed(filepath, group="mm10")
}

\arguments{
  \item{filepath, group}{
    See \code{?\link{TENxMatrix}} for a description of these arguments.
  }
}

\details{
  A TENxMatrixSeed object supports the same limited set of methods as
  an \link{H5SparseMatrixSeed} object. See \code{?\link{H5SparseMatrixSeed}}
  for the details.
}

\value{
  \code{TENxMatrixSeed()} returns a TENxMatrixSeed object.
}

\section{TENxMatrixSeed vs TENxMatrix objects}{
  In order to have access to the full set of operations that are available
  for \link[DelayedArray]{DelayedMatrix} objects, a TENxMatrixSeed object
  first needs to be wrapped in a \link[DelayedArray]{DelayedMatrix} object,
  typically by calling the \code{\link[DelayedArray]{DelayedArray}()}
  constructor on it.

  This is what the \code{\link{TENxMatrix}()} constructor function does.

  Note that the result of this wrapping is a \link{TENxMatrix} object,
  which is just a TENxMatrixSeed object wrapped in a
  \link[DelayedArray]{DelayedMatrix} object.
}

\seealso{
  \itemize{
    \item \link{TENxMatrix} objects.

    \item \link{H5SparseMatrixSeed} objects.

    \item The \code{\link[TENxBrainData]{TENxBrainData}} dataset (in the
          \pkg{TENxBrainData} package).

    \item \code{\link{h5ls}} to list the content of an HDF5 file.
  }
}

\examples{
## The 1.3 Million Brain Cell Dataset from 10x Genomics is available
## via ExperimentHub:
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]

## 'fname' is an HDF5 file. Use h5ls() to list its content:
h5ls(fname)

## The 1.3 Million Brain Cell Dataset is represented by the "mm10"
## group. We point the TENxMatrixSeed() constructor to this group
## to create a TENxMatrixSeed object representing the dataset:
seed <- TENxMatrixSeed(fname, "mm10")
seed
path(seed)
dim(seed)
is_sparse(seed)
sparsity(seed)

DelayedArray(seed)
stopifnot(class(DelayedArray(seed)) == "TENxMatrix")
}
\keyword{classes}
\keyword{methods}
