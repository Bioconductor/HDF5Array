\name{H5ADMatrix-class}
\docType{class}

\alias{class:H5ADMatrix}
\alias{H5ADMatrix-class}
\alias{H5ADMatrix}

\alias{DelayedArray,Dense_H5ADMatrixSeed-method}
\alias{DelayedArray,Sparse_H5ADMatrixSeed-method}

\alias{sparsity,H5ADMatrix-method}
\alias{read_sparse_block,H5ADMatrix-method}
\alias{extractNonzeroDataByCol,H5ADMatrix-method}

\alias{coerce,H5ADMatrix,dgCMatrix-method}
\alias{coerce,H5ADMatrix,sparseMatrix-method}

\title{.h5ad datasets as DelayedMatrix objects}

\description{
  *** WORK-IN-PROGRESS! ***

  The H5ADMatrix class is a \link[DelayedArray]{DelayedMatrix} subclass
  for representing and operating on the matrix stored in an \code{.h5ad}
  file.

  All the operations available for \link[DelayedArray]{DelayedMatrix}
  objects work on H5ADMatrix objects.
}

\usage{
## Constructor function:
H5ADMatrix(filepath)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the \code{.h5ad} file.
  }
}

\value{
  \code{H5ADMatrix()} returns an H5ADMatrix object.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing conventional (a.k.a.
          dense) HDF5 datasets as \link[DelayedArray]{DelayedArray} objects.

    \item \link[DelayedArray]{DelayedMatrix} objects in the \pkg{DelayedArray}
          package.

    \item The \link{H5ADMatrixSeed} helper class.
  }
}

\examples{
library(zellkonverter)
h5ad_file <- system.file("extdata", "krumsiek11.h5ad",
                         package="zellkonverter")
H5ADMatrix(h5ad_file)
}
\keyword{classes}
\keyword{methods}
