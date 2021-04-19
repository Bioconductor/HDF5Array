\name{H5SparseMatrix-class}
\docType{class}

\alias{class:H5SparseMatrix}
\alias{H5SparseMatrix-class}
\alias{H5SparseMatrix}

\alias{DelayedArray,H5SparseMatrixSeed-method}

\alias{sparsity,H5SparseMatrix-method}
\alias{read_sparse_block,H5SparseMatrix-method}
\alias{extractNonzeroDataByCol,H5SparseMatrix-method}
\alias{extractNonzeroDataByRow,H5SparseMatrix-method}

\alias{coerce,H5SparseMatrix,dgCMatrix-method}
\alias{coerce,H5SparseMatrix,sparseMatrix-method}

\title{HDF5 sparse matrices as DelayedMatrix objects}

\description{
  The H5SparseMatrix class is a \link[DelayedArray]{DelayedMatrix} subclass
  for representing and operating on an HDF5 sparse matrix stored in
  CSR/CSC/Yale format.

  All the operations available for \link[DelayedArray]{DelayedMatrix}
  objects work on H5SparseMatrix objects.
}

\usage{
## Constructor function:
H5SparseMatrix(filepath, group)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the HDF5 file (\code{.h5} or
    \code{.h5ad}) where the sparse matrix is located.
  }
  \item{group}{
    The name of the group in the HDF5 file where the sparse matrix is stored.
  }
}

\value{
  An H5SparseMatrix object.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing conventional (a.k.a.
          dense) HDF5 datasets as \link[DelayedArray]{DelayedArray} objects.

    \item \link{H5ADMatrix} objects for representing h5ad central
          matrices (or matrices in the \code{/layers} group)
          as \link[DelayedArray]{DelayedMatrix} objects.

    \item \link{TENxMatrix} objects for representing 10x Genomics
          datasets as \link[DelayedArray]{DelayedMatrix} objects.

    \item \link[DelayedArray]{DelayedMatrix} objects in the \pkg{DelayedArray}
          package.

    \item The \link{H5SparseMatrixSeed} helper class.

    \item \code{\link{h5ls}} to list the content of an HDF5 file (\code{.h5}
          or \code{.h5ad}).
  }
}

\examples{
library(zellkonverter)
h5ad_file <- system.file("extdata", "example_anndata.h5ad",
                         package="zellkonverter")
h5ls(h5ad_file)

M <- H5SparseMatrix(h5ad_file, "/obsp/connectivities")
M

is(M, "DelayedMatrix")  # TRUE
dim(M)
seed(M)
path(M)
is_sparse(M)  # TRUE
}
\keyword{classes}
\keyword{methods}
