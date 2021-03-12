\name{H5ADMatrix-class}
\docType{class}

\alias{class:H5ADMatrix}
\alias{H5ADMatrix-class}
\alias{H5ADMatrix}

\alias{DelayedArray,H5ADMatrixSeed-method}

\alias{sparsity,H5ADMatrix-method}
\alias{read_sparse_block,H5ADMatrix-method}
\alias{extractNonzeroDataByCol,H5ADMatrix-method}
\alias{extractNonzeroDataByRow,H5ADMatrix-method}

\alias{coerce,H5ADMatrix,dgCMatrix-method}
\alias{coerce,H5ADMatrix,sparseMatrix-method}

\title{h5ad datasets as DelayedMatrix objects}

\description{
  \code{h5ad} files are HDF5 files used for on-disk representation
  of AnnData Python objects. See \url{https://anndata.readthedocs.io/}
  for more information.

  The H5ADMatrix class is a \link[DelayedArray]{DelayedMatrix} subclass
  for representing and operating on the matrix stored in an \code{h5ad}
  file.

  All the operations available for \link[DelayedArray]{DelayedMatrix}
  objects work on H5ADMatrix objects.
}

\usage{
## Constructor function:
H5ADMatrix(filepath, name="X")
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the \code{h5ad} file.
  }
  \item{name}{
    The name of the dataset or group in the \code{h5ad} file that stores
    the matrix.
  }
}

\value{
  \code{H5ADMatrix()} returns an H5ADMatrix object.
}

\references{
  \url{https://anndata.readthedocs.io/} for AnnData Python objects
  and the \code{h5ad} format.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing conventional (a.k.a.
          dense) HDF5 datasets as \link[DelayedArray]{DelayedArray} objects.

    \item \link[DelayedArray]{DelayedMatrix} objects in the \pkg{DelayedArray}
          package.

    \item The \link{H5ADMatrixSeed} helper class.

    \item \code{\link[zellkonverter]{readH5AD}} and
          \code{\link[zellkonverter]{writeH5AD}} in
          the \pkg{zellkonverter} package for
          importing/exporting an \code{h5ad} file as/from a
          \link[SingleCellExperiment]{SingleCellExperiment} object.
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
