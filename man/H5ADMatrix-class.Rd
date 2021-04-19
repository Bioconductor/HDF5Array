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

\title{h5ad central matrices (or matrices in the /layers group)
       as DelayedMatrix objects}

\description{
  \code{h5ad} files are HDF5 files used for on-disk representation
  of AnnData Python objects. At the very minimum, they contain a central
  data matrix, named \code{X}, of shape #observations x #variables, and
  possibly additional data matrices (stored in the HDF5 group \code{/layers})
  that share the shape and dimnames of \code{X}.
  See \url{https://anndata.readthedocs.io/} for more information.

  The H5ADMatrix class is a \link[DelayedArray]{DelayedMatrix} subclass
  for representing and operating on the central matrix of an \code{h5ad}
  file, or any matrix in its \code{/layers} group.

  All the operations available for \link[DelayedArray]{DelayedMatrix}
  objects work on H5ADMatrix objects.
}

\usage{
## Constructor function:
H5ADMatrix(filepath, layer=NULL)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the \code{h5ad} file.
  }
  \item{layer}{
    \code{NULL} (the default) or the name of a matrix in the \code{/layers}
    group. By default (i.e. when \code{layer} is not specified)
    \code{H5ADMatrix()} returns the central matrix (\code{X}).
  }
}

\value{
  \code{H5ADMatrix()} returns an H5ADMatrix object of shape #variables x
  #observations. Note that in Python and HDF5 the shape of this matrix is
  considered to be #observations x #variables, but in R it is transposed.
  This follows the widely adopted convention of transposing HDF5 matrices
  when they get loaded into R.
}

\references{
  \url{https://anndata.readthedocs.io/} for AnnData Python objects
  and the \code{h5ad} format.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing conventional (a.k.a.
          dense) HDF5 datasets as \link[DelayedArray]{DelayedArray} objects.

    \item \link{H5SparseMatrix} objects for representing HDF5 sparse matrices
          as \link[DelayedArray]{DelayedMatrix} objects.

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
X <- H5ADMatrix(h5ad_file)
X
}
\keyword{classes}
\keyword{methods}
