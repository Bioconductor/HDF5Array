\name{H5ADMatrix-class}
\docType{class}

\alias{class:H5ADMatrix}
\alias{H5ADMatrix-class}
\alias{H5ADMatrix}

\alias{DelayedArray,H5ADMatrixSeed-method}

\alias{nzcount,H5ADMatrix-method}
\alias{extractNonzeroDataByCol,H5ADMatrix-method}
\alias{extractNonzeroDataByRow,H5ADMatrix-method}

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
    Can also be an \link[h5mread]{H5File} object.

    Note that you must create and use an \link[h5mread]{H5File} object if
    the HDF5 file to access is stored in an Amazon S3 bucket.
    See \code{?\link[h5mread]{H5File}} in the \pkg{h5mread} package for
    how to do this.

    Also please note that \link[h5mread]{H5File} objects must NOT be used
    in the context of parallel evaluation at the moment.
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

    \item \link[SparseArray]{SparseArray} objects in the \pkg{SparseArray}
          package.
  }
}

\examples{
library(zellkonverter)
h5ad_file <- system.file("extdata", "krumsiek11.h5ad",
                         package="zellkonverter")
X <- H5ADMatrix(h5ad_file)
X

class(X)  # H5ADMatrix
is(X, "DelayedMatrix")  # TRUE

class(seed(X))  # Dense_H5ADMatrixSeed

dim(X)
path(X)
is_sparse(X)  # FALSE

## Use coercion to load the full dataset into memory:
as.matrix(X)          # as ordinary array (usually not recommended)

\dontrun{
  ## Works only if H5ADMatrix object is sparse!
  as(X, "dgCMatrix")    # as dgCMatrix
  as(X, "SparseArray")  # as SparseArray object (most efficient)
  SparseArray(X)        # equivalent to 'as(X, "SparseArray")'
}
}
\keyword{classes}
\keyword{methods}
