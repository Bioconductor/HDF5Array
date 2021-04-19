\name{H5ADMatrixSeed-class}
\docType{class}

\alias{class:H5ADMatrixSeed}
\alias{H5ADMatrixSeed-class}
\alias{H5ADMatrixSeed}

\alias{class:Dense_H5ADMatrixSeed}
\alias{Dense_H5ADMatrixSeed-class}
\alias{Dense_H5ADMatrixSeed}

\alias{class:CSC_H5ADMatrixSeed}
\alias{CSC_H5ADMatrixSeed-class}
\alias{CSC_H5ADMatrixSeed}

\alias{class:CSR_H5ADMatrixSeed}
\alias{CSR_H5ADMatrixSeed-class}
\alias{CSR_H5ADMatrixSeed}

\alias{dimnames,Dense_H5ADMatrixSeed-method}

\alias{t.CSC_H5ADMatrixSeed}
\alias{t,CSC_H5ADMatrixSeed-method}
\alias{t.CSR_H5ADMatrixSeed}
\alias{t,CSR_H5ADMatrixSeed-method}

\title{H5ADMatrixSeed objects}

\description{
  H5ADMatrixSeed is a low-level helper class used to represent a pointer
  to the central matrix stored of an \code{h5ad} file, or to one of the
  matrices in the \code{/layers} group.

  It is a virtual class with three concrete subclasses: Dense_H5ADMatrixSeed,
  CSC_H5ADMatrixSeed, and CSR_H5ADMatrixSeed:
  \itemize{
    \item The Dense_H5ADMatrixSeed class is used when the matrix is stored
          as a conventional HDF5 dataset in the \code{h5ad} file. It is
          a direct entension of the \link{HDF5ArraySeed} class.

    \item The CSC_H5ADMatrixSeed or CSR_H5ADMatrixSeed classes is used
          when the matrix is stored in the \emph{Compressed Sparse Column}
          or \emph{Compressed Sparse Row} format in the \code{h5ad} file.
          CSC_H5ADMatrixSeed is a direct entension of
          \link{CSC_H5SparseMatrixSeed}, and CSR_H5ADMatrixSeed a
          direct entension of \link{CSR_H5SparseMatrixSeed}.
  }

  Note that an H5ADMatrixSeed derivative is not intended to be used directly.
  Most end users will typically create and manipulate a higher-level
  \link{H5ADMatrix} object instead. See \code{?\link{H5ADMatrix}} for
  more information.
}

\usage{
## Constructor function:
H5ADMatrixSeed(filepath, layer=NULL)
}

\arguments{
  \item{filepath, layer}{
    See \code{?\link{H5ADMatrix}} for a description of these arguments.
  }
}

\details{
  Dense_H5ADMatrixSeed objects support the same limited set of methods as
  \link{HDF5ArraySeed} objects, and CSC_H5ADMatrixSeed and CSR_H5ADMatrixSeed
  objects support the same limited set of methods as \link{H5SparseMatrixSeed}
  objects.
  See \code{?\link{HDF5ArraySeed}} and \code{?\link{H5SparseMatrixSeed}}
  for the details.
}

\value{
  \code{H5ADMatrixSeed()} returns an H5ADMatrixSeed derivative
  (Dense_H5ADMatrixSeed or CSC_H5ADMatrixSeed or CSR_H5ADMatrixSeed)
  of shape #variables x #observations.
}

\section{H5ADMatrixSeed vs H5ADMatrix objects}{
  In order to have access to the full set of operations that are available
  for \link[DelayedArray]{DelayedMatrix} objects, an H5ADMatrixSeed
  derivative first needs to be wrapped in a \link[DelayedArray]{DelayedMatrix}
  object, typically by calling the \code{\link[DelayedArray]{DelayedArray}()}
  constructor on it.

  This is what the \code{\link{H5ADMatrix}()} constructor function does.

  Note that the result of this wrapping is an \link{H5ADMatrix} object,
  which is just an H5ADMatrixSeed derivative wrapped in a
  \link[DelayedArray]{DelayedMatrix} object.
}

\references{
  \url{https://anndata.readthedocs.io/} for AnnData Python objects
  and the \code{h5ad} format.
}

\seealso{
  \itemize{
    \item \link{H5ADMatrix} objects.

    \item \link{HDF5ArraySeed} and \link{H5SparseMatrixSeed} objects.

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
