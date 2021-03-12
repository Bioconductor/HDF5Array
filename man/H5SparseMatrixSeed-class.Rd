\name{H5SparseMatrixSeed-class}
\docType{class}

\alias{class:H5SparseMatrixSeed}
\alias{H5SparseMatrixSeed-class}
\alias{H5SparseMatrixSeed}

\alias{class:CSC_H5SparseMatrixSeed}
\alias{CSC_H5SparseMatrixSeed-class}
\alias{CSC_H5SparseMatrixSeed}

\alias{class:CSR_H5SparseMatrixSeed}
\alias{CSR_H5SparseMatrixSeed-class}
\alias{CSR_H5SparseMatrixSeed}

\alias{t.CSC_H5SparseMatrixSeed}
\alias{t,CSC_H5SparseMatrixSeed-method}
\alias{t.CSR_H5SparseMatrixSeed}
\alias{t,CSR_H5SparseMatrixSeed-method}

\alias{path,H5SparseMatrixSeed-method}
\alias{path<-,H5SparseMatrixSeed-method}
\alias{dim,H5SparseMatrixSeed-method}
\alias{dimnames,H5SparseMatrixSeed-method}
\alias{extract_array,H5SparseMatrixSeed-method}
\alias{chunkdim,CSC_H5SparseMatrixSeed-method}
\alias{chunkdim,CSR_H5SparseMatrixSeed-method}

\alias{sparsity}
\alias{sparsity,H5SparseMatrixSeed-method}
\alias{is_sparse,H5SparseMatrixSeed-method}
\alias{extract_sparse_array,H5SparseMatrixSeed-method}
\alias{read_sparse_block,H5SparseMatrixSeed-method}
\alias{extractNonzeroDataByCol}
\alias{extractNonzeroDataByCol,CSC_H5SparseMatrixSeed-method}
\alias{extractNonzeroDataByRow}
\alias{extractNonzeroDataByRow,CSR_H5SparseMatrixSeed-method}

\alias{coerce,CSC_H5SparseMatrixSeed,dgCMatrix-method}
\alias{coerce,CSC_H5SparseMatrixSeed,sparseMatrix-method}
\alias{coerce,CSR_H5SparseMatrixSeed,dgCMatrix-method}
\alias{coerce,CSR_H5SparseMatrixSeed,sparseMatrix-method}

\alias{show,H5SparseMatrixSeed-method}

\title{H5SparseMatrixSeed objects}

\description{
  H5SparseMatrixSeed is a low-level helper class for representing a
  pointer to an HDF5 sparse matrix stored in CSR/CSC/Yale format.

  It is a virtual class with two concrete subclasses: CSC_H5SparseMatrixSeed
  for the \emph{Compressed Sparse Column} format and CSR_H5SparseMatrixSeed
  for the \emph{Compressed Sparse Row} format.
  The former is used by 10x Genomics (e.g. "1.3 Million Brain Cell Dataset").
  \code{h5ad} files can use one or the other format to store a sparse matrix.

  Note that an H5SparseMatrixSeed derivative is not intended to be used
  directly. Most end users will typically create and manipulate a
  higher-level object like a \link{TENxMatrix} or \link{H5ADMatrix}
  object instead. See \code{?\link{TENxMatrix}} and \code{?\link{H5ADMatrix}}
  for more information.
}

\usage{
## --- Constructor function ---

H5SparseMatrixSeed(filepath, group)

## --- Accessors --------------

\S4method{path}{H5SparseMatrixSeed}(object)

\S4method{dim}{H5SparseMatrixSeed}(x)

\S4method{dimnames}{H5SparseMatrixSeed}(x)

\S4method{chunkdim}{CSC_H5SparseMatrixSeed}(x)
\S4method{chunkdim}{CSR_H5SparseMatrixSeed}(x)

## --- Data extraction --------

\S4method{extract_array}{H5SparseMatrixSeed}(x, index)

\S4method{extract_sparse_array}{H5SparseMatrixSeed}(x, index)

\S4method{read_sparse_block}{H5SparseMatrixSeed}(x, viewport)

\S4method{extractNonzeroDataByCol}{CSC_H5SparseMatrixSeed}(x, j)
\S4method{extractNonzeroDataByRow}{CSR_H5SparseMatrixSeed}(x, i)

## --- Other methods ----------

\S4method{is_sparse}{H5SparseMatrixSeed}(x)

\S4method{sparsity}{H5SparseMatrixSeed}(x)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the HDF5 file (\code{.h5} or
    \code{h5ad}) where the sparse matrix is located.
  }
  \item{group}{
    The name of the group in the HDF5 file that stores the sparse matrix.
  }
  \item{object, x}{
    An H5SparseMatrixSeed derivative.
  }
  \item{index}{
    See \code{?\link[DelayedArray]{extract_array}} in the \pkg{DelayedArray}
    package.
  }
  \item{viewport}{
    See \code{?\link[DelayedArray]{read_block}} in the \pkg{DelayedArray}
    package.
  }
  \item{j}{
    An integer vector containing valid column indices.
  }
  \item{i}{
    An integer vector containing valid row indices.
  }
}

\details{
  Note that the implementation of CSC_H5SparseMatrixSeed and
  CSR_H5SparseMatrixSeed objects follows the usual convention of
  transposing the HDF5 matrix when loading it into R.
  This means that a CSC_H5SparseMatrixSeed object represents a
  sparse matrix stored in the CSR format (Compressed Sparse Row)
  at the HDF5 level, and a CSR_H5SparseMatrixSeed object represents a
  sparse matrix stored in the CSC format (Compressed Sparse Column)
  at the HDF5 level.

  Finally note that H5SparseMatrixSeed derivatives support a very limited
  set of methods:
  \itemize{
    \item \code{path()}: Returns the path to the HDF5 file where the sparse
          matrix is located.
    \item \code{dim()}, \code{dimnames()}.
    \item \code{chunkdim()}, \code{extract_array()},
          \code{extract_sparse_array()}, \code{read_sparse_block()},
          \code{is_sparse()}: These generics are defined and documented
          in the \pkg{DelayedArrray} package.
    \item \code{sparsity()}: Returns the number of zero-valued matrix elements
          in the object divided by its total number of elements (a.k.a. its
          length).
    \item \code{extractNonzeroDataByCol()}: Works on CSC_H5SparseMatrixSeed
          objects only. Returns a \link[IRanges]{NumericList} or
          \link[IRanges]{IntegerList} object \emph{parallel} to \code{j},
          that is, with one list element per column index in \code{j}.
          The row indices of the values are not returned. Furthermore, the
          values within a given list element can be returned in **any order**.
          In particular, do NOT assume that they are ordered by ascending
          row index.
    \item \code{extractNonzeroDataByRow()}: Works on CSR_H5SparseMatrixSeed
          objects only. Returns a \link[IRanges]{NumericList} or
          \link[IRanges]{IntegerList} object \emph{parallel} to \code{i},
          that is, with one list element per row index in \code{i}.
          The column indices of the values are not returned. Furthermore, the
          values within a given list element can be returned in **any order**.
          In particular, do NOT assume that they are ordered by ascending
          column index.
  }

  In order to have access to the full set of operations that
  are available for \link[DelayedArray]{DelayedMatrix} objects,
  an H5SparseMatrixSeed derivative would first need to be wrapped in
  a \link[DelayedArray]{DelayedMatrix} object, typically by calling
  the \code{\link[DelayedArray]{DelayedArray}()} constructor on it.
}

\value{
  \code{H5SparseMatrixSeed()} returns an H5SparseMatrixSeed derivative
  (CSC_H5SparseMatrixSeed or CSR_H5SparseMatrixSeed object).
}

\references{
  \url{https://en.wikipedia.org/wiki/Sparse_matrix} for a description
  of the CSR/CSC/Yale format (section "Compressed sparse row (CSR, CRS
  or Yale format)").
}

\seealso{
  \link{TENxMatrixSeed} and \link{H5ADMatrixSeed} objects.
}

\examples{
showClass("H5SparseMatrixSeed")
}
\keyword{classes}
\keyword{methods}
