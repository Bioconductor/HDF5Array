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

\alias{show,TENxMatrixSeed-method}

% TENxMatrix

\alias{class:TENxMatrix}
\alias{TENxMatrix-class}
\alias{TENxMatrix}

\alias{DelayedArray,TENxMatrixSeed-method}

% extractNonZeroValuesByCol

\alias{extractNonZeroValuesByCol}
\alias{extractNonZeroValuesByCol,TENxMatrixSeed-method}
\alias{extractNonZeroValuesByCol,TENxMatrix-method}

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

## A convenient data extractor:
extractNonZeroValuesByCol(x, j)
}

\arguments{
  \item{filepath}{
    The path (as a single character string) to the HDF5 file where the
    10x Genomics dataset is located.
  }
  \item{group}{
    The name of the group in the HDF5 file containing the 10x Genomics data.
  }
  \item{x}{
    A TENxMatrix or TENxMatrixSeed object.
  }
  \item{j}{
    An integer vector containing valid column indices.
  }
}

\value{
  \code{TENxMatrix}: A TENxMatrix object.

  \code{TENxMatrixSeed}: A TENxMatrixSeed object.

  \code{extractNonZeroValuesByCol}: A \link[IRanges]{NumericList} or
  \link[IRanges]{IntegerList} object \emph{parallel} to \code{j} i.e.
  with one list element per column index in \code{j}.
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
  }
}

\examples{
## ---------------------------------------------------------------------
## THE "1.3 Million Brain Cell Dataset" AS A DelayedArray OBJECT
## ---------------------------------------------------------------------
## The 1.3 Million Brain Cell Dataset from 10x Genomics is available via
## ExperimentHub:
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]

## The structure of the file can be seen using the h5ls() command from
## the rhdf5 package:
library(rhdf5)
h5ls(fname)

## The 1.3 Million Brain Cell Dataset is represented by the "mm10"
## group. We point the TENxMatrix() constructor to this group to
## create a TENxMatrix object representing the dataset:
oneM <- TENxMatrix(fname, "mm10")
oneM

is(oneM, "DelayedArray")  # TRUE

## extractNonZeroValuesByCol() is an efficient way to extract the
## non-zero values in a compact form:
counts_by_col <- extractNonZeroValuesByCol(oneM, 1:2500)
counts_by_col

## This provides a more efficient way to compute column summaries (such
## as the sum or the square sums) than using something like colSums(oneM)
## colSums(oneM * oneM):
head(lengths(counts_by_col))      # nb of non-zero counts per column
head(sum(counts_by_col))          # total count per column
head(sum(counts_by_col * counts_by_col))  # total square count per column
}
\keyword{methods}
\keyword{classes}
