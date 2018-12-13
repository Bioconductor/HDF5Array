\name{TENxMatrixSeed-class}
\docType{class}

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

% sparsity

\alias{sparsity,TENxMatrixSeed-method}
\alias{is_sparse,TENxMatrixSeed-method}
\alias{extract_sparse_array,TENxMatrixSeed-method}
\alias{read_sparse_block,TENxMatrixSeed-method}
\alias{extractNonzeroDataByCol,TENxMatrixSeed-method}

\title{TENxMatrixSeed objects}

\description{
  TENxMatrixSeed is a low-level helper class for representing a pointer
  to an HDF5-based sparse matrix like one used by 10x Genomics for the
  1.3 Million Brain Cell Dataset.
  TENxMatrixSeed objects are not intended to be used directly. Most
  end-users should create and manipulate \link{TENxMatrix} objects
  instead. See \code{?\link{TENxMatrix}} for more information.
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
  No operation can be performed directly on a TENxMatrixSeed object.
  It first needs to be wrapped in a \link[DelayedArray]{DelayedMatrix}
  object. The result of this wrapping is a \link{TENxMatrix} object
  (a \link{TENxMatrix} object is just a TENxMatrixSeed object wrapped
  in a \link[DelayedArray]{DelayedMatrix} object).
}

\value{
  \code{TENxMatrixSeed()} returns a TENxMatrixSeed object.

  See \code{?\link{TENxMatrix}} for the value returned by \code{sparsity()}
  and \code{extractNonzeroDataByCol()}.
}

\seealso{
  \itemize{
    \item \link{TENxMatrix} objects.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which
          TENxMatrixSeed objects are implemented.

    \item The \code{\link[TENxBrainData]{TENxBrainData}} dataset (in the
          \pkg{TENxBrainData} package).
  }
}

\examples{
## The 1.3 Million Brain Cell Dataset from 10x Genomics is available
## via ExperimentHub:
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]

## The structure of this HDF5 file can be seen using the h5ls() command
## from the rhdf5 package:
library(rhdf5)
h5ls(fname)

## The 1.3 Million Brain Cell Dataset is represented by the "mm10"
## group. We point the TENxMatrixSeed() constructor to this group
## to create a TENxMatrixSeed object representing the dataset:
x <- TENxMatrixSeed(fname, "mm10")
x
path(x)
dim(x)
sparsity(x)
}
\keyword{methods}
\keyword{classes}
