\name{HDF5Matrix-class}
\docType{class}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{coerce,HDF5Array,HDF5Matrix-method}
\alias{coerce,matrix,HDF5Matrix-method}

\alias{show,HDF5Matrix-method}

\title{HDF5Matrix objects}

\description{
  The HDF5Matrix class is a container for representing and manipulating an
  on-disk HDF5 dataset with 2 dimensions, or a 2-dimensional slice of an
  on-disk HDF5 dataset with 2 or more dimensions.

  HDF5Matrix objects can be manipulated via a \emph{restricted matrix-like
  interface} i.e. they support common \link[base]{matrix} operations like
  \code{dim}, \code{nrow}, \code{ncol}, \code{dimnames}, \code{rownames},
  \code{colnames}, \code{[}, \code{rowSums}, \code{colSums}, etc... as long
  as the operation doesn't need to alter the content of the underlying
  (on-disk) dataset or write a new dataset to the disk.
}

\usage{
HDF5Matrix(file, group, name)  # constructor function
}

\arguments{
  \item{file}{
    The path to the HDF5 file where the dataset is located.
  }
  \item{group, name}{
    The HDF5 group and name of the dataset, as reported by
    \code{\link[rhdf5]{h5ls}}.
  }
}

\value{
  A HDF5Matrix object.
}

\section{Accessors}{
}

\section{Coercion}{
}

\section{Subsetting}{
}

\seealso{
  \itemize{
    \item \link[base]{matrix} objects in base R.

    \item \link{HDF5Array} objects.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which the
          \pkg{HDF5Array} package is implemented.
  }
}

\examples{
## ---------------------------------------------------------------------
## CONSTRUCTION
## ---------------------------------------------------------------------


## ---------------------------------------------------------------------
## COERCION FROM HDF5Array TO HDF5Matrix
## ---------------------------------------------------------------------
library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

## Pick up "Coverages" dataset for Human chromosome 16:
cov0 <- HDF5Array(tally_file, "/ExampleStudy/16", "Coverages")
dimnames(cov0) <- list(paste0("s", 1:6), c("+", "-"))

pcov <- as(cov0[ , 1, ], "HDF5Matrix")  # coverage on plus strand
mcov <- as(cov0[ , 2, ], "HDF5Matrix")  # coverage on minus strand

nrow(pcov)  # nb of samples
ncol(pcov)  # length of Human chromosome 16

## ---------------------------------------------------------------------
## HDF5Matrix OBJECTS AS ASSAYS OF A SummarizedExperiment OBJECT
## ---------------------------------------------------------------------
library(SummarizedExperiment)

## The convention for a SummarizedExperiment object is to have 1 column
## per sample so first we need to transpose 'pcov' and 'mcov':
pcov <- t(pcov)
mcov <- t(mcov)
se <- SummarizedExperiment(list(pcov=pcov, mcov=mcov))
se
stopifnot(validObject(se, complete=TRUE))

## A GPos object can be used to represent the genomic positions along
## the dataset:
gpos <- GPos(GRanges("16", IRanges(1, nrow(se))))
gpos
rowRanges(se) <- gpos
se
stopifnot(validObject(se))
}
\keyword{methods}
\keyword{classes}
