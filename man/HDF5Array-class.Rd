\name{HDF5Array-class}
\docType{class}

% HDF5Dataset

\alias{class:HDF5Dataset}
\alias{HDF5Dataset-class}
\alias{HDF5Dataset}

\alias{dim,HDF5Dataset-method}

\alias{coerce,array,HDF5Dataset-method}

% HDF5Array and HDF5Matrix

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{coerce,HDF5Dataset,HDF5Array-method}
\alias{coerce,array,HDF5Array-method}
\alias{coerce,HDF5Array,HDF5Matrix-method}
\alias{coerce,array,HDF5Matrix-method}

\alias{coerce,DelayedArray,HDF5Array-method}
\alias{coerce,DelayedArray,HDF5Matrix-method}
\alias{coerce,DelayedMatrix,HDF5Matrix-method}

\title{HDF5 datasets as array-like objects}

\description{
  We provide 3 classes for representing an (on-disk) HDF5 dataset as an
  array-like object in R:
  \itemize{
    \item Two high-level classes, HDF5Array and HDF5Matrix, that extend
          \link{DelayedArray} and \link{DelayedMatrix}, respectively.
          All the operations available on \link{DelayedArray} and
          \link{DelayedMatrix} objects work on HDF5Array and HDF5Matrix
          objects.

    \item A low-level class, HDF5Dataset, only for pointing to an HDF5
          dataset. No operation can be performed directly on an HDF5Dataset
          object. It needs to be wrapped in a DelayedArray object first.
          An HDF5Array (or HDF5Matrix) object is just an HDF5Dataset object
          wrapped in a \link{DelayedArray} (or \link{DelayedMatrix}) object.
  }
}

\usage{
## Constructor functions
HDF5Array(file, group, name, type=NA)
HDF5Matrix(file, group, name, type=NA)
HDF5Dataset(file, group, name, type=NA)
}

\arguments{
  \item{file}{
    The path to the HDF5 file where the dataset is located.
  }
  \item{group, name}{
    The HDF5 group and name of the dataset, as reported by
    \code{\link[rhdf5]{h5ls}}.
  }
  \item{type}{
    \code{NA} or the \emph{R atomic type} (specified as a single string)
    corresponding to the type of the HDF5 dataset.
  }
}

\value{
  An HDF5Array or HDF5Matrix object.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \link{HDF5_output_settings} to manage settings for writing
          output to HDF5.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which the
          \pkg{HDF5Array} package is implemented.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## CONSTRUCTION
## ---------------------------------------------------------------------
library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

## Pick up "Coverages" dataset for Human chromosome 16:
cov0 <- HDF5Array(tally_file, "/ExampleStudy/16", "Coverages")
cov0

## ---------------------------------------------------------------------
## dim/dimnames
## ---------------------------------------------------------------------
dim(cov0)

dimnames(cov0)
dimnames(cov0) <- list(paste0("s", 1:6), c("+", "-"))
dimnames(cov0)

## ---------------------------------------------------------------------
## SLICING (A.K.A. SUBSETTING)
## ---------------------------------------------------------------------
cov1 <- cov0[ , , 29000001:29000007]
cov1

dim(cov1)
as.array(cov1)
stopifnot(identical(dim(as.array(cov1)), dim(cov1)))
stopifnot(identical(dimnames(as.array(cov1)), dimnames(cov1)))

cov2 <- cov0[ , "+", 29000001:29000007]
cov2
as.matrix(cov2)

## ---------------------------------------------------------------------
## COERCION FROM HDF5Array TO DelayedMatrix
## ---------------------------------------------------------------------
library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

## Pick up "Coverages" dataset for Human chromosome 16:
cov0 <- HDF5Array(tally_file, "/ExampleStudy/16", "Coverages")
dimnames(cov0) <- list(paste0("s", 1:6), c("+", "-"))

pcov <- as(cov0[ , 1, ], "DelayedMatrix")  # coverage on plus strand
mcov <- as(cov0[ , 2, ], "DelayedMatrix")  # coverage on minus strand

nrow(pcov)  # nb of samples
ncol(pcov)  # length of Human chromosome 16

## ---------------------------------------------------------------------
## DelayedMatrix OBJECTS AS ASSAYS OF A SummarizedExperiment OBJECT
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
