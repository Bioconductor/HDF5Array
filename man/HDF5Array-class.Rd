\name{HDF5Array-class}
\docType{class}

% HDF5Dataset

\alias{class:HDF5Dataset}
\alias{HDF5Dataset-class}
\alias{HDF5Dataset}

\alias{dim,HDF5Dataset-method}
\alias{writeHDF5Dataset}

% HDF5Array and HDF5Matrix

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{coerce,HDF5Array,HDF5Matrix-method}

\alias{coerce,DelayedArray,HDF5Array-method}
\alias{coerce,DelayedArray,HDF5Matrix-method}
\alias{coerce,DelayedMatrix,HDF5Matrix-method}

\title{HDF5 datasets as array-like objects}

\description{
  We provide 2 classes for representing an (on-disk) HDF5 dataset as an
  array-like object in R:
  \itemize{
    \item HDF5Array: A high-level class HDF5Array that extends
          \link{DelayedArray}. All the operations available on
          \link{DelayedArray} objects work on HDF5Array objects.

    \item HDF5Dataset: A low-level class for pointing to an HDF5 dataset. No
          operation can be performed directly on an HDF5Dataset object. It
          needs to be wrapped in a \link{DelayedArray} or HDF5Array object
          first. An HDF5Array object is just an HDF5Dataset object wrapped
          in a \link{DelayedArray} object.
  }
}

\usage{
## Constructor functions
HDF5Array(file, name, type=NA)
HDF5Dataset(file, name, type=NA)

## Write an array-like object to an HDF5 file
writeHDF5Dataset(x, file, name)
}

\arguments{
  \item{file}{
    For \code{HDF5Array} and \code{HDF5Dataset}: The path (as a single
    character string) to the HDF5 file where the dataset is located.
    Alternatively \code{file} can be a \link{DelayedArray} object or
    an ordinary array, in which case the object is written to disk
    as a new HDF5 dataset by calling \code{writeHDF5Dataset} internally.
    Note that, when given a \link{DelayedArray} object, \code{writeHDF5Dataset}
    \emph{realizes} it on disk, that is, all the delayed operations carried
    by the object are executed while the object is written to disk.
    See "On-disk realization of a DelayedArray object as an HDF5 dataset"
    section below for more information.

    For \code{writeHDF5Dataset}: The path (as a single character string)
    to the (new or existing) HDF5 file where to write the dataset.
  }
  \item{name}{
    For \code{HDF5Array} and \code{HDF5Dataset}:
    The name of the dataset in the HDF5 file.

    For \code{writeHDF5Dataset}:
    The name of the HDF5 dataset to write.
  }
  \item{type}{
    \code{NA} or the \emph{R atomic type} (specified as a single string)
    corresponding to the type of the HDF5 dataset.
  }
  \item{x}{
    The array-like object to write to an HDF5 file. Can be an ordinary array,
    or a \link{DelayedArray} or \link{HDF5Dataset} object at the moment.

    If \code{x} is a \link{DelayedArray} object, \code{writeHDF5Dataset}
    \emph{realizes} it on disk, that is, all the delayed operations carried
    by the object are executed while the object is written to disk.
    See "On-disk realization of a DelayedArray object as an HDF5 dataset"
    section below for more information.

    If \code{x} is an \link{HDF5Dataset} object, \code{writeHDF5Dataset}
    first wraps it in a \link{DelayedArray} object in order to trigger the
    block-processing strategy used by the on-disk realization mechanism.
  }
}

\details{
  The \code{HDF5Array} and \code{HDF5Dataset} constructor functions can be
  used either to point to an existing HDF5 dataset or to create a new one
  (see description of the \code{file} argument above).

  When used to create a new HDF5 dataset, the location where to write the
  dataset can be controlled with the \code{\link{setHDF5DumpFile}} and
  \code{\link{setHDF5DumpName}} utility functions.
}

\value{
  An HDF5Array object for \code{HDF5Array()}.

  An HDF5Dataset object for \code{HDF5Dataset()}.

  An invisible HDF5Dataset object for \code{writeHDF5Dataset()}.
}

\section{On-disk realization of a DelayedArray object as an HDF5 dataset}{
  When passed an \link{DelayedArray} object, \code{writeHDF5Dataset}
  \emph{realizes} it on disk, that is, all the delayed operations carried
  by the object are executed on-the-fly while the object is written to disk.
  This uses a block-processing strategy so that the full object is not
  realized at once in memory. Instead the object is processed block by block
  i.e. the blocks are realized in memory and written to disk one at a time.

  In other words, \code{writeHDF5Dataset(x, ...)} is semantically equivalent
  to \code{writeHDF5Dataset(as.array(x), ...)}, except that \code{as.array(x)}
  is not called because this would realize the full object at once in memory.

  See \code{?\link{DelayedArray}} for general information about
  \link{DelayedArray} objects.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \link{setHDF5DumpFile} to control the location of the new
          HDF5 datasets created by \code{HDF5Array} and \code{HDF5Dataset}.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which HDF5Array
          objects are implemented.

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
cov0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")
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
cov1 <- drop(cov0[ , , 29000001:29000007])
cov1

dim(cov1)
as.array(cov1)
stopifnot(identical(dim(as.array(cov1)), dim(cov1)))
stopifnot(identical(dimnames(as.array(cov1)), dimnames(cov1)))

cov2 <- drop(cov0[ , "+", 29000001:29000007])
cov2
as.matrix(cov2)

## ---------------------------------------------------------------------
## DelayedMatrix OBJECTS AS ASSAYS OF A SummarizedExperiment OBJECT
## ---------------------------------------------------------------------
library(SummarizedExperiment)

pcov <- drop(cov0[ , 1, ])  # coverage on plus strand
mcov <- drop(cov0[ , 2, ])  # coverage on minus strand

nrow(pcov)  # nb of samples
ncol(pcov)  # length of Human chromosome 16

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

## ---------------------------------------------------------------------
## writeHDF5Dataset()
## ---------------------------------------------------------------------
out_file <- tempfile()
writeHDF5Dataset(cov1, out_file, "cov1")
h5ls(out_file)
}
\keyword{methods}
\keyword{classes}
