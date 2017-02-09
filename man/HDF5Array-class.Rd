\name{HDF5Array-class}
\docType{class}

% HDF5ArraySeed

\alias{class:HDF5ArraySeed}
\alias{HDF5ArraySeed-class}
\alias{HDF5ArraySeed}

\alias{dim,HDF5ArraySeed-method}
\alias{subset_seed_as_array,HDF5ArraySeed-method}

% HDF5Array and HDF5Matrix

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{coerce,HDF5Array,HDF5Matrix-method}
\alias{matrixClass,HDF5Array-method}

\alias{coerce,ANY,HDF5Matrix-method}

\title{HDF5 datasets as array-like objects}

\description{
  We provide 2 classes for representing an (on-disk) HDF5 dataset as an
  array-like object in R:
  \itemize{
    \item HDF5Array: A high-level class HDF5Array that extends
          \link{DelayedArray}. All the operations available on
          \link{DelayedArray} objects work on HDF5Array objects.

    \item HDF5ArraySeed: A low-level class for pointing to an HDF5 dataset. No
          operation can be performed directly on an HDF5ArraySeed object. It
          needs to be wrapped in a \link{DelayedArray} or HDF5Array object
          first. An HDF5Array object is just an HDF5ArraySeed object wrapped
          in a \link{DelayedArray} object.
  }
}

\usage{
## Constructor functions:
HDF5Array(file, name, type=NA)
HDF5ArraySeed(file, name, type=NA)
}

\arguments{
  \item{file}{
    The path (as a single character string) to the HDF5 file where the dataset
    is located.
  }
  \item{name}{
    The name of the dataset in the HDF5 file.
  }
  \item{type}{
    \code{NA} or the \emph{R atomic type} (specified as a single string)
    corresponding to the type of the HDF5 dataset.
  }
}

\value{
  An HDF5Array object for \code{HDF5Array()}.

  An HDF5ArraySeed object for \code{HDF5ArraySeed()}.
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \code{\link{writeHDF5Array}} for writting an array-like object
          to an HDF5 file and other low-level utilities to control the
          location of automatically created HDF5 datasets.

    \item \code{\link[SummarizedExperiment]{saveHDF5SummarizedExperiment}} and
          \code{\link[SummarizedExperiment]{loadHDF5SummarizedExperiment}} in
          the \pkg{SummarizedExperiment} package for saving/loading a
          HDF5-based SummarizedExperiment object to/from disk.

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
## SummarizedExperiment OBJECTS WITH DELAYED ASSAYS
## ---------------------------------------------------------------------

## DelayedArray objects can be used inside a SummarizedExperiment object
## to hold the assay data and to delay operations on them.
 
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
assays(se)$pcov
assays(se)$mcov
}
\keyword{methods}
\keyword{classes}
