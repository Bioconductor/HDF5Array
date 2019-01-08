\name{HDF5Array-class}
\docType{class}

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\alias{DelayedArray,HDF5ArraySeed-method}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{matrixClass,HDF5Array-method}
\alias{coerce,HDF5Array,HDF5Matrix-method}
\alias{coerce,HDF5Matrix,HDF5Array-method}
\alias{coerce,ANY,HDF5Matrix-method}

\title{HDF5 datasets as DelayedArray objects}

\description{
  The HDF5Array class is a \link[DelayedArray]{DelayedArray} subclass
  for representing a conventional (i.e. dense) HDF5 dataset.

  All the operations available for \link[DelayedArray]{DelayedArray}
  objects work on HDF5Array objects.
}

\usage{
## Constructor function:
HDF5Array(filepath, name, type=NA)
}

\arguments{
  \item{filepath}{
    The path (as a single character string) to the HDF5 file where the
    dataset is located.
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
  An HDF5Array object.
}

\note{
  The 1.3 Million Brain Cell Dataset and other datasets published by
  10x Genomics use an HDF5-based sparse matrix representation instead of
  the conventional (i.e. dense) HDF5 representation.

  If your dataset uses the conventional (i.e. dense) HDF5 representation,
  use the \code{HDF5Array()} constructor.

  If your dataset uses the HDF5-based sparse matrix representation from
  10x Genomics, use the \code{\link{TENxMatrix}()} constructor.
}

\seealso{
  \itemize{
    \item \link{TENxMatrix} objects for representing 10x Genomics
          datasets as \link[DelayedArray]{DelayedMatrix} objects.

    \item \link[DelayedArray]{DelayedArray} objects in the \pkg{DelayedArray}
          package.

    \item \code{\link{writeHDF5Array}} for writing an array-like object
          to an HDF5 file.

    \item \link{HDF5-dump-management} for controlling the location and
          physical properties of automatically created HDF5 datasets.

    \item \code{\link{saveHDF5SummarizedExperiment}} and
          \code{\link{loadHDF5SummarizedExperiment}} in this
          package (the \pkg{HDF5Array} package) for saving/loading
          an HDF5-based \link[SummarizedExperiment]{SummarizedExperiment}
          object to/from disk.

    \item The \link{HDF5ArraySeed} helper class.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## CONSTRUCTION
## ---------------------------------------------------------------------
library(h5vcData)
tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")

library(rhdf5)  # for h5ls()
h5ls(tally_file)

## Pick up "Coverages" dataset for Human chromosome 16:
cov0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")
cov0

is(cov0, "DelayedArray")  # TRUE
seed(cov0)
path(cov0)
chunkdim(cov0)

## ---------------------------------------------------------------------
## dim/dimnames
## ---------------------------------------------------------------------
dim(cov0)

dimnames(cov0)
dimnames(cov0) <- list(paste0("s", 1:6), c("+", "-"), NULL)
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
## SummarizedExperiment OBJECTS WITH DELAYED ASSAYS
## ---------------------------------------------------------------------

## DelayedArray objects can be used inside a SummarizedExperiment object
## to hold the assay data and to delay operations on them.
 
library(SummarizedExperiment)

pcov <- cov0[ , 1, ]  # coverage on plus strand
mcov <- cov0[ , 2, ]  # coverage on minus strand

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
\keyword{classes}
\keyword{methods}
