\name{HDF5Array-class}
\docType{class}

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\alias{DelayedArray,HDF5ArraySeed-method}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{is_sparse<-,HDF5Array-method}

\alias{matrixClass,HDF5Array-method}
\alias{coerce,HDF5Array,HDF5Matrix-method}
\alias{coerce,HDF5Matrix,HDF5Array-method}
\alias{coerce,ANY,HDF5Matrix-method}

\title{HDF5 datasets as DelayedArray objects}

\description{
  The HDF5Array class is a \link[DelayedArray]{DelayedArray} subclass
  for representing and operating on a conventional (a.k.a. dense) HDF5
  dataset.

  All the operations available for \link[DelayedArray]{DelayedArray}
  objects work on HDF5Array objects.
}

\usage{
## Constructor function:
HDF5Array(filepath, name, as.sparse=FALSE, type=NA)
}

\arguments{
  \item{filepath}{
    The path (as a single string or \link{H5File} object) to the HDF5 file
    (\code{.h5} or \code{.h5ad}) where the dataset is located.

    Note that you must create and use an \link{H5File} object if the HDF5
    file to access is stored in an Amazon S3 bucket. See \code{?\link{H5File}}
    for how to do this.

    Also please note that \link{H5File} objects must NOT be used in the
    context of parallel evaluation at the moment.
  }
  \item{name}{
    The name of the dataset in the HDF5 file.
  }
  \item{as.sparse}{
    Whether the HDF5 dataset should be flagged as sparse or not, that is,
    whether it should be considered sparse (and treated as such) or not.
    Note that HDF5 doesn't natively support sparse storage at the moment
    so HDF5 datasets cannot be stored in a sparse format, only in a dense
    one. However a dataset stored in a dense format can still contain a lot
    of zeros. Using \code{as.sparse=TRUE} on such dataset will enable
    some optimizations that can lead to a lower memory footprint (and
    possibly better performance) when operating on the HDF5Array.

    IMPORTANT NOTE: If the dataset is in the 10x Genomics format (i.e. if
    it uses the HDF5-based sparse matrix representation from 10x Genomics),
    you should use the \code{\link{TENxMatrix}()} constructor instead of
    the \code{HDF5Array()} constructor.
  }
  \item{type}{
    By default the \code{\link[DelayedArray]{type}} of the returned
    object is inferred from the H5 datatype of the HDF5 dataset.
    This can be overridden by specifying the \code{type} argument.
    The specified type must be an \emph{R atomic type} (e.g.
    \code{"integer"}) or \code{"list"}.
  }
}

\value{
  An HDF5Array (or HDF5Matrix) object. (Note that HDF5Matrix extends HDF5Array.)
}

\note{
  The "1.3 Million Brain Cell Dataset" and other datasets published by
  10x Genomics use an HDF5-based sparse matrix representation instead
  of the conventional (a.k.a. dense) HDF5 representation.

  If your dataset uses the conventional (a.k.a. dense) HDF5 representation,
  use the \code{HDF5Array()} constructor documented here.

  But if your dataset uses the HDF5 sparse matrix representation from
  10x Genomics, use the \code{\link{TENxMatrix}()} constructor instead.
}

\seealso{
  \itemize{
    \item \link{H5File} objects.

    \item \link{H5SparseMatrix} objects for representing HDF5 sparse matrices
          as \link[DelayedArray]{DelayedMatrix} objects.

    \item \link{H5ADMatrix} objects for representing h5ad central
          matrices (or matrices in the \code{/layers} group)
          as \link[DelayedArray]{DelayedMatrix} objects.

    \item \link{TENxMatrix} objects for representing 10x Genomics
          datasets as \link[DelayedArray]{DelayedMatrix} objects.

    \item \link{ReshapedHDF5Array} objects for representing HDF5 datasets
          as \link[DelayedArray]{DelayedArray} objects with a user-supplied
          upfront virtual reshaping.

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

    \item \code{\link{h5ls}} to list the content of an HDF5 file (\code{.h5}
          or \code{.h5ad}).
  }
}

\examples{
## ---------------------------------------------------------------------
## A. CONSTRUCTION
## ---------------------------------------------------------------------

## With a local file:
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

HDF5Array(toy_h5, "M2")
HDF5Array(toy_h5, "M2", type="integer")
HDF5Array(toy_h5, "M2", type="complex")

## With a file stored in an Amazon S3 bucket:
if (Sys.info()[["sysname"]] != "Darwin") {
    public_S3_url <-
     "https://rhdf5-public.s3.eu-central-1.amazonaws.com/rhdf5ex_t_float_3d.h5"
    h5file <- H5File(public_S3_url, s3=TRUE)
    h5ls(h5file)

    HDF5Array(h5file, "a1")
}

## ---------------------------------------------------------------------
## B. BASIC MANIPULATION
## ---------------------------------------------------------------------

library(h5vcData)
tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

## Pick up "Coverages" dataset for Human chromosome 16:
name <- "/ExampleStudy/16/Coverages"
cvg <- HDF5Array(tally_file, name)
cvg

is(cvg, "DelayedArray")  # TRUE
seed(cvg)
path(cvg)
chunkdim(cvg)

## The data in the dataset looks sparse. In this case it is recommended
## to set 'as.sparse' to TRUE when constructing the HDF5Array object.
## This will make block processing (used in operations like sum()) more
## memory efficient and likely faster:
cvg0 <- HDF5Array(tally_file, name, as.sparse=TRUE)
is_sparse(cvg0)  # TRUE

## Note that we can also flag the HDF5Array object as sparse after
## creation:
is_sparse(cvg) <- TRUE
cvg  # same as 'cvg0'

## dim/dimnames:

dim(cvg0)

dimnames(cvg0)
dimnames(cvg0) <- list(paste0("s", 1:6), c("+", "-"), NULL)
dimnames(cvg0)

## ---------------------------------------------------------------------
## C. SLICING (A.K.A. SUBSETTING)
## ---------------------------------------------------------------------

cvg1 <- cvg0[ , , 29000001:29000007]
cvg1

dim(cvg1)
as.array(cvg1)
stopifnot(identical(dim(as.array(cvg1)), dim(cvg1)))
stopifnot(identical(dimnames(as.array(cvg1)), dimnames(cvg1)))

cvg2 <- cvg0[ , "+", 29000001:29000007]
cvg2
as.matrix(cvg2)

## ---------------------------------------------------------------------
## D. SummarizedExperiment OBJECTS WITH DELAYED ASSAYS
## ---------------------------------------------------------------------

## DelayedArray objects can be used inside a SummarizedExperiment object
## to hold the assay data and to delay operations on them.
 
library(SummarizedExperiment)

pcvg <- cvg0[ , 1, ]  # coverage on plus strand
mcvg <- cvg0[ , 2, ]  # coverage on minus strand

nrow(pcvg)  # nb of samples
ncol(pcvg)  # length of Human chromosome 16

## The convention for a SummarizedExperiment object is to have 1 column
## per sample so first we need to transpose 'pcvg' and 'mcvg':
pcvg <- t(pcvg)
mcvg <- t(mcvg)
se <- SummarizedExperiment(list(pcvg=pcvg, mcvg=mcvg))
se
stopifnot(validObject(se, complete=TRUE))

## A GPos object can be used to represent the genomic positions along
## the dataset:
gpos <- GPos(GRanges("16", IRanges(1, nrow(se))))
gpos
rowRanges(se) <- gpos
se
stopifnot(validObject(se))
assays(se)$pcvg
assays(se)$mcvg
}
\keyword{classes}
\keyword{methods}
