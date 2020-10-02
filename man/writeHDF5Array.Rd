\name{writeHDF5Array}

\alias{class:HDF5RealizationSink}
\alias{HDF5RealizationSink-class}
\alias{HDF5RealizationSink}

\alias{dimnames,HDF5RealizationSink-method}
\alias{type,HDF5RealizationSink-method}
\alias{chunkdim,HDF5RealizationSink-method}
\alias{is_sparse,HDF5RealizationSink-method}
\alias{write_block,HDF5RealizationSink-method}

\alias{coerce,HDF5RealizationSink,HDF5ArraySeed-method}
\alias{coerce,HDF5RealizationSink,HDF5Array-method}
\alias{coerce,HDF5RealizationSink,DelayedArray-method}

\alias{writeHDF5Array}

\alias{coerce,ANY,HDF5Array-method}
\alias{coerce,DelayedArray,HDF5Array-method}
\alias{coerce,DelayedMatrix,HDF5Matrix-method}

\title{Write an array-like object to an HDF5 file}

\description{
  A function for writing an array-like object to an HDF5 file.
}

\usage{
writeHDF5Array(x, filepath=NULL, name=NULL,
                  H5type=NULL, chunkdim=NULL, level=NULL, as.sparse=NA,
                  with.dimnames=FALSE, verbose=NA)
}

\arguments{
  \item{x}{
    The array-like object to write to an HDF5 file.

    If \code{x} is a \link{DelayedArray} object, \code{writeHDF5Array}
    \emph{realizes} it on disk, that is, all the delayed operations carried
    by the object are executed while the object is written to disk.
    See "On-disk realization of a DelayedArray object as an HDF5 dataset"
    section below for more information.
  }
  \item{filepath}{
    \code{NULL} or the path (as a single string) to the (new or existing)
    HDF5 file where to write the dataset.
    If \code{NULL}, then the dataset will be written to the current \emph{HDF5
    dump file} i.e. to the file whose path is \code{\link{getHDF5DumpFile}}.
  }
  \item{name}{
    \code{NULL} or the name of the HDF5 dataset to write.
    If \code{NULL}, then the name returned by \code{\link{getHDF5DumpName}}
    will be used.
  }
  \item{H5type}{
    The H5 datatype to use for the HDF5 dataset to be written to the HDF5
    file is automatically inferred from the type of \code{x} (\code{type(x)}).
    Advanced users can override this by specifying the H5 datatype they want
    via the \code{H5type} argument.

    See \code{rhdf5::h5const("H5T")} for a list of available H5 datatypes.
    See References section below for the link to the HDF Group's Support
    Portal where H5 predefined datatypes are documented.

    A typical use case is to use a datatype that is smaller than the
    automatic one in order to reduce the size of the dataset on disk.
    For example you could use \code{"H5T_IEEE_F32LE"} when \code{type(x)}
    is \code{"double"} and you don't care about preserving the precision
    of 64-bit floating-point numbers (the automatic H5 datatype used for
    \code{"double"} is \code{"H5T_IEEE_F64LE"}).
    Another example is to use \code{"H5T_STD_U16LE"} when \code{x} contains
    small non-negative integer values like counts (the automatic H5
    datatype used for \code{"integer"} is \code{"H5T_STD_I32LE"}).
  }
  \item{chunkdim}{
    The dimensions of the chunks to use for writing the data to disk.
    By default (i.e. when \code{chunkdim} is set to \code{NULL}),
    \code{getHDF5DumpChunkDim(dim(x))} will be used.
    See \code{?\link{getHDF5DumpChunkDim}} for more information.

    Set \code{chunkdim} to 0 to write \emph{unchunked data} (a.k.a.
    \emph{contiguous data}).
  }
  \item{level}{
    The compression level to use for writing the data to disk.
    By default, \code{getHDF5DumpCompressionLevel()} will be used.
    See \code{?\link{getHDF5DumpCompressionLevel}} for more information.
  }
  \item{as.sparse}{
    Whether the data in the returned \link{HDF5Array} object should be
    flagged as sparse or not. If set to \code{NA} (the default), then
    \code{is_sparse(x)} is used.

    IMPORTANT NOTE: This only controls the \code{as.sparse} flag of
    the returned \link{HDF5Array} object. See man page of the
    \code{\link{HDF5Array}()} constructor for more information.
    In particular this does NOT affect how the data will be laid out
    to the HDF5 file in any way (HDF5 doesn't natively support sparse
    storage at the moment). In other words, the data will always be
    stored in a dense format, even when \code{as.sparse} is set to
    \code{TRUE}.
  }
  \item{with.dimnames}{
    By default the dimnames on \code{x} are not written to the HDF5 file.
    Set \code{with.dimnames} to \code{TRUE} to also have them written.

    Note that \code{\link{h5writeDimnames}} is used internally to write
    the dimnames to disk. Setting \code{with.dimnames} to \code{FALSE}
    and calling \code{\link{h5writeDimnames}} is another way to write
    the dimnames on \code{x} to disk that gives more control.
    See \code{?\link{h5writeDimnames}} for more information.
  }
  \item{verbose}{
    Whether block processing progress should be displayed or not.
    If set to \code{NA} (the default), verbosity is controlled
    by \code{DelayedArray:::get_verbose_block_processing()}.
    Setting \code{verbose} to \code{TRUE} or \code{FALSE} overrides this.
  }
}

\details{
  Please note that, depending on the size of the data to write to disk
  and the performance of the disk, \code{writeHDF5Array()} can take a
  long time to complete. Use \code{verbose=TRUE} to see its progress.

  Use \code{\link{setHDF5DumpFile}} and \code{\link{setHDF5DumpName}} to
  control the location of automatically created HDF5 datasets.

  Use \code{\link{setHDF5DumpChunkLength}},
  \code{\link{setHDF5DumpChunkShape}}, and
  \code{\link{setHDF5DumpCompressionLevel}}, to control the
  physical properties of automatically created HDF5 datasets.
}

\value{
  An \link{HDF5Array} object pointing to the newly written HDF5 dataset
  on disk.
}

\section{On-disk realization of a DelayedArray object as an HDF5 dataset}{
  When passed a \link{DelayedArray} object, \code{writeHDF5Array}
  \emph{realizes} it on disk, that is, all the delayed operations carried
  by the object are executed on-the-fly while the object is written to disk.
  This uses a block-processing strategy so that the full object is not
  realized at once in memory. Instead the object is processed block by block
  i.e. the blocks are realized in memory and written to disk one at a time.

  In other words, \code{writeHDF5Array(x, ...)} is semantically equivalent
  to \code{writeHDF5Array(as.array(x), ...)}, except that \code{as.array(x)}
  is not called because this would realize the full object at once in memory.

  See \code{?\link{DelayedArray}} for general information about
  \link{DelayedArray} objects.
}

\references{
  Documentation of the H5 predefined datatypes on the HDF Group's Support
  Portal: \url{https://portal.hdfgroup.org/display/HDF5/Predefined+Datatypes}
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects.

    \item \code{\link{h5writeDimnames}} for writing the dimnames of an
          HDF5 dataset to disk.

    \item \code{\link{saveHDF5SummarizedExperiment}} and
          \code{\link{loadHDF5SummarizedExperiment}} in this
          package (the \pkg{HDF5Array} package) for saving/loading
          an HDF5-based \link[SummarizedExperiment]{SummarizedExperiment}
          object to/from disk.

    \item \link{HDF5-dump-management} to control the location and
          physical properties of automatically created HDF5 datasets.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## WRITE AN ORDINARY ARRAY TO AN HDF5 FILE
## ---------------------------------------------------------------------
m <- matrix(runif(364, min=-1), nrow=26,
            dimnames=list(letters, LETTERS[1:14]))

h5file <- tempfile(fileext=".h5")
M1 <- writeHDF5Array(m, h5file, name="M1", chunkdim=c(5, 5))
M1
chunkdim(M1)

## By default, writeHDF5Array() does not write the dimnames to the HDF5
## file so they are lost:
dimnames(M1)   # no dimnames

## Set 'with.dimnames' to TRUE to write them to the file:
M1b <- writeHDF5Array(m, h5file, name="M1b", with.dimnames=TRUE)
dimnames(M1b)  # dimnames are back

## With sparse data:
sm <- rsparsematrix(20, 8, density=0.1)
M2 <- writeHDF5Array(sm, h5file, name="M2", chunkdim=c(5, 5))
M2
is_sparse(M2)  # TRUE

## ---------------------------------------------------------------------
## WRITE A DelayedArray OBJECT TO AN HDF5 FILE
## ---------------------------------------------------------------------
M3 <- log(t(DelayedArray(m)) + 1)
M3 <- writeHDF5Array(M3, h5file, name="M3", chunkdim=c(5, 5))
M3
chunkdim(M3)

library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

cvg0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")

cvg1 <- cvg0[ , , 29000001:29000007]

writeHDF5Array(cvg1, h5file, "cvg1")
h5ls(h5file)
}
\keyword{methods}
