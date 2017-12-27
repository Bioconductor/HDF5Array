\name{writeHDF5Array}

\alias{class:HDF5RealizationSink}
\alias{HDF5RealizationSink-class}
\alias{HDF5RealizationSink}

\alias{dimnames,HDF5RealizationSink-method}
\alias{chunk_dim,HDF5RealizationSink-method}
\alias{write_block_to_sink,HDF5RealizationSink-method}

\alias{coerce,HDF5RealizationSink,HDF5ArraySeed-method}
\alias{coerce,HDF5RealizationSink,HDF5Array-method}
\alias{coerce,HDF5RealizationSink,DelayedArray-method}

\alias{writeHDF5Array}
\alias{writeHDF5Dataset}

\alias{coerce,ANY,HDF5Array-method}
\alias{coerce,DelayedArray,HDF5Array-method}
\alias{coerce,DelayedMatrix,HDF5Matrix-method}

\title{Write an array-like object to an HDF5 file}

\description{
  A function for writing an array-like object to an HDF5 file.
}

\usage{
writeHDF5Array(x, filepath=NULL, name=NULL, chunk_dim=NULL, level=NULL,
               verbose=FALSE)
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
    dump file} i.e. the path returned by \code{\link{getHDF5DumpFile}}
    will be used.
  }
  \item{name}{
    \code{NULL} or the name of the HDF5 dataset to write.
    If \code{NULL}, then the name returned by \code{\link{getHDF5DumpName}}
    will be used.
  }
  \item{chunk_dim}{
    The dimensions of the chunks to use for writing the data to disk.
    By default, \code{getHDF5DumpChunkDim(dim(x), type(x))} is used.
    See \code{?\link{getHDF5DumpChunkDim}} for more information.
  }
  \item{level}{
    The compression level to use for writing the data to disk.
    By default, \code{getHDF5DumpCompressionLevel()} is used.
    See \code{?\link{getHDF5DumpCompressionLevel}} for more information.
  }
  \item{verbose}{
    Set to \code{TRUE} to make the function display progress.
  }
}

\details{
  Please note that, depending on the size of the data to write to disk
  and the performance of the disk, \code{writeHDF5Array} can take a long
  time to complete. Use \code{verbose=TRUE} to see its progress.

  Use \code{\link{setHDF5DumpFile}} and \code{\link{setHDF5DumpName}} to
  control the location of automatically created HDF5 datasets.
}

\value{
  An HDF5Array object pointing to the newly written HDF5 dataset on disk.
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

\seealso{
  \itemize{
    \item \code{\link[SummarizedExperiment]{saveHDF5SummarizedExperiment}} and
          \code{\link[SummarizedExperiment]{loadHDF5SummarizedExperiment}} in
          the \pkg{SummarizedExperiment} package for saving/loading a
          HDF5-based SummarizedExperiment object to/from disk.

    \item \link{HDF5-dump-management} for controlling the location of
          automatically created HDF5 datasets.

    \item \link{HDF5Array} objects.

    \item \link{DelayedArray} objects.

    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.
  }
}

\examples{
library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

cov0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")

cov1 <- drop(cov0[ , , 29000001:29000007])

out_file <- tempfile()
writeHDF5Array(cov1, out_file, "cov1")
h5ls(out_file)
}
\keyword{methods}
