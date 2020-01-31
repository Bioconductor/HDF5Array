\name{writeHDF5Array}

\alias{class:HDF5RealizationSink}
\alias{HDF5RealizationSink-class}
\alias{HDF5RealizationSink}

\alias{dimnames,HDF5RealizationSink-method}
\alias{type,HDF5RealizationSink-method}
\alias{chunkdim,HDF5RealizationSink-method}
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
writeHDF5Array(x, filepath=NULL, name=NULL, chunkdim=NULL, level=NULL,
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
    dump file} i.e. to the file whose path is \code{\link{getHDF5DumpFile}}.
  }
  \item{name}{
    \code{NULL} or the name of the HDF5 dataset to write.
    If \code{NULL}, then the name returned by \code{\link{getHDF5DumpName}}
    will be used.
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

\seealso{
  \itemize{
    \item \link{HDF5Array} objects.

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
m <- matrix(runif(350, min=-1), nrow=25)
out_file <- tempfile()

M1 <- writeHDF5Array(m, out_file, name="M1", chunkdim=c(5, 5))
M1
chunkdim(M1)

## ---------------------------------------------------------------------
## WRITE A DelayedArray OBJECT TO AN HDF5 FILE
## ---------------------------------------------------------------------
M2 <- log(t(DelayedArray(m)) + 1)
M2 <- writeHDF5Array(M2, out_file, name="M2", chunkdim=c(5, 5))
M2
chunkdim(M2)

library(rhdf5)
library(h5vcData)

tally_file <- system.file("extdata", "example.tally.hfs5",
                          package="h5vcData")
h5ls(tally_file)

cvg0 <- HDF5Array(tally_file, "/ExampleStudy/16/Coverages")

cvg1 <- cvg0[ , , 29000001:29000007]

writeHDF5Array(cvg1, out_file, "cvg1")
h5ls(out_file)
}
\keyword{methods}
