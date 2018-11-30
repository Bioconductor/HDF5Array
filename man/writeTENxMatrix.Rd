\name{writeTENxMatrix}

\alias{class:TENxRealizationSink}
\alias{TENxRealizationSink-class}
\alias{TENxRealizationSink}

\alias{dimnames,TENxRealizationSink-method}
\alias{type,TENxRealizationSink-method}
\alias{write_sparse_block,TENxRealizationSink-method}
\alias{write_block,TENxRealizationSink-method}
\alias{close,TENxRealizationSink-method}

\alias{coerce,TENxRealizationSink,TENxMatrixSeed-method}
\alias{coerce,TENxRealizationSink,TENxMatrix-method}
\alias{coerce,TENxRealizationSink,DelayedArray-method}

\alias{writeTENxMatrix}

\alias{coerce,ANY,TENxMatrix-method}
\alias{coerce,DelayedArray,TENxMatrix-method}
\alias{coerce,DelayedMatrix,TENxMatrix-method}

\title{Write a matrix-like object as an HDF5-based sparse matrix}

\description{
  The 1.3 Million Brain Cell Dataset and other datasets published by
  10x Genomics use an HDF5-based sparse matrix representation instead
  of the conventional (i.e. dense) HDF5 representation.

  \code{writeTENxMatrix} writes a matrix-like object to this format.

  IMPORTANT NOTE: Only use \code{writeTENxMatrix} if the matrix-like
  object to write is sparse, that is, if most of its elements are zero.
  Using \code{writeTENxMatrix} on dense data is very inefficient!
  In this case, you should use \code{\link{writeHDF5Array}} instead.
}

\usage{
writeTENxMatrix(x, filepath=NULL, group=NULL, level=NULL, verbose=FALSE)
}

\arguments{
  \item{x}{
    The matrix-like object to write to an HDF5 file.

    The object to write should typically be sparse, that is, most of its
    elements should be zero.

    If \code{x} is a \link{DelayedArray} object, \code{writeTENxMatrix}
    \emph{realizes} it on disk, that is, all the delayed operations carried
    by the object are executed while the object is written to disk.
  }
  \item{filepath}{
    \code{NULL} or the path (as a single string) to the (new or existing)
    HDF5 file where to write the data.
    If \code{NULL}, then the data will be written to the current \emph{HDF5
    dump file} i.e. to the file whose path is \code{\link{getHDF5DumpFile}}.
  }
  \item{group}{
    \code{NULL} or the name of the HDF5 group where to write the data.
    If \code{NULL}, then the name returned by \code{\link{getHDF5DumpName}}
    will be used.
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
  and the performance of the disk, \code{writeTENxMatrix} can take a long
  time to complete. Use \code{verbose=TRUE} to see its progress.

  Use \code{\link{setHDF5DumpFile}} and \code{\link{setHDF5DumpName}} to
  control the location of automatically created HDF5 datasets.
}

\value{
  A \link{TENxMatrix} object pointing to the newly written HDF5 data on disk.
}

\seealso{
  \itemize{
    \item \link{TENxMatrix} objects.

    \item The \code{\link[TENxBrainData]{TENxBrainData}} dataset (in the
          \pkg{TENxBrainData} package).

    \item \link{HDF5-dump-management} to control the location and
          physical properties of automatically created HDF5 datasets.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A SIMPLE EXAMPLE
## ---------------------------------------------------------------------
m0 <- matrix(0L, nrow=25, ncol=12,
             dimnames=list(letters[1:25], LETTERS[1:12]))
m0[cbind(2:24, c(12:1, 2:12))] <- 100L + sample(55L, 23, replace=TRUE)
out_file <- tempfile()
M0 <- writeTENxMatrix(m0, out_file, group="m0")
M0
sparsity(M0)

path(M0)  # same as 'out_file'

## Use the h5ls() command from the rhdf5 package to see the structure of
## the file:
library(rhdf5)
h5ls(path(M0))

## ---------------------------------------------------------------------
## USING THE "1.3 Million Brain Cell Dataset"
## ---------------------------------------------------------------------

## The 1.3 Million Brain Cell Dataset from 10x Genomics is available via
## ExperimentHub:
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]
oneM <- TENxMatrix(fname, "mm10")  # see ?TENxMatrix for the details
oneM

## Note that the following transformation preserves sparsity:
M2 <- log(oneM + 1)  # delayed
M2                   # a DelayedMatrix instance

## In order to reduce computation times, we'll write only the first
## 5000 columns of M2 to disk:
out_file <- tempfile()
M3 <- writeTENxMatrix(M2[ , 1:5000], out_file, group="mm10", verbose=TRUE)
M3                   # a TENxMatrix instance
}
\keyword{methods}
