\name{DelayedArray-utils}

\alias{DelayedArray-utils}

\alias{length,ArrayBlocks-method}

% DelayedArray methods

\alias{is.na,DelayedArray-method}
\alias{!,DelayedArray-method}

\alias{+,DelayedArray,missing-method}
\alias{-,DelayedArray,missing-method}

\alias{pmax2}
\alias{pmin2}
\alias{pmax2,ANY,ANY-method}
\alias{pmin2,ANY,ANY-method}
\alias{pmax2,DelayedArray,vector-method}
\alias{pmin2,DelayedArray,vector-method}
\alias{pmax2,vector,DelayedArray-method}
\alias{pmin2,vector,DelayedArray-method}
\alias{pmax2,DelayedArray,DelayedArray-method}
\alias{pmin2,DelayedArray,DelayedArray-method}

\alias{anyNA,DelayedArray-method}

\alias{mean.DelayedArray}
\alias{mean,DelayedArray-method}

\alias{apply}
\alias{apply,DelayedArray-method}

% DelayedMatrix methods

\alias{rowSums,DelayedMatrix-method}
\alias{colSums,DelayedMatrix-method}
\alias{rowMeans,DelayedMatrix-method}
\alias{colMeans,DelayedMatrix-method}
\alias{\%*\%,DelayedMatrix,matrix-method}
\alias{\%*\%,matrix,DelayedMatrix-method}
\alias{\%*\%,DelayedMatrix,DelayedMatrix-method}

\title{Common operations on DelayedArray objects}

\description{
  Common operations on \link{DelayedArray} objects.
}

\details{
  Operations currently supported on all \link{DelayedArray} objects:
  \itemize{
    \item \code{is.na}, \code{!}
    \item The members of the \code{\link[methods]{Math}} group.
    \item The members of the \code{\link[methods]{Ops}} group.
    \item \code{pmax2}, \code{pmin2}.
    \item \code{anyNA}.
    \item The members of the \code{\link[methods]{Summary}} group.
    \item \code{mean}.
    \item \code{apply}.
  }

  Operations currently supported on \link{DelayedMatrix} objects only:
  \itemize{
    \item \code{rowSums}, \code{colSums}, \code{rowMeans}, \code{colMeans}.
    \item Matrix multiplication (\%*\%) of an ordinary matrix by a
          DelayedMatrix object.
  }
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{setHDF5OutputFile} to control the location of automatically
          created HDF5 datasets.

    \item \link{HDF5Array} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## Matrix multiplication writes a new HDF5 dataset to disk and returns
## an HDF5Matrix object that points to this new dataset.
m <- matrix(runif(60), ncol=12)
M <- DelayedArray(matrix(runif(240), nrow=12))

getHDF5OutputFile()  # HDF5 file where the new dataset will be written
lsHDF5OutputFile()
P <- m \%*\% M
P
getHDF5OutputFile()
lsHDF5OutputFile()
}
\keyword{methods}
