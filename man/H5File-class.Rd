\name{H5File-class}
\docType{class}

\alias{class:H5FileID}
\alias{H5FileID-class}
\alias{H5FileID}

\alias{open.H5FileID}
\alias{close.H5FileID}
\alias{show,H5FileID-method}

\alias{class:H5File}
\alias{H5File-class}
\alias{H5File}

\alias{open.H5File}
\alias{close.H5File}

\alias{class:character_OR_H5File}
\alias{character_OR_H5File-class}
\alias{character_OR_H5File}

\alias{class:H5DSetDescriptor}
\alias{H5DSetDescriptor-class}
\alias{H5DSetDescriptor}
\alias{destroy_H5DSetDescriptor}
\alias{show,H5DSetDescriptor-method}

\title{H5File objects}

\description{
  *** WORK-IN-PROGRESS! ***

  The H5File class provides a formal representation of an HDF5 file.
}

\usage{
## Constructor function:
H5File(filepath, s3=FALSE, s3credentials=NULL)
}

\arguments{
  \item{filepath}{
    A single string specifying the path or URL to an HDF5 file.
  }
  \item{s3}{
    *** NOT SUPPORTED YET! ***

    \code{TRUE} or \code{FALSE}. Should the file argument be
    treated as a URL to an Amazon S3 bucket, rather than a local file path?
  }
  \item{s3credentials}{
    *** NOT SUPPORTED YET! ***

    A list of length three, providing the credentials for accessing files
    in a private Amazon S3 bucket.
  }
}

\details{
  COMING SOON...
}

\value{
  An H5File object.
}

\seealso{
  \itemize{
    \item The \link{HDF5Array} class for representing and operating on
          a conventional (a.k.a. dense) HDF5 dataset.

    \item The \link{TENxMatrix} class for representing and operating on
          a 10x Genomics dataset.

    \item The \link{H5ADMatrix} class for representing and operating on
          the matrix stored in an \code{.h5ad} file.

    \item The \code{\link{h5mread}} function in this package (\pkg{HDF5Array})
          that is used internally by \link{HDF5Array}, \link{TENxMatrix},
          and \link{H5ADMatrix} objects, for (almost) all their data
          reading needs.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
library(rhdf5)  # for h5ls()
h5ls(toy_h5)

h5file <- H5File(toy_h5)
h5file

h5mread(h5file, "M2", list(1:10, 1:6))
get_h5mread_returned_type(h5file, "M2")
}
\keyword{classes}
\keyword{methods}
