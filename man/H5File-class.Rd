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

\alias{path,H5File-method}
\alias{open.H5File}
\alias{close.H5File}
\alias{show,H5File-method}
\alias{coerce,H5File,H5IdComponent-method}

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
  The H5File class provides a formal representation of an HDF5 file.
}

\usage{
## Constructor function:
H5File(filepath, s3=FALSE, s3credentials=NULL, .no_rhdf5_h5id=FALSE)
}

\arguments{
  \item{filepath}{
    A single string specifying the path or URL to an HDF5 file.
  }
  \item{s3}{
    \code{TRUE} or \code{FALSE}. Should the \code{filepath} argument be
    treated as the URL to a file stored in an Amazon S3 bucket, rather than
    the path to a local file?
  }
  \item{s3credentials}{
    A list of length 3, providing the credentials for accessing
    files stored in a private Amazon S3 bucket.
    See \code{?\link[rhdf5]{H5P_FILE_ACCESS}} in the \pkg{rhdf5} package
    for more information.
  }
  \item{.no_rhdf5_h5id}{
    For internal use.
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
    \item \link[rhdf5]{H5P_FILE_ACCESS} in the \pkg{rhdf5} package for
          detailed information about how to pass your S3 credentials
          to the \code{s3credentials} argument.

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

    \item \code{\link{h5ls}} to list the content of an HDF5 file.
  }
}

\examples{
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5file <- H5File(toy_h5)
h5ls(h5file)
path(h5file)

h5mread(h5file, "M2", list(1:10, 1:6))
get_h5mread_returned_type(h5file, "M2")

## With an URL to a file stored in an Amazon S3 bucket:
public_S3_url <-
 "https://rhdf5-public.s3.eu-central-1.amazonaws.com/rhdf5ex_t_float_3d.h5"
h5file2 <- H5File(public_S3_url, s3=TRUE)
h5ls(h5file2)

h5mread(h5file2, "a1")
get_h5mread_returned_type(h5file2, "a1")
}
\keyword{classes}
\keyword{methods}
