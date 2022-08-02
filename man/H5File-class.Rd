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
  The H5File class provides a formal representation of an HDF5
  file (local or remote).
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
    See \code{?\link[rhdf5]{H5Pset_fapl_ros3}} in the \pkg{rhdf5} package
    for more information.
  }
  \item{.no_rhdf5_h5id}{
    For internal use only. Don't use.
  }
}

\details{
  IMPORTANT NOTE ABOUT H5File OBJECTS AND PARALLEL EVALUATION

  The short story is that H5File objects cannot be used in the context
  of parallel evaluation at the moment.

  Here is why:

  H5File objects contain an identifier to an open connection to the HDF5
  file. This identifier becomes invalid in the 2 following situations:
  \itemize{
    \item After serialization/deserialization, that is, after loading
          a serialized H5File object with \code{\link{readRDS}()} or
          \code{\link{load}()}.
    \item In the context of parallel evaluation, when using the
          \link[BiocParallel]{SnowParam} parallelization backend.
          This is because, unlike the \link[BiocParallel]{MulticoreParam}
          backend which used a system \code{fork}, the
          \link[BiocParallel]{SnowParam} backend uses
          serialization/deserialization to transmit the object
          to the workers.
  }
  In both cases, the connection to the file is lost and any attempt to
  read data from the H5File object will fail. Note that the above also
  happens to any H5File object that got serialized indirectly i.e. as
  part of a bigger object. For example, if an \link{HDF5Array} object
  was constructed from an H5File object, then it contains the H5File
  object and therefore
  \code{\link[DelayedArray]{blockApply}(..., BPPARAM=SnowParam(4))}
  cannot be used on it.

  Furthermore, even if sometimes an H5File object \emph{seems} to work
  fine with the \link[BiocParallel]{MulticoreParam} parallelization
  backend, this is highly unreliable and must be avoided.
}

\value{
  An H5File object.
}

\seealso{
  \itemize{
    \item \link[rhdf5]{H5Pset_fapl_ros3} in the \pkg{rhdf5} package for
          detailed information about how to pass your S3 credentials
          to the \code{s3credentials} argument.

    \item The \link{HDF5Array} class for representing and operating on
          a conventional (a.k.a. dense) HDF5 dataset.

    \item The \link{H5SparseMatrix} class for representing and operating on
          an HDF5 sparse matrix.

    \item The \link{H5ADMatrix} class for representing and operating on
          the central matrix of an \code{h5ad} file, or any matrix in
          its \code{/layers} group.

    \item The \link{TENxMatrix} class for representing and operating on
          a 10x Genomics dataset.

    \item The \code{\link{h5mread}} function in this package (\pkg{HDF5Array})
          that is used internally by \link{HDF5Array}, \link{TENxMatrix},
          and \link{H5ADMatrix} objects, for (almost) all their data
          reading needs.

    \item \code{\link{h5ls}} to list the content of an HDF5 file.

    \item \code{\link[BiocParallel]{bplapply}},
          \code{\link[BiocParallel]{MulticoreParam}}, and
          \code{\link[BiocParallel]{SnowParam}}, in the
          \pkg{BiocParallel} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. BASIC USAGE
## ---------------------------------------------------------------------

## With a local file:
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5file1 <- H5File(toy_h5)
h5ls(h5file1)
path(h5file1)

h5mread(h5file1, "M2", list(1:10, 1:6))
get_h5mread_returned_type(h5file1, "M2")

## With a file stored in an Amazon S3 bucket:
if (Sys.info()[["sysname"]] != "Darwin") {
  public_S3_url <-
   "https://rhdf5-public.s3.eu-central-1.amazonaws.com/rhdf5ex_t_float_3d.h5"
  h5file2 <- H5File(public_S3_url, s3=TRUE)
  h5ls(h5file2)

  h5mread(h5file2, "a1")
  get_h5mread_returned_type(h5file2, "a1")
}

## ---------------------------------------------------------------------
## B. H5File OBJECTS AND PARALLEL EVALUATION
## ---------------------------------------------------------------------
## H5File objects cannot be used in the context of parallel evaluation
## at the moment!

library(BiocParallel)

FUN1 <- function(i, h5file, name)
    sum(HDF5Array::h5mread(h5file, name, list(i, NULL)))

FUN2 <- function(i, h5file, name)
    sum(HDF5Array::h5mread(h5file, name, list(i, NULL, NULL)))

## With the SnowParam parallelization backend, the H5File object
## does NOT work on the workers:
\dontrun{
## ERROR!
res1 <- bplapply(1:150, FUN1, h5file1, "M2", BPPARAM=SnowParam(3))
## ERROR!
res2 <- bplapply(1:5, FUN2, h5file2, "a1", BPPARAM=SnowParam(3))
}

## With the MulticoreParam parallelization backend, the H5File object
## might seem to work on the workers. However this is highly unreliable
## and must be avoided:
\dontrun{
if (.Platform$OS.type != "windows") {
  ## UNRELIABLE!
  res1 <- bplapply(1:150, FUN1, h5file1, "M2", BPPARAM=MulticoreParam(3))
  ## UNRELIABLE!
  res2 <- bplapply(1:5, FUN2, h5file2, "a1", BPPARAM=MulticoreParam(3))
}
}
}
\keyword{classes}
\keyword{methods}
