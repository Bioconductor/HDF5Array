\name{h5ls}

\alias{h5ls}

\title{A wrapper to rhdf5::h5ls() that works on H5File objects}

\description{
  Like \code{rhdf5::\link[rhdf5]{h5ls}()}, but works on an
  \link{H5File} object.
}

\usage{
h5ls(file, recursive=TRUE, all=FALSE, datasetinfo=TRUE,
           index_type=h5default("H5_INDEX"), order=h5default("H5_ITER"),
           s3=FALSE, s3credentials=NULL, native=FALSE)
}

\arguments{
  \item{file, recursive, all, datasetinfo, index_type, order, s3, s3credentials, native}{

    See \code{?rhdf5::\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package
    for a description of these arguments.

    Note that the only difference with \code{rhdf5::\link[rhdf5]{h5ls}()}
    is that, with \code{HDF5Array::h5ls()}, \code{file} can be an
    \link{H5File} object.
  }
}

\value{
  See \code{?rhdf5::\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
}

\seealso{
  \itemize{
    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item \link{H5File} objects.
  }
}

\examples{
toy_h5 <- system.file("extdata", "toy.h5", package="HDF5Array")
h5ls(toy_h5)

h5file <- H5File(toy_h5)
h5ls(h5file)

## See '?H5File' for more examples.
}
\keyword{utilities}
