\name{HDF5Array-class}
\docType{class}

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\alias{dim,HDF5Array-method}
\alias{dimnames,HDF5Array-method}
\alias{dimnames<-,HDF5Array-method}

\alias{as.array,HDF5Array-method}
\alias{as.matrix,HDF5Array-method}

\alias{show,HDF5Array-method}

\alias{[,HDF5Array-method}

\title{HDF5Array objects}

\description{
  The HDF5Array class is a container for representing and manipulating an
  on-disk HDF5 dataset.

  HDF5Array objects can be manipulated via a \emph{restricted assay-like
  interface} interface i.e. they support common \link[base]{array} operations
  like \code{dim}, \code{dimnames}, \code{[}, etc... as long as the operation
  doesn't need to alter the content of the underlying (on-disk) dataset or
  write a new dataset to the disk.
}

\usage{
HDF5Array(file, group, name)  # constructor function
}

\arguments{
  \item{file}{
  }
  \item{group}{
  }
  \item{name}{
  }
}

\value{
  A HDF5Array object.
}

\section{Accessors}{
}

\section{Coercion}{
}

\section{Subsetting}{
}

\seealso{
  \itemize{
    \item \link[base]{array} objects in base R.

    \item \link{HDF5Matrix} objects.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which the
          \pkg{HDF5Array} package is building.
  }
}

\examples{
library(rhdf5)
library(h5vcData)

tallyFile <- system.file("extdata", "example.tally.hfs5",
                         package="h5vcData")
h5ls(tallyFile)

cov0 <- HDF5Array(tallyFile, "/ExampleStudy/16", "Coverages")
cov0

cov1 <- cov0[ , , 29000001:29000007]
cov1
as.array(cov1)

cov2 <- cov0[ , 1, 29000001:29000007]
cov2
as.matrix(cov2)
}
\keyword{methods}
\keyword{classes}
