\name{HDF5Array-class}
\docType{class}

\alias{class:HDF5Array}
\alias{HDF5Array-class}
\alias{HDF5Array}

\title{HDF5Array objects}

\description{
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
