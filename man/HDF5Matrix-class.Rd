\name{HDF5Matrix-class}
\docType{class}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\title{HDF5Matrix objects}

\description{
}

\usage{
HDF5Matrix(file, group, name)  # constructor function
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
  A HDF5Matrix object.
}

\section{Accessors}{
}

\section{Coercion}{
}

\section{Subsetting}{
}

\seealso{
  \itemize{
    \item \link[base]{matrix} objects in base R.

    \item \link{HDF5Array} objects.
  }
}

\examples{
library(rhdf5)
library(h5vcData)

tallyFile <- system.file("extdata", "example.tally.hfs5",
                         package="h5vcData")
h5ls(tallyFile)

cov0 <- HDF5Array(tallyFile, "/ExampleStudy/16", "Coverages")
cov1 <- cov0[ , , 29000001:29000007]
cov2 <- cov0[ , 1, 29000001:29000007]

pcov <- as(cov0[ , 1, ], "HDF5Matrix")  # coverage on plus strand
mcov <- as(cov0[ , 2, ], "HDF5Matrix")  # coverage on minus strand

library(SummarizedExperiment)
se <- SummarizedExperiment(list(pcov=pcov, mcov=mcov))
validObject(se, complete=TRUE)
}
\keyword{methods}
\keyword{classes}
