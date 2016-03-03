\name{HDF5Matrix-class}
\docType{class}

\alias{class:HDF5Matrix}
\alias{HDF5Matrix-class}
\alias{HDF5Matrix}

\alias{dim,HDF5Matrix-method}
\alias{dimnames,HDF5Matrix-method}
\alias{dimnames<-,HDF5Matrix-method}

\alias{coerce,HDF5Array,HDF5Matrix-method}

\alias{show,HDF5Matrix-method}

\alias{[,HDF5Matrix-method}

\title{HDF5Matrix objects}

\description{
  The HDF5Matrix class is a container for representing and manipulating an
  on-disk HDF5 dataset with 2 dimensions, or a 2-dimensional slice of an
  on-disk HDF5 dataset with 2 or more dimensions.

  HDF5Matrix objects can be manipulated via a \emph{restricted matrix-like
  interface} i.e. they support common \link[base]{matrix} operations like
  \code{dim}, \code{nrow}, \code{ncol}, \code{dimnames}, \code{rownames},
  \code{colnames}, \code{[}, \code{rowSums}, \code{colSums}, etc... as long
  as the operation doesn't need to alter the content of the underlying
  (on-disk) dataset or write a new dataset to the disk.
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
