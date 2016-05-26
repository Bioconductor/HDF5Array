\name{DelayedArray-class}
\docType{class}

\alias{class:DelayedArray}
\alias{DelayedArray-class}
\alias{DelayedArray}

\alias{class:DelayedMatrix}
\alias{DelayedMatrix-class}
\alias{DelayedMatrix}

\alias{coerce,DelayedArray,DelayedMatrix-method}

\alias{dim,DelayedArray-method}
\alias{dim<-,DelayedArray-method}
\alias{length,DelayedArray-method}
\alias{isEmpty,DelayedArray-method}
\alias{dimnames,DelayedArray-method}
\alias{dimnames<-,DelayedArray-method}
\alias{names,DelayedArray-method}
\alias{names<-,DelayedArray-method}

\alias{drop,DelayedArray-method}

\alias{[,DelayedArray-method}
\alias{t,DelayedArray-method}

\alias{as.array.DelayedArray}
\alias{as.array,DelayedArray-method}
\alias{as.vector.DelayedArray}
\alias{as.vector,DelayedArray-method}
\alias{as.matrix.DelayedArray}
\alias{as.matrix,DelayedArray-method}

\alias{[[,DelayedArray-method}

\alias{show,DelayedArray-method}

\alias{c,DelayedArray-method}

\title{DelayedArray objects}

\description{
  Wrapping an array-like object (typically an on-disk object) in a
  DelayedArray object allows one to perform common array operations on it
  without loading the object in memory. In order to reduce memory usage and
  optimize performance, operations on the object are either delayed or
  executed using a block processing mechanism.
}

\usage{
DelayedArray(x)  # constructor function
}

\arguments{
  \item{x}{
    An array-like object.
  }
}

\section{Details}{
  To \emph{realize} a DelayedArray object (i.e. to trigger execution of the
  delayed operations carried by the object and return the result as an
  ordinary array), call \code{as.array} on it. However this realizes the
  object in memory and could require too much memory. Big DelayedArray
  objects are preferrably realized on disk e.g. by calling the
  \code{\link{HDF5Dataset}} constructor on it (other on-disk backends can be
  supported). In that case, the full object is not realized at once in memory,
  but split into small blocks first, and the blocks are realized and written
  to disk one at a time.
}

\section{Accessors}{
  DelayedArray objects support the same set of getters as ordinary arrays
  i.e. \code{dim()}, \code{length()}, and \code{dimnames()}.

  Only \code{dimnames()} is supported as a setter.
}

\section{Subsetting}{
  A DelayedArray object can be subsetted like an ordinary object but with
  the following differences:
  \itemize{
    \item The \code{drop} argument of the \code{[} operator is ignored i.e.
          subsetting a DelayedArray object always returns a DelayedArray
          object with the same number of dimensions. You need to call
          \code{drop()} on the subsetted object to actually drop its
          ineffective dimensions (i.e. the dimensions equal to 1).

    \item Linear subsetting (a.k.a. 1D-style subsetting, that is, subsetting
          with a single subscript \code{i}) is not supported.
  }

  Subsetting with \code{[[} is supported but only the linear form of it.

  DelayedArray objects don't support subassignment (\code{[<-} or \code{[[<-}).
}

\seealso{
  \itemize{
    \item \link{DelayedArray-utils} for common operations on DelayedArray
          objects.

    \item \link{setHDF5DumpFile} to control the location of automatically
          created HDF5 datasets.

    \item \link{HDF5Array} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## WITH AN ORDINARY array OBJECT
## ---------------------------------------------------------------------
a <- array(runif(1500000), c(10000, 30, 5))
A <- DelayedArray(a)
A

toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- toto(a)
head(b)

B <- toto(A)  # very fast! (operations are delayed)
B             # still 3 dimensions (subsetting a DelayedArray object
              # never drops dimensions)
B <- drop(B)
B

cs <- colSums(b)
CS <- colSums(B)
stopifnot(identical(cs, CS))

## ---------------------------------------------------------------------
## WITH A HDF5Dataset OBJECT
## ---------------------------------------------------------------------
h5a <- HDF5Dataset(a)    # create the dataset
h5a

A2 <- DelayedArray(h5a)  # wrap the dataset in a DelayedArray object
A2

B2 <- toto(A2)  # very fast! (operations are delayed)
B2 <- drop(B2)

CS2 <- colSums(B2)
stopifnot(identical(cs, CS2))

## ---------------------------------------------------------------------
## STORE THE RESULT IN A NEW HDF5Dataset OBJECT
## ---------------------------------------------------------------------
b2 <- HDF5Dataset(B2)  # "realize" B2 on disk (as an HDF5 dataset)

## If this is just an intermediate result, you can either keep going
## with B2 or replace it with b2 wrapped in a DelayedArray object etc...
B2 <- DelayedArray(b2)  # semantically equivalent to the previous B2
}
\keyword{methods}
\keyword{classes}
