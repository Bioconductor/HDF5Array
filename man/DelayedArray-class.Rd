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
\alias{length,DelayedArray-method}
\alias{isEmpty,DelayedArray-method}
\alias{dimnames,DelayedArray-method}
\alias{dimnames<-,DelayedArray-method}
\alias{names,DelayedArray-method}
\alias{names<-,DelayedArray-method}

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
  DelayedArray object allows one to perform common array operations in an
  efficient way. In order to optimize performance and reduce memory usage,
  operations on DelayedArray objects are implemented as delayed operations
  and/or via a block processing strategy.
}

\usage{
DelayedArray(x)
}

\arguments{
  \item{x}{
    An array-like object.
  }
}

\section{Coercion}{
  To \emph{realize} a DelayedArray object (i.e. to trigger execution of the
  delayed operations that it carries and turn the object back into an ordinary
  array), call \code{as.array} on it. Alternatively, write it to disk by
  coercing it to HDF5Array (other on-disk backends can be supported). In that
  case, the full object is not realized at once in memory, but split into
  blocks first, and blocks are realized one at a time and written to disk.
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
          object with the same number of dimensions.
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

    \item \link{HDF5_output_settings} to manage settings for writing
          output to HDF5.

    \item \link{HDF5Array} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## WITH AN ORDINARY array OBJECT
## ---------------------------------------------------------------------

a <- array(runif(15000000), c(300000, 10, 5))
A <- DelayedArray(a)
A

toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- toto(a)
head(b)

B <- toto(A)  # very fast! (operations are delayed)
B             # still 3 dimensions (subsetting never drops dimensions)
B <- drop(B)
B

cs <- colSums(b)
CS <- colSums(B)
stopifnot(identical(cs, CS))

## ---------------------------------------------------------------------
## WITH A HDF5Dataset OBJECT
## ---------------------------------------------------------------------

a2 <- HDF5Dataset(a)    # create the dataset
a2

A2 <- DelayedArray(a2)  # wrap the dataset in a DelayedArray object
A2

B2 <- toto(A2)  # very fast! (operations are delayed)
B2 <- drop(B2)

CS2 <- colSums(B2)
stopifnot(identical(cs, CS2))
}
\keyword{methods}
\keyword{classes}
