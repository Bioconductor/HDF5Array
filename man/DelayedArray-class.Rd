\name{DelayedArray-class}
\docType{class}

% DelayedArray
\alias{class:DelayedArray}
\alias{DelayedArray-class}

\alias{DelayedArray}
\alias{coerce,ANY,DelayedArray-method}

\alias{dim,DelayedArray-method}
\alias{length,DelayedArray-method}
\alias{isEmpty,DelayedArray-method}
\alias{dimnames,DelayedArray-method}

\alias{dimnames<-,DelayedArray-method}
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

% DelayedMatrix

\alias{class:DelayedMatrix}
\alias{DelayedMatrix-class}
\alias{DelayedMatrix}

\alias{coerce,DelayedArray,DelayedMatrix-method}
\alias{coerce,array,DelayedMatrix-method}

\title{DelayedArray objects}

\description{
  Wrapping an array-like object (typically an on-disk object) in a
  DelayedArray object allows one to perform common array operations in an
  efficient way. Operations on DelayedArray objects are implemented as delayed
  operations and/or via a block processing strategy in order to optimize
  performance and reduce memory usage.
}

\section{Coercion}{
  To wrap an array-like object \code{a} in a DelayedArray object, coerce it
  to DelayedArray with \code{as(a, "DelayedArray")}, or to DelayedMatrix with
  \code{as(a, "DelayedMatrix")} if it's a matrix-like object.

  To \emph{realize} it (i.e. trigger execution of the delayed operations
  and turn it back into an ordinary array), call \code{as.array} on it.
  Alternatively, write it to disk by coercing it to HDF5Array (other
  on-disk backends can be supported). In that case, the full object is
  not realized at once in memory, but split into blocks first, and blocks
  are realized one at a time and written to disk.
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
    \item \link[base]{array} objects in base R.

    \item \link{HDF5Dataset} and \link{HDF5Array} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## WITH AN ORDINARY array OBJECT
## ---------------------------------------------------------------------

a <- array(runif(15000000), c(300000, 10, 5))
A <- as(a, "DelayedArray")
A

toto <- function(x) (5 * x[ , , 1] ^ 3 + 1L) * log(x[, , 2])
b <- toto(a)
head(b)

B <- toto(A)  # very fast! (operations are delayed)
B             # still 3 dimensions (subsetting never drops dimensions)
as.array(B[1:6, , ])

## The as.array method for DelayedArray (as.array.DelayedArray) has
## a 'drop' argument that is FALSE by default. Set it to TRUE to drop
## the ineffective dimensions:
as.array(B[1:6, , ], drop=TRUE)  # same as drop(as.array(B[1:6, , ]))

M <- as(B, "DelayedMatrix")
M

cs <- colSums(b)
CS <- colSums(M)
stopifnot(identical(cs, CS))

## ---------------------------------------------------------------------
## WITH A HDF5Dataset OBJECT
## ---------------------------------------------------------------------

}
\keyword{methods}
\keyword{classes}
