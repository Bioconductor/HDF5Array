\name{write_h5dimnames}

\alias{write_h5dimnames}
\alias{read_h5dimnames}

\title{Write/read the dimnames of an HDF5 dataset}

\description{
  \code{write_h5dimnames} and \code{read_h5dimnames} can be used to
  write/read the dimnames of an HDF5 dataset to/from disk.
}

\usage{
write_h5dimnames(dimnames, filepath, name, group=NA, dimscales=NULL)

read_h5dimnames(filepath, name)
}

\arguments{
  \item{dimnames}{
    The dimnames to write to the HDF5 file. Must be supplied as a list
    (possibly named) with one list element per dimension in the HDF5 dataset
    specified via the \code{name} argument.
    Each list element in \code{dimnames} must be an atomic vector or a
    \code{NULL}. When not a \code{NULL}, its length must equal the extent
    of the corresponding dimension in the HDF5 dataset.
  }
  \item{filepath}{
    The path (as a single string) to the HDF5 file where the dimnames
    should be written to or read from.
  }
  \item{name}{
    The name of the dataset in the HDF5 file for which the dimnames should
    be written or read.
  }
  \item{group}{
    \code{NA} (the default) or the name of the HDF5 group where to write
    the dimnames.
    If set to \code{NA} then the group name is automatically generated
    from \code{name}.
    If set to the empty string (\code{""}) then no group will be used.

    Except when \code{group} is set to the empty string, the names in
    \code{dimscales} (see below) must be relative to the group.
  }
  \item{dimscales}{
    \code{NULL} (the default) or a character vector containing the names
    of the HDF5 datasets (one per list element in \code{dimnames}) where
    to write the dimnames. Names associated with \code{NULL} list elements
    in \code{dimnames} are ignored and should typically be NAs.

    If set to \code{NULL} then the names are automatically set to numbers
    indicating the associated dimensions (\code{"1"} for the first dimension,
    \code{"2"} for the second, etc...)

    Note that in HDF5 jargon datasets that represent attributes along the
    dimensions of other datasets are called \emph{Dimension Scales dataset}
    hence the name of this argument.
  }
}

\value{
  \code{write_h5dimnames} returns nothing.

  \code{read_h5dimnames} returns a list (possibly named) with one list
  element per dimension in HDF5 dataset \code{name} and containing its
  dimnames retrieved from the file.
}

\seealso{
  \itemize{
    \item \code{\link[rhdf5]{h5write}} in the \pkg{rhdf5} package that
          \code{write_h5dimnames} uses internally to write the dimnames
          to the HDF5 file.

    \item \code{\link{h5mread}} in this package (\pkg{HDF5Array}) that
          \code{read_h5dimnames} uses internally to read the dimnames
          from the HDF5 file.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.
  }
}

\examples{
library(rhdf5)  # for h5write() and h5ls()

m0 <- matrix(1:60, ncol=5)
colnames(m0) <- LETTERS[1:5]
h5file <- tempfile(fileext=".h5")
h5write(m0, h5file, "m0")  # h5write() ignores the dimnames
h5ls(h5file)

write_h5dimnames(dimnames(m0), h5file, "m0")
h5ls(h5file)

read_h5dimnames(h5file, "m0")

## Reconstruct 'm0' from HDF5 file.
m <- h5mread(h5file, "m0")
dimnames(m) <- read_h5dimnames(h5file, "m0")
stopifnot(identical(m0, m))
}
\keyword{methods}
