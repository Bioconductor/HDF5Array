\name{write_h5dimnames}

\alias{write_h5dimnames}
\alias{read_h5dimnames}

\title{Write/read the dimnames of an HDF5 dataset}

\description{
  Low-level functions \code{write_h5dimnames} and \code{read_h5dimnames}
  can be used to write/read the dimnames of an HDF5 dataset to/from
  the HDF5 file.

  Note that \code{write_h5dimnames} is used internally by
  \code{\link{writeHDF5Array}(x, ..., with.dimnames=TRUE)}
  to write the dimnames of \code{x} together with the array
  data to the HDF5 file.
}

\usage{
write_h5dimnames(dimnames, filepath, name, group=NA, dimscales=NULL)

read_h5dimnames(filepath, name, as.character=FALSE)
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
    dimensions of other datasets are called \emph{Dimension Scales dataset},
    hence the name of this argument.
  }
  \item{as.character}{
    Even though the dimnames of an HDF5 dataset are usually stored as
    datasets of type \code{"character"} (H5 datatype \code{"H5T_STRING"})
    in the HDF5 file, this is not a requirement.
    By default \code{read_h5dimnames} will return them \emph{as-is}.
    Set \code{as.character} to \code{TRUE} to make sure that they are
    returned as character vectors. See example below.
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
    \item \code{\link{writeHDF5Array}} for a high-level function to write
          an array-like object and its dimnames to an HDF5 file.

    \item \code{\link[rhdf5]{h5write}} in the \pkg{rhdf5} package that
          \code{write_h5dimnames} uses internally to write the dimnames
          to the HDF5 file.

    \item \code{\link{h5mread}} in this package (\pkg{HDF5Array}) that
          \code{read_h5dimnames} uses internally to read the dimnames
          from the HDF5 file.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item \link{HDF5Array} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------
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
m1 <- h5mread(h5file, "m0")
dimnames(m1) <- read_h5dimnames(h5file, "m0")
stopifnot(identical(m0, m1))

## ---------------------------------------------------------------------
## USE write_h5dimnames() AFTER A CALL TO writeHDF5Array()
## ---------------------------------------------------------------------

## After calling writeHDF5Array(x, ..., with.dimnames=FALSE) the
## dimnames on 'x' can still be written to the HDF5 file by doing the
## following:

## 1. Write 'm0' to the HDF5 file and ignore the dimnames (for now):
writeHDF5Array(m0, h5file, "m2")

## 2. Use write_h5dimnames() to write 'dimnames(m0)' to the file and
##    associate them with the "m2" dataset:
write_h5dimnames(dimnames(m0), h5file, "m2")

## 3. Use the HDF5Array() constructor to make an HDF5Array object that
##    points to the "m2" dataset:
HDF5Array(h5file, "m2")

## Note that at step 2. you can use the extra arguments of
## write_h5dimnames() to take full control of where the dimnames should
## be stored in the file:
writeHDF5Array(m0, h5file, "m3")
write_h5dimnames(dimnames(m0), h5file, "m3",
                 group="a_secret_place", dimscales=c("NA", "m3_dim2"))
h5ls(h5file)
## read_h5dimnames() and HDF5Array() still "finds" the dimnames:
read_h5dimnames(h5file, "m3")
HDF5Array(h5file, "m3")

## Sanity checks:
stopifnot(identical(dimnames(m0), read_h5dimnames(h5file, "m3")))
stopifnot(identical(dimnames(m0), dimnames(HDF5Array(h5file, "m3"))))

## ---------------------------------------------------------------------
## USE A NON-CHARACTER TYPE TO STORE THE DIMNAMES
## ---------------------------------------------------------------------

writeHDF5Array(m0, h5file, "m4")
dimnames <- list(1001:1012, as.raw(11:15))
write_h5dimnames(dimnames, h5file, "m4")
h5ls(h5file)

read_h5dimnames(h5file, "m4")
read_h5dimnames(h5file, "m4", as.character=TRUE)

## Sanity checks:
stopifnot(identical(dimnames, read_h5dimnames(h5file, "m4")))
dimnames(m0) <- dimnames
stopifnot(identical(
    dimnames(m0),
    read_h5dimnames(h5file, "m4", as.character=TRUE)
))
}
\keyword{methods}
