\name{h5writeDimnames}

\alias{h5writeDimnames}
\alias{h5readDimnames}
\alias{set_h5dimnames}
\alias{get_h5dimnames}

\title{Write/read the dimnames of an HDF5 dataset}

\description{
  \code{h5writeDimnames} and \code{h5readDimnames} can be used to
  write/read the dimnames of an HDF5 dataset to/from the HDF5 file.

  Note that \code{h5writeDimnames} is used internally by
  \code{\link{writeHDF5Array}(x, ..., with.dimnames=TRUE)}
  to write the dimnames of \code{x} to the HDF5 file together
  with the array data.

  \code{set_h5dimnames} and \code{get_h5dimnames} are low-level
  utilities that can be used to attach existing HDF5 datasets
  along the dimensions of a given HDF5 dataset, or to retrieve
  the names of the HDF5 datasets that are attached along the
  dimensions of a given HDF5 dataset.
}

\usage{
h5writeDimnames(dimnames, filepath, name, group=NA, h5dimnames=NULL)
h5readDimnames(filepath, name, as.character=FALSE)

set_h5dimnames(filepath, name, h5dimnames, dry.run=FALSE)
get_h5dimnames(filepath, name)
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
    For \code{h5writeDimnames} and \code{h5readDimnames}: The path (as a
    single string) to the HDF5 file where the dimnames should be written
    to or read from.

    For \code{set_h5dimnames} and \code{get_h5dimnames}: The path (as a
    single string) to the HDF5 file where to set or get the
    \emph{h5dimnames}.
  }
  \item{name}{
    For \code{h5writeDimnames} and \code{h5readDimnames}: The name of the
    dataset in the HDF5 file for which the dimnames should be written or
    read.

    For \code{set_h5dimnames} and \code{get_h5dimnames}: The name of the
    dataset in the HDF5 file for which to set or get the \emph{h5dimnames}.
  }
  \item{group}{
    \code{NA} (the default) or the name of the HDF5 group where to write
    the dimnames.
    If set to \code{NA} then the group name is automatically generated
    from \code{name}.
    If set to the empty string (\code{""}) then no group will be used.

    Except when \code{group} is set to the empty string, the names in
    \code{h5dimnames} (see below) must be relative to the group.
  }
  \item{h5dimnames}{
    For \code{h5writeDimnames}: \code{NULL} (the default) or a character
    vector containing the names of the HDF5 datasets (one per list
    element in \code{dimnames}) where to write the dimnames.
    Names associated with \code{NULL} list elements in \code{dimnames}
    are ignored and should typically be NAs.

    If set to \code{NULL} then the names are automatically set to numbers
    indicating the associated dimensions (\code{"1"} for the first dimension,
    \code{"2"} for the second, etc...)

    For \code{set_h5dimnames}: A character vector containing the names
    of the HDF5 datasets to attach as dimnames of the dataset specified
    in \code{name}. The vector must have one element per dimension in
    dataset \code{name}. NAs are allowed and indicate dimensions along
    which nothing should be attached.
  }
  \item{as.character}{
    Even though the dimnames of an HDF5 dataset are usually stored as
    datasets of type \code{"character"} (H5 datatype \code{"H5T_STRING"})
    in the HDF5 file, this is not a requirement.
    By default \code{h5readDimnames} will return them \emph{as-is}.
    Set \code{as.character} to \code{TRUE} to make sure that they are
    returned as character vectors. See example below.
  }
  \item{dry.run}{
    When set to \code{TRUE}, \code{set_h5dimnames} doesn't make any
    change to the HDF5 file but will still raise errors if the operation
    cannot be done.
  }
}

\value{
  \code{h5writeDimnames} and \code{set_h5dimnames} return nothing.

  \code{h5readDimnames} returns a list (possibly named) with one list
  element per dimension in HDF5 dataset \code{name} and containing its
  dimnames retrieved from the file.

  \code{get_h5dimnames} returns a character vector containing the names
  of the HDF5 datasets that are currently set as the dimnames of the
  dataset specified in \code{name}. The vector has one element per
  dimension in dataset \code{name}. NAs in the vector indicate dimensions
  along which nothing is set.
}

\seealso{
  \itemize{
    \item \code{\link{writeHDF5Array}} for a high-level function to write
          an array-like object and its dimnames to an HDF5 file.

    \item \code{\link[rhdf5]{h5write}} in the \pkg{rhdf5} package that
          \code{h5writeDimnames} uses internally to write the dimnames
          to the HDF5 file.

    \item \code{\link{h5mread}} in this package (\pkg{HDF5Array}) that
          \code{h5readDimnames} uses internally to read the dimnames
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
h5write(m0, h5file, "M0")  # h5write() ignores the dimnames
h5ls(h5file)

h5writeDimnames(dimnames(m0), h5file, "M0")
h5ls(h5file)

get_h5dimnames(h5file, "M0")
h5readDimnames(h5file, "M0")

## Reconstruct 'm0' from HDF5 file:
m1 <- h5mread(h5file, "M0")
dimnames(m1) <- h5readDimnames(h5file, "M0")
stopifnot(identical(m0, m1))

## Create an HDF5Array object that points to HDF5 dataset M0:
HDF5Array(h5file, "M0")

## Sanity checks:
stopifnot(identical(dimnames(m0), h5readDimnames(h5file, "M0")))
stopifnot(identical(dimnames(m0), dimnames(HDF5Array(h5file, "M0"))))

## ---------------------------------------------------------------------
## SHARED DIMNAMES
## ---------------------------------------------------------------------
## If a collection of HDF5 datasets share the same dimnames, the
## dimnames only need to be written once in the HDF5 file. Then they
## can be attached to the individual datasets with set_h5dimnames():

h5write(array(runif(240), c(12, 5:4)), h5file, "A1")
set_h5dimnames(h5file, "A1", get_h5dimnames(h5file, "M0"))
get_h5dimnames(h5file, "A1")
h5readDimnames(h5file, "A1")
HDF5Array(h5file, "A1")

h5write(matrix(sample(letters, 60, replace=TRUE), ncol=5), h5file, "A2")
set_h5dimnames(h5file, "A2", get_h5dimnames(h5file, "M0"))
get_h5dimnames(h5file, "A2")
h5readDimnames(h5file, "A2")
HDF5Array(h5file, "A2")

## Sanity checks:
stopifnot(identical(dimnames(m0), h5readDimnames(h5file, "A1")[1:2]))
stopifnot(identical(dimnames(m0), h5readDimnames(h5file, "A2")))

## ---------------------------------------------------------------------
## USE h5writeDimnames() AFTER A CALL TO writeHDF5Array()
## ---------------------------------------------------------------------
## After calling writeHDF5Array(x, ..., with.dimnames=FALSE) the
## dimnames on 'x' can still be written to the HDF5 file by doing the
## following:

## 1. Write 'm0' to the HDF5 file and ignore the dimnames (for now):
writeHDF5Array(m0, h5file, "M2")

## 2. Use h5writeDimnames() to write 'dimnames(m0)' to the file and
##    associate them with the "M2" dataset:
h5writeDimnames(dimnames(m0), h5file, "M2")

## 3. Use the HDF5Array() constructor to make an HDF5Array object that
##    points to the "M2" dataset:
HDF5Array(h5file, "M2")

## Note that at step 2. you can use the extra arguments of
## h5writeDimnames() to take full control of where the dimnames
## should be stored in the file:
writeHDF5Array(m0, h5file, "M3")
h5writeDimnames(dimnames(m0), h5file, "M3",
                 group="a_secret_place", h5dimnames=c("NA", "M3_dim2"))
h5ls(h5file)
## h5readDimnames() and HDF5Array() still "finds" the dimnames:
h5readDimnames(h5file, "M3")
HDF5Array(h5file, "M3")

## Sanity checks:
stopifnot(identical(dimnames(m0), h5readDimnames(h5file, "M3")))
stopifnot(identical(dimnames(m0), dimnames(HDF5Array(h5file, "M3"))))

## ---------------------------------------------------------------------
## STORE THE DIMNAMES AS NON-CHARACTER TYPES
## ---------------------------------------------------------------------
writeHDF5Array(m0, h5file, "M4")
dimnames <- list(1001:1012, as.raw(11:15))
h5writeDimnames(dimnames, h5file, "M4")
h5ls(h5file)

h5readDimnames(h5file, "M4")
h5readDimnames(h5file, "M4", as.character=TRUE)

## Sanity checks:
stopifnot(identical(dimnames, h5readDimnames(h5file, "M4")))
dimnames(m0) <- dimnames
stopifnot(identical(
    dimnames(m0),
    h5readDimnames(h5file, "M4", as.character=TRUE)
))
}
\keyword{methods}
