\name{saveHDF5SummarizedExperiment}

\alias{saveHDF5SummarizedExperiment}
\alias{quickResaveHDF5SummarizedExperiment}
\alias{loadHDF5SummarizedExperiment}

\title{Save/load an HDF5-based SummarizedExperiment object}

\description{
  \code{saveHDF5SummarizedExperiment} and
  \code{loadHDF5SummarizedExperiment} can be used to save/load an HDF5-based
  \link[SummarizedExperiment]{SummarizedExperiment} object to/from disk.
}

\usage{
saveHDF5SummarizedExperiment(x, dir="my_h5_se", prefix="", replace=FALSE,
                                chunkdim=NULL, level=NULL, verbose=FALSE)

loadHDF5SummarizedExperiment(dir="my_h5_se", prefix="")

quickResaveHDF5SummarizedExperiment(x, verbose=FALSE)
}

\arguments{
  \item{x}{
    A \link[SummarizedExperiment]{SummarizedExperiment} object or derivative.

    For \code{quickResaveHDF5SummarizedExperiment} the object must have been
    previously saved with \code{saveHDF5SummarizedExperiment} (and has been
    possibly modified since then).
  }
  \item{dir}{
    The path (as a single string) to the directory where to save the
    HDF5-based \link[SummarizedExperiment]{SummarizedExperiment} object
    or to load it from.

    When saving, the directory will be created if it doesn't already exist.
    If the directory already exists and no prefix is specified and
    \code{replace} is set to \code{TRUE}, then it's replaced with an
    empty directory.
  }
  \item{prefix}{
    An optional prefix to add to the names of the files created
    inside \code{dir}. Allows saving more than one object in the same
    directory.
  }
  \item{replace}{
    When no prefix is specified, should a pre-existing directory be replaced
    with a new empty one? The content of the pre-existing directory will
    be lost!
  }
  \item{chunkdim, level}{
    The dimensions of the chunks and the compression level to use for
    writing the assay data to disk.
    Passed to the internal calls to \code{writeHDF5Array}.
    See \code{?\link{writeHDF5Array}} for more information.
  }
  \item{verbose}{
    Set to \code{TRUE} to make the function display progress.
  }
}

\details{
  These functions use functionalities from the \pkg{SummarizedExperiment}
  package internally and so require this package to be installed.

  \code{saveHDF5SummarizedExperiment} creates the directory specified
  thru the \code{dir} argument and then populates it with the HDF5 datasets
  (one per assay in \code{x}) plus a serialized version of \code{x} that
  contains pointers to these datasets. This directory provides a
  self-contained HDF5-based representation of \code{x} that can then
  be loaded back in R with \code{loadHDF5SummarizedExperiment}.
  Note that this directory is \emph{relocatable} i.e. it can be moved
  (or copied) to a different place, on the same or a different computer,
  before calling \code{loadHDF5SummarizedExperiment} on it. For convenient
  sharing with collaborators, it is suggested to turn it into a tarball
  (with Unix command \code{tar}), or zip file, before the transfer.
  Please keep in mind that \code{saveHDF5SummarizedExperiment} and
  \code{loadHDF5SummarizedExperiment} don't know how to produce/read
  tarballs or zip files at the moment, so the process of packaging/extracting
  the tarball or zip file is entirely the user responsibility. It is
  typically done from outside R.

  Finally please note that, depending on the size of the data to write to
  disk and the performance of the disk, \code{saveHDF5SummarizedExperiment}
  can take a long time to complete. Use \code{verbose=TRUE} to see its
  progress.

  \code{loadHDF5SummarizedExperiment} is generally very fast, even if
  the assay data is big, because all the assays in the returned object
  are \link{HDF5Array} objects pointing to the on-disk HDF5 datasets
  located in \code{dir}. \link{HDF5Array} objects are typically
  light-weight in memory.

  \code{quickResaveHDF5SummarizedExperiment} preserves the HDF5 file and
  datasets that the assays in \code{x} are already pointing to (and which
  were created by an earlier call to \code{saveHDF5SummarizedExperiment}).
  All it does is re-serialize \code{x} on top of the \code{.rds} file
  that is associated with this HDF5 file (and which was created by an
  earlier call to \code{saveHDF5SummarizedExperiment} or
  \code{quickResaveHDF5SummarizedExperiment}). Because the delayed
  operations possibly carried by the assays in \code{x} are not realized,
  this is very fast.
}

\value{
  \code{saveHDF5SummarizedExperiment} returns an invisible
  \link[SummarizedExperiment]{SummarizedExperiment} object where all
  the assays are \link{HDF5Array} objects pointing to the HDF5 datasets
  saved in \code{dir}. It's in fact the same obect as the object that
  would be returned by calling \code{loadHDF5SummarizedExperiment} on
  \code{dir}.
}

\note{
  The files created by \code{saveHDF5SummarizedExperiment} should not
  be renamed.
  The \emph{directory} created by \code{saveHDF5SummarizedExperiment} is
  relocatable i.e. it can be renamed and/or moved around, but not the files
  in it.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link[SummarizedExperiment]{SummarizedExperiment} and
          \link[SummarizedExperiment]{RangedSummarizedExperiment}
          objects in the \pkg{SummarizedExperiment} package.

    \item The \code{\link{writeHDF5Array}} function which
          \code{saveHDF5SummarizedExperiment} uses internally to write
          the assay data to disk.
  }
}

\examples{
## ---------------------------------------------------------------------
## saveHDF5SummarizedExperiment() / loadHDF5SummarizedExperiment()
## ---------------------------------------------------------------------
library(SummarizedExperiment)

nrow <- 200
ncol <- 6
counts <- matrix(as.integer(runif(nrow * ncol, 1, 1e4)), nrow)
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
se0 <- SummarizedExperiment(assays=list(counts=counts), colData=colData)
se0

## Save 'se0' as an HDF5-based SummarizedExperiment object:
dir <- tempfile("h5_se0_")
h5_se0 <- saveHDF5SummarizedExperiment(se0, dir)
list.files(dir)

h5_se0
assay(h5_se0, withDimnames=FALSE)   # HDF5Matrix object

h5_se0b <- loadHDF5SummarizedExperiment(dir)
h5_se0b
assay(h5_se0b, withDimnames=FALSE)  # HDF5Matrix object

## Sanity checks:
stopifnot(is(assay(h5_se0, withDimnames=FALSE), "HDF5Matrix"))
stopifnot(identical(assay(se0), as.matrix(assay(h5_se0))))
stopifnot(is(assay(h5_se0b, withDimnames=FALSE), "HDF5Matrix"))
stopifnot(identical(assay(se0), as.matrix(assay(h5_se0b))))

## ---------------------------------------------------------------------
## More sanity checks
## ---------------------------------------------------------------------

## Make a copy of directory 'dir':
somedir <- tempfile("somedir")
dir.create(somedir)
file.copy(dir, somedir, recursive=TRUE)
dir2 <- list.files(somedir, full.names=TRUE)

## 'dir2' contains a copy of 'dir'. Call loadHDF5SummarizedExperiment()
## on it.
h5_se0c <- loadHDF5SummarizedExperiment(dir2)

stopifnot(is(assay(h5_se0c, withDimnames=FALSE), "HDF5Matrix"))
stopifnot(identical(assay(se0), as.matrix(assay(h5_se0c))))

## ---------------------------------------------------------------------
## Using a prefix
## ---------------------------------------------------------------------

se1 <- se0[51:100, ]
saveHDF5SummarizedExperiment(se1, dir, prefix="xx_")
list.files(dir)
loadHDF5SummarizedExperiment(dir, prefix="xx_")

## ---------------------------------------------------------------------
## quickResaveHDF5SummarizedExperiment()
## ---------------------------------------------------------------------

se2 <- loadHDF5SummarizedExperiment(dir, prefix="xx_")
se2 <- se2[1:14, ]
assay1 <- assay(se2, withDimnames=FALSE)
assays(se2) <- c(assays(se2), list(score=assay1/100))
rowRanges(se2) <- GRanges("chr1", IRanges(1:14, width=5))
rownames(se2) <- letters[1:14]
se2

## This will replace saved 'se1'!
quickResaveHDF5SummarizedExperiment(se2, verbose=TRUE)
list.files(dir)
loadHDF5SummarizedExperiment(dir, prefix="xx_")
}
