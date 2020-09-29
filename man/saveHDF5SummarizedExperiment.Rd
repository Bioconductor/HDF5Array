\name{saveHDF5SummarizedExperiment}

\alias{saveHDF5SummarizedExperiment}
\alias{quickResaveHDF5SummarizedExperiment}
\alias{loadHDF5SummarizedExperiment}

\title{Save/load an HDF5-based SummarizedExperiment object}

\description{
  \code{saveHDF5SummarizedExperiment} and
  \code{loadHDF5SummarizedExperiment} can be used to save/load an HDF5-based
  \link[SummarizedExperiment]{SummarizedExperiment} object to/from disk.

  NOTE: These functions use functionalities from the
  \pkg{SummarizedExperiment} package internally and so require this
  package to be installed.
}

\usage{
saveHDF5SummarizedExperiment(x, dir="my_h5_se", prefix="", replace=FALSE,
                                chunkdim=NULL, level=NULL, as.sparse=NA,
                                verbose=FALSE)

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
  \item{as.sparse}{
    Whether the assay data should be flagged as sparse or not. If set to
    \code{NA} (the default), then the specific \code{as.sparse} value to
    use for each assay is determined by calling \code{is_sparse()} on them.

    Passed to the internal calls to \code{writeHDF5Array}.
    See \code{?\link{writeHDF5Array}} for more information and an
    IMPORTANT NOTE.
  }
  \item{verbose}{
    Set to \code{TRUE} to make the function display progress.
  }
}

\details{
  \describe{
    \item{\code{saveHDF5SummarizedExperiment()}:}{
      Creates the directory specified thru the \code{dir} argument and
      populates it with the HDF5 datasets (one per assay in \code{x})
      plus a serialized version of \code{x} that contains pointers to
      these datasets. This directory provides a self-contained HDF5-based
      representation of \code{x} that can then be loaded back in R with
      \code{loadHDF5SummarizedExperiment}.

      Note that this directory is \emph{relocatable} i.e. it can be moved
      (or copied) to a different place, on the same or a different
      computer, before calling \code{loadHDF5SummarizedExperiment} on it.
      For convenient sharing with collaborators, it is suggested to turn
      it into a tarball (with Unix command \code{tar}), or zip file,
      before the transfer.

      Please keep in mind that \code{saveHDF5SummarizedExperiment} and
      \code{loadHDF5SummarizedExperiment} don't know how to produce/read
      tarballs or zip files at the moment, so the process of
      packaging/extracting the tarball or zip file is entirely the user
      responsibility. This is typically done from outside R.

      Finally please note that, depending on the size of the data to
      write to disk and the performance of the disk,
      \code{saveHDF5SummarizedExperiment} can take a long time to complete.
      Use \code{verbose=TRUE} to see its progress.
    }
    \item{\code{loadHDF5SummarizedExperiment()}:}{
       Typically very fast, even if the assay data is big, because all
       the assays in the returned object are \link{HDF5Array} objects
       pointing to the on-disk HDF5 datasets located in \code{dir}.
       \link{HDF5Array} objects are typically light-weight in memory.
    }
    \item{\code{quickResaveHDF5SummarizedExperiment()}:}{
      Preserves the HDF5 file and datasets that the assays in \code{x}
      are already pointing to (and which were created by an earlier call
      to \code{saveHDF5SummarizedExperiment}). All it does is re-serialize
      \code{x} on top of the \code{.rds} file that is associated with
      this HDF5 file (and which was created by an earlier call to
      \code{saveHDF5SummarizedExperiment} or
      \code{quickResaveHDF5SummarizedExperiment}). Because the delayed
      operations possibly carried by the assays in \code{x} are not
      realized, this is very fast.
    }
  }
}

\value{
  \code{saveHDF5SummarizedExperiment} returns an invisible
  \link[SummarizedExperiment]{SummarizedExperiment} object that is the
  same as what \code{loadHDF5SummarizedExperiment} will return when loading
  back the object.
  All the assays in the object are \link{HDF5Array} objects pointing to
  datasets in the HDF5 file saved in \code{dir}.
}

\section{Difference between saveHDF5SummarizedExperiment() and saveRDS()}{
  Roughly speaking, \code{saveRDS()} only serializes the part of an object
  that resides in memory (the reality is a little bit more nuanced, but
  discussing the full details is not important here, and would only distract
  us). For most objects in R, that's the whole object, so \code{saveRDS()}
  does the job.

  However some objects are pointing to on-disk data.
  For example: a \link[GenomicFeatures]{TxDb} object (the
  \link[GenomicFeatures]{TxDb} class is implemented and documented
  in the \pkg{GenomicFeatures} package) points to an SQLite db;
  an \link{HDF5Array} object points to a dataset in an HDF5 file;
  a \link[SummarizedExperiment]{SummarizedExperiment} derivative can
  have one or more of its assays that point to datasets (one per assay)
  in an HDF5 file. These objects have 2 parts: one part is in memory, and
  one part is on disk. The 1st part is sometimes called the \emph{object
  shell} and is generally thin (i.e. it has a small memory footprint).
  The 2nd part is the data and is typically big. The object shell and data
  are linked together via some kind of pointer stored in the shell (e.g.
  an SQLite connection, or a path to a file, etc...). Note that this is a
  \emph{one way link} in the sense that the object shell "knows" where
  to find the on-disk data but the on-disk data knows nothing about the
  object shell (and is completely agnostic about what kind of object
  shell could be pointing to it). Furthermore, at any given time on a
  given system, there could be more than one object shell pointing to
  the same on-disk data. These object shells could exist in the same R
  session or in sessions in other languages (e.g. Python). These various
  sessions could be run by the same or by different users.

  Using \code{saveRDS()} on such object will only serialize the shell part
  so will produce a small \code{.rds} file that contains the serialized
  object shell but not the object data.

  This is problematic because:
  \enumerate{
    \item If you later unserialize the object (with \code{readRDS()})
          on the same system where you originally serialized it, it is
          possible that you will get back an object that is fully functional
          and semantically equivalent to the original object. But here is
          the catch: this will be the case ONLY if the data is still at
          the original location and has not been modified (i.e. nobody
          wrote or altered the data in the SQLite db or HDF5 file in the
          mean time), and if the serialization/unserialization cycle
          didn't break the link between the object shell and the data
          (this serialization/unserialization cycle is known to break open
          SQLite connections).

    \item After serialization the object shell and data are stored in
          separate files (in the new \code{.rds} file for the shell,
          still in the original SQLite or HDF5 file for the data),
          typically in very different places on the file system. But
          these 2 files are not relocatable, that is, moving or copying
          them to another system or sending them to collaborators will
          typically break the link between them. Concretely this means
          that the object obtained by using \code{readRDS()} on the
          destination system will be broken.
  }

  \code{saveHDF5SummarizedExperiment()} addresses these issues by saving
  the object shell and assay data in a folder that is relocatable.

  Note that it only works on \link[SummarizedExperiment]{SummarizedExperiment}
  derivatives. What it does exactly is (1) write all the assay data to an
  HDF5 file, and (2) serialize the object shell, which in this case is
  everything in the object that is not the assay data. The 2 files
  (HDF5 and \code{.rds}) are written to the directory specified by the user.
  The resulting directory contains a full representation of the object and
  is relocatable, that is, it can be moved or copied to another place on
  the system, or to another system (possibly after making a tarball of it),
  where \code{loadHDF5SummarizedExperiment()} can then be used to load the
  object back in R.
}

\note{
  The files created by \code{saveHDF5SummarizedExperiment} in the
  user-specified directory \code{dir} should not be renamed.

  The user-specified \emph{directory} created by
  \code{saveHDF5SummarizedExperiment} is relocatable i.e. it can
  be renamed and/or moved around, but not the individual files in it.
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

    \item \code{base::\link[base]{saveRDS}}
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
assays(se2, withDimnames=FALSE) <- c(assays(se2), list(score=assay1/100))
rowRanges(se2) <- GRanges("chr1", IRanges(1:14, width=5))
rownames(se2) <- letters[1:14]
se2

## This will replace saved 'se1'!
quickResaveHDF5SummarizedExperiment(se2, verbose=TRUE)
list.files(dir)
loadHDF5SummarizedExperiment(dir, prefix="xx_")
}
