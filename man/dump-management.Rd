\name{HDF5-dump-management}

\alias{HDF5-dump-management}
\alias{dump-management}

\alias{setHDF5DumpDir}
\alias{getHDF5DumpDir}

\alias{setHDF5DumpFile}
\alias{getHDF5DumpFile}
\alias{lsHDF5DumpFile}

\alias{setHDF5DumpName}
\alias{getHDF5DumpName}

\alias{getHDF5DumpChunkDim}

\alias{setHDF5DumpCompressionLevel}
\alias{getHDF5DumpCompressionLevel}

\alias{appendDatasetCreationToHDF5DumpLog}
\alias{showHDF5DumpLog}

\title{HDF5 dump management}

\description{
  A set of utilities to control the location of automatically created HDF5
  datasets.
}

\usage{
setHDF5DumpDir(dir)
setHDF5DumpFile(filepath)
setHDF5DumpName(name)
setHDF5DumpCompressionLevel(level=6L)

getHDF5DumpDir()
getHDF5DumpFile(for.use=FALSE)
getHDF5DumpName(for.use=FALSE)
getHDF5DumpCompressionLevel()

lsHDF5DumpFile()

showHDF5DumpLog()

## For developers:
getHDF5DumpChunkDim(dim, type, ratio=75)
appendDatasetCreationToHDF5DumpLog(filepath, name, dim, type,
                                   chunkdim, level)
}

\arguments{
  \item{dir}{
    The path (as a single string) to the current \emph{HDF5 dump directory},
    that is, to the (new or existing) directory where \emph{HDF5 dump files}
    with automatic names will be created. This is ignored if the user
    specified an \emph{HDF5 dump file} with \code{setHDF5DumpFile}.
    If \code{dir} is missing, then the \emph{HDF5 dump directory} is set
    back to its default value i.e. to some directory under \code{tempdir()}
    (call \code{getHDF5DumpDir()} to get the exact path).
  }
  \item{filepath}{
    For \code{setHDF5DumpFile}:
    The path (as a single string) to the current \emph{HDF5 dump file},
    that is, to the (new or existing) HDF5 file where the \emph{next
    automatic HDF5 datasets} will be written. If \code{filepath} is
    missing, then a new file with an automatic name will be created
    (in \code{getHDF5DumpDir()}) and used for each new dataset.

    For \code{appendDatasetCreationToHDF5DumpLog}:
    See the Note TO DEVELOPERS below.
  }
  \item{name}{
    For \code{setHDF5DumpName}:
    The name of the \emph{next automatic HDF5 dataset} to be written to
    the current \emph{HDF5 dump file}.

    For \code{appendDatasetCreationToHDF5DumpLog}:
    See the Note TO DEVELOPERS below.
  }
  \item{level}{
    For \code{setHDF5DumpCompressionLevel}:
    The compression level to use for writting \emph{automatic HDF5 datasets}
    to disk.
    See the \code{level} argument in \code{?rhdf5::\link{h5createDataset}}
    (in the \pkg{rhdf5} package) for more information about this.

    For \code{appendDatasetCreationToHDF5DumpLog}:
    See the Note TO DEVELOPERS below.
  }
  \item{for.use}{
    Whether the returned file or dataset name is for use by the caller or not.
    See below for the details.
  }
  \item{dim}{
    The dimensions of the HDF5 dataset to be written to disk, that is,
    an integer vector of length one or more giving the maximal indices in each
    dimension.
    See the \code{dims} argument in \code{?rhdf5::\link{h5createDataset}}
    (in the \pkg{rhdf5} package) for more information about this.
  }
  \item{type}{
    The type (a.k.a. storage mode) of the data to be written to disk. Can be
    obtained with \code{\link[DelayedArray]{type}()} on an array-like object
    (which is equivalent to \code{storage.mode()} or \code{typeof()} on an
    ordinary array). This is typically what an application writing datasets
    to the \emph{HDF5 dump} should pass to the \code{storage.mode} argument
    of its call to \code{rhdf5::\link{h5createDataset}}.
    See the Note TO DEVELOPERS below for more information.
  }
  \item{ratio}{
    The number of chunks per block. By default, a ratio of 75 is used.
  }
  \item{chunkdim}{
    The dimensions of the chunks.
  }
}

\details{
  Calling \code{getHDF5DumpFile()} and \code{getHDF5DumpName()} with no
  argument should be \emph{informative} only i.e. it's a mean for the user
  to know where the \emph{next automatic HDF5 dataset} will be written.
  Since a given file/name combination can be used only once, the user should
  be careful to not use that combination to explicitely create an HDF5 dataset
  because that would get in the way of the creation of the \emph{next
  automatic HDF5 dataset}.
  See the Note TO DEVELOPERS below if you actually need to use this file/name
  combination.

  \code{lsHDF5DumpFile()} is a just convenience wrapper for
  \code{rhdf5::\link[rhdf5]{h5ls}(getHDF5DumpFile())}.
}

\value{
  \code{getHDF5DumpDir} returns the absolute path to the directory
  where \emph{HDF5 dump files} with automatic names will be created.
  Only meaningful if the user did NOT specify an \emph{HDF5 dump file}
  with \code{setHDF5DumpFile}.

  \code{getHDF5DumpFile} returns the absolute path to the HDF5 file where
  the \emph{next automatic HDF5 dataset} will be written.

  \code{getHDF5DumpName} returns the name of the \emph{next automatic HDF5
  dataset}.

  \code{getHDF5DumpCompressionLevel} returns the compression level currently
  used for writting \emph{automatic HDF5 datasets} to disk.

  \code{showHDF5DumpLog} returns the dump log in an invisible data frame.

  \code{getHDF5DumpChunkDim} returns the dimension of the chunks for the
  specified \code{ratio}. By default a ratio of 75 is used i.e. the
  \emph{automatic HDF5 dataset} are written to disk using 75 chunks per block.
}

\note{
  TO DEVELOPERS:

  If your application needs to write its own dataset to the \emph{HDF5 dump}
  then it should:
  \enumerate{
    \item Get a file/name combination by calling
          \code{getHDF5DumpFile(for.use=TRUE)} and
          \code{getHDF5DumpName(for.use=TRUE)}.

    \item [OPTIONAL] Call \code{getHDF5DumpChunkDim(dim, type)} (possibly
          with a non-default \code{ratio}) to get reasonable chunk dimensions
          to use for writting the data to disk. Or choose your own chunk
          dimensions.

    \item Add an entry to the dump log by calling
          \code{appendDatasetCreationToHDF5DumpLog}.
          Typically, this should be done right after creating the dataset
          (e.g. with \code{rhdf5::h5createDataset}) and before starting to
          write the data to disk. The values passed to
          \code{appendDatasetCreationToHDF5DumpLog} via the \code{filepath},
          \code{name}, \code{dim}, \code{type}, \code{chunkdim}, and
          \code{level} arguments should be those that were passed to
          \code{rhdf5::h5createDataset} via the \code{file}, \code{dataset},
          \code{dims}, \code{storage.mode}, \code{chunk}, and \code{level}
          arguments, respectively.
          Note that \code{appendDatasetCreationToHDF5DumpLog} uses a lock
          mechanism so is safe to use in the context of parallel execution.
  }
  This is actually what the coercion method to \link{HDF5Array} does
  internally.
}

\seealso{
  \itemize{
    \item \code{\link{writeHDF5Array}} for writting an array-like object
          to an HDF5 file.

    \item \link{HDF5Array} objects.

    \item The \code{\link[rhdf5]{h5ls}} function in the \pkg{rhdf5} package,
          on which \code{lsHDF5DumpFile} is based.

    \item \code{\link[DelayedArray]{type}} in the \pkg{DelayedArray} package.
  }
}

\examples{
getHDF5DumpDir()
getHDF5DumpFile()

## Use setHDF5DumpFile() to change the current HDF5 dump file.
## If the specified file exists, then it must be in HDF5 format or
## an error will be raised. If it doesn't exist, then it will be
## created.
#setHDF5DumpFile("path/to/some/HDF5/file")

lsHDF5DumpFile()

a <- array(1:600, c(150, 4))
A <- as(a, "HDF5Array")
lsHDF5DumpFile()
A

b <- array(runif(6000), c(4, 2, 150))
B <- as(b, "HDF5Array")
lsHDF5DumpFile()
B

C <- (log(2 * A + 0.88) - 5)^3 * t(drop(B[ , 1, ]))
as(C, "HDF5Array")  # realize C on disk
lsHDF5DumpFile()

## Matrix multiplication is not delayed: the output matrix is realized
## block by block. The current "realization backend" controls where
## realization happens e.g. in memory if set to NULL or in an HDF5 file
## if set to "HDF5Array". See '?realize' in the DelayedArray package for
## more information about "realization backends".
setRealizationBackend("HDF5Array")
m <- matrix(runif(20), nrow=4)
P <- C \%*\% m
lsHDF5DumpFile()

## See all the HDF5 datasets created in the current session so far:
showHDF5DumpLog()

## Wrap the call in suppressMessages() if you are only interested in the
## data frame version of the dump log:
dump_log <- suppressMessages(showHDF5DumpLog())
dump_log
}
\keyword{methods}
