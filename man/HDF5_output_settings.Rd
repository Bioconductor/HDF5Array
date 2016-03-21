\name{HDF5_output_settings}

\alias{HDF5_output_settings}

\alias{setHDF5OutputFile}
\alias{getHDF5OutputFile}
\alias{lsHDF5OutputFile}

\alias{setHDF5OutputName}
\alias{getHDF5OutputName}

\title{Manage settings for writing output to HDF5}

\description{
  Utility functions to manage settings for writing output to HDF5.
}

\usage{
setHDF5OutputFile(file=paste0(tempfile(), ".h5"))
getHDF5OutputFile()
lsHDF5OutputFile()

setHDF5OutputName(name)
getHDF5OutputName()
}

\arguments{
  \item{file}{
    The path to the \emph{current output HDF5 file}, that is, the HDF5 file
    where to write all newly created datasets.
  }
  \item{name}{
    The name of the \emph{next} dataset to be written to the current output
    HDF5 file. This is for a one-time use only.
  }
}

\seealso{
  \itemize{
    \item \link{DelayedArray} objects.

    \item \link{HDF5Array} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
getHDF5OutputFile()
lsHDF5OutputFile()
}
\keyword{methods}
