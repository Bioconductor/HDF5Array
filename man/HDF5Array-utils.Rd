\name{HDF5Array-utils}

\alias{HDF5Array-utils}
\alias{.shorten_assay2h5_links}
\alias{.restore_absolute_assay2h5_links}
\alias{validate_HDF5ArraySeed_dataset}
\alias{.write_h5_assays}
\alias{.create_dir}
\alias{.replace_dir}
\alias{.check_and_delete_files}
\alias{.stop_if_bad_dir}

\title{Utilities for HDF5Array-dependent packages}

\description{
    \code{.shorten_assay2h5_links} and \code{.restore_absolute_assay2h5_links}
    are used to serialize the HDF5-based assays in the object and to load
    \code{\link{HDF5ArraySeed}} object and restore the file paths to their
    absolute form.  The \code{HDF5ArraySeed} object is also validated during
    link-restore with \code{validate_HDF5ArraySeed_dataset}.

    \code{.write_h5_assays} is a low-level function for taking the assays
    part of the object along with a location and additional arguments to
    save the assays into the H5 file with \code{\link{writeHDF5Array}}.

    Finally, \code{.create_dir}, \code{.replace_dir},
    \code{.check_and_delete_files}, and \code{.stop_if_bad_dir} are all
    file and directory management utitilty functions for ensuring
    consistency in writing of the file.
}

\keyword{internal}
