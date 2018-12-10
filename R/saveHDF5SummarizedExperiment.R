### =========================================================================
### Save/load an HDF5-based SummarizedExperiment object
### -------------------------------------------------------------------------


.HOWTO_INSTALL_PKG <- c(
    "You need the SummarizedExperiment package for this. ",
    "Please install\n  it with:\n",
    "      if (!requireNamespace(\"BiocManager\", quietly=TRUE))\n",
    "          install.packages(\"BiocManager\")\n",
    "      BiocManager::install(\"SummarizedExperiment\")"
)

.load_SummarizedExperiment_package <- function()
{
    if (!requireNamespace("SummarizedExperiment", quietly=TRUE))
        stop(.HOWTO_INSTALL_PKG)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .write_HDF5SummarizedExperiment() / .read_HDF5SummarizedExperiment
###

.write_h5_assays <- function(assays, h5_path, chunkdim, level, verbose)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        h5_name <- sprintf("assay%03d", i)
        if (verbose)
            message("Start writing assay ", i, "/", nassay, " to '",
                    h5_path, "':")
        a <- writeHDF5Array(a, h5_path, h5_name, chunkdim, level,
                            verbose=verbose)
        if (verbose)
            message("Finished writing assay ", i, "/", nassay, " to '",
                    h5_path, "'.")
        assays[[i]] <- a
    }
    assays
}

.shorten_h5_paths <- function(assays)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        a@seed@filepath <- basename(a@seed@filepath)
        assays[[i]] <- a
    }
    assays
}

.write_HDF5SummarizedExperiment <- function(x,
                                            rds_path="se.rds",
                                            h5_path="assays.h5",
                                            chunkdim=NULL, level=NULL,
                                            verbose=FALSE)
{
    .load_SummarizedExperiment_package()

    if (!is(x, "SummarizedExperiment"))
        stop("'x' must be a SummarizedExperiment object")

    if (!isSingleString(rds_path) || rds_path == "")
        stop(wmsg("'rds_path' must be a a non-empty string ",
                  "specifying the path to the RDS file ",
                  "where to write the ", class(x), " object"))

    if (!isSingleString(h5_path) || h5_path == "")
        stop(wmsg("'h5_path' must be a a non-empty string ",
                  "specifying the path to the HDF5 file ",
                  "where to write the assays of the ", class(x), " object"))

    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")

    x@assays <- .write_h5_assays(x@assays, h5_path, chunkdim, level, verbose)
    ans <- x
    x@assays <- .shorten_h5_paths(x@assays)
    saveRDS(x, file=rds_path)
    invisible(ans)
}

.ckeck_and_fix_HDF5Array_assay <- function(x, i, dir)
{
    .load_SummarizedExperiment_package()

    stopifnot(is(x, "SummarizedExperiment"))
    a <- assay(x, i, withDimnames=FALSE)
    if (!is(a, "HDF5Array"))
        stop(wmsg("all assays in the SummarizedExperiment object ",
                  "to load are expected to be HDF5Array objects"))

    ## Check 'a@seed@filepath' and 'a@seed@name'.
    h5_path <- file.path(dir, a@seed@filepath)
    if (!file.exists(h5_path))
        stop(wmsg("assay ", i, " in the SummarizedExperiment object to ",
                  "load points to an HDF5 file that does not exist: ",
                  h5_path))
    h5_content <- try(h5ls(h5_path), silent=TRUE)
    if (inherits(h5_content, "try-error"))
        stop(wmsg("assay ", i, " in the SummarizedExperiment object to ",
                  "load points to an invalid HDF5 file: ", h5_path))
    if (!(a@seed@name %in% h5_content[ , "name"]))
        stop(wmsg("assay ", i, " in the SummarizedExperiment object to ",
                  "load points to an HDF5 dataset (", a@seed@name, ") ",
                  "that does not exist in HDF5 file: ", h5_path))

    ## Fix 'a@seed@filepath'.
    a@seed@filepath <- file_path_as_absolute(h5_path)
    assay(x, i, withDimnames=FALSE) <- a
    x
}

.read_HDF5SummarizedExperiment <- function(filepath)
{
    .load_SummarizedExperiment_package()

    ans <- updateObject(readRDS(filepath), check=FALSE)
    if (!is(ans, "SummarizedExperiment"))
        stop(wmsg("the object serialized in \"", filepath, "\" is not ",
                  "a SummarizedExperiment object or derivative"))
    dir <- dirname(filepath)
    for (i in seq_along(assays(ans)))
        ans <- .ckeck_and_fix_HDF5Array_assay(ans, i, dir)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### saveHDF5SummarizedExperiment() / loadHDF5SummarizedExperiment()
###

.create_dir <- function(dir)
{
    if (file.exists(dir))
        stop(wmsg("\"", dir, "\" already exists and is a file, ",
                  "not a directory"))
    if (!suppressWarnings(dir.create(dir)))
        stop("cannot create directory \"", dir, "\"")
}

.replace_dir <- function(dir, replace)
{
    if (!replace)
        stop(wmsg("Directory \"", dir, "\" already exists. ",
                  "Use 'replace=TRUE' to replace it. ",
                  "Its content will be lost!"))
    if (unlink(dir, recursive=TRUE) != 0L)
        stop("failed to delete directory \"", dir, "\"")
    if (!suppressWarnings(dir.create(dir)))
        stop("cannot create directory \"", dir, "\"")
}

.check_and_delete_files <- function(rds_path, h5_path, replace)
{
    if (dir.exists(rds_path) || dir.exists(h5_path))
        stop(wmsg("\"", rds_path, "\" and/or \"", h5_path, "\" ",
                  "already exist and are directories, not files"))
    if (!(file.exists(rds_path) || file.exists(h5_path)))
        return()
    if (!replace)
        stop(wmsg("Files \"", rds_path, "\" and/or \"", h5_path, "\" ",
                  "already exist. Use a different 'prefix' or use ",
                  "'replace=TRUE' if you really want to replace them."))
    if (unlink(rds_path, recursive=TRUE) != 0L)
        stop("failed to delete file \"", rds_path, "\"")
    if (unlink(h5_path, recursive=TRUE) != 0L)
        stop("failed to delete file \"", h5_path, "\"")
}

### Save all the assays in HDF5 format, including in-memory assays.
### Delayed assays with delayed operations on them are realized while they
### are written to disk..
saveHDF5SummarizedExperiment <- function(x, dir="my_h5_se", prefix="",
                                         replace=FALSE,
                                         chunkdim=NULL, level=NULL,
                                         verbose=FALSE)
{
    .load_SummarizedExperiment_package()

    if (!is(x, "SummarizedExperiment"))
        stop("'x' must be a SummarizedExperiment object")

    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(x),
                  " object (the directory will be created if needed)"))

    if (!isSingleString(prefix))
        stop(wmsg("'prefix' must be a single string"))

    if (!isTRUEorFALSE(replace))
        stop(wmsg("'replace' must be TRUE or FALSE"))

    if (!dir.exists(dir)) {
        .create_dir(dir)
    } else if (prefix == "") {
        .replace_dir(dir, replace)
    }
    rds_path <- file.path(dir, paste0(prefix, "se.rds"))
    h5_path <- file.path(dir, paste0(prefix, "assays.h5"))
    if (prefix != "")
        .check_and_delete_files(rds_path, h5_path, replace)

    .write_HDF5SummarizedExperiment(x, rds_path=rds_path,
                                       h5_path=h5_path,
                                       chunkdim=chunkdim, level=level,
                                       verbose=verbose)
}

.THE_EXPECTED_STUFF <- c(
    "an HDF5-based SummarizedExperiment object ",
    "previously saved with saveHDF5SummarizedExperiment",
    "()"
)

.stop_if_bad_dir <- function(dir, prefix)
{
    if (prefix == "") {
        msg <- c("directory \"", dir, "\" does not seem to contain ",
                 .THE_EXPECTED_STUFF)
    } else {
        msg <- c("Directory \"", dir, "\" does not seem to contain ",
                 head(.THE_EXPECTED_STUFF, n=-1L),
                 "(..., prefix=\"", prefix, "\"). ",
                 "Make sure you're using the same 'prefix' ",
                 "that was used when the object was saved.")
    }
    stop(wmsg(msg))
}

### Does a lot of checking and tries to fail graciously if the content
### of 'dir' doesn't look as expected.
loadHDF5SummarizedExperiment <- function(dir="my_h5_se", prefix="")
{
    .load_SummarizedExperiment_package()

    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))

    if (!isSingleString(prefix))
        stop(wmsg("'prefix' must be a single string"))

    if (!dir.exists(dir)) {
        if (file.exists(dir))
            stop(wmsg("\"", dir, "\" is a file, not a directory"))
        stop(wmsg("directory \"", dir, "\" not found"))
    }

    rds_path <- file.path(dir, paste0(prefix, "se.rds"))
    if (!file.exists(rds_path))
        .stop_if_bad_dir(dir, prefix)

    ans <- try(.read_HDF5SummarizedExperiment(rds_path), silent=TRUE)
    if (inherits(ans, "try-error"))
        .stop_if_bad_dir(dir, prefix)
    ans
}

