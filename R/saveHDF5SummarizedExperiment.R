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

.SE_RDS_BASENAME <- "se.rds"
.ASSAYS_H5_BASENAME <- "assays.h5"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level manipulation of the assay-to-HDF5 links
###

### Return an error if any of the seeds in any of the assays is not an
### HDF5ArraySeed object.
.get_unique_assay2h5_links <- function(assays)
{
    h5_paths <- lapply(seq_along(assays),
        function(i) {
            a <- assays[[i]]
            ok <- unlist(seedApply(a, is, "HDF5ArraySeed"))
            if (!all(ok))
                stop(wmsg("assay ", i, " in the SummarizedExperiment ",
                          "object to load is not HDF5-based"))
            unique(unlist(seedApply(a, path)))
        })
    unique(unlist(h5_paths))
}

### Assume that all the assays are fully HDF5-based (i.e. that all their
### seeds are HDF5ArraySeed objects). This means that all the seeds should
### have a 'filepath' slot. This is NOT checked!
### Note that shortening the path stored in the HDF5ArraySeed objects breaks
### these objects so we must use direct slot access instead of the path()
### setter to do this. This is because the latter is intended for the end user
### so it makes sure that the path replacement is not breaking the object.
.shorten_assay2h5_links <- function(assays)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        assays[[i]] <- modify_seeds(assays[[i]],
            function(x) {
                x@filepath <- basename(x@filepath)
                x
            })
    }
    assays
}

### Check that all the assays are fully HDF5-based (i.e. that all their seeds
### are HDF5ArraySeed objects), and, via validate_HDF5ArraySeed_dataset(),
### that all the HDF5ArraySeed objects point to HDF5 datasets that are
### accessible and "as expected".
### Restore all the file paths to their absolute canonical form.
.restore_absolute_assay2h5_links <- function(assays, dir)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        assays[[i]] <- modify_seeds(assays[[i]],
            function(x) {
                if (!is(x, "HDF5ArraySeed"))
                    stop(wmsg("assay ", i, " in the SummarizedExperiment ",
                              "object to load is not HDF5-based"))
                h5_path <- file.path(dir, x@filepath)
                ## file_path_as_absolute() will fail if the file does
                ## not exist.
                if (!file.exists(h5_path))
                    stop(wmsg("assay ", i, " in the SummarizedExperiment ",
                              "object to load points to an HDF5 file ",
                              "that does not exist: ", h5_path))
                x@filepath <- file_path_as_absolute(h5_path)
                ## Check that 'x' points to an HDF5 dataset that is accessible
                ## and "as expected".
                msg <- validate_HDF5ArraySeed_dataset(x)
                if (!isTRUE(msg))
                    stop(wmsg("assay ", i, " in the SummarizedExperiment ",
                              "object to load ", msg))
                x
            })
    }
    assays
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .write_HDF5SummarizedExperiment() / .read_HDF5SummarizedExperiment()
###

.serialize_HDF5SummarizedExperiment <- function(x, rds_path, verbose)
{
    x@assays <- .shorten_assay2h5_links(x@assays)
    if (verbose)
        message("Serialize ", class(x), " object to ",
                ifelse(file.exists(rds_path), "existing ", ""),
                "RDS file:\n  ", rds_path)
    saveRDS(x, file=rds_path)
}

.write_h5_assays <- function(assays, h5_path, chunkdim, level, verbose)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        h5_name <- sprintf("assay%03d", i)
        if (verbose)
            message("Start writing assay ", i, "/", nassay, " to ",
                    "HDF5 file:\n  ", h5_path)
        a <- writeHDF5Array(a, h5_path, h5_name, chunkdim, level,
                            verbose=verbose)
        if (verbose)
            message("Finished writing assay ", i, "/", nassay, " to ",
                    "HDF5 file:\n  ", h5_path, "\n")
        assays[[i]] <- a
    }
    assays
}

.write_HDF5SummarizedExperiment <- function(x,
                                            rds_path=.SE_RDS_BASENAME,
                                            h5_path=.ASSAYS_H5_BASENAME,
                                            chunkdim=NULL, level=NULL,
                                            verbose=FALSE)
{
    .load_SummarizedExperiment_package()

    if (!is(x, "SummarizedExperiment"))
        stop(wmsg("'x' must be a SummarizedExperiment object"))

    if (!isSingleString(rds_path) || rds_path == "")
        stop(wmsg("'rds_path' must be a a non-empty string ",
                  "specifying the path to the RDS file ",
                  "where to write the ", class(x), " object"))

    if (!isSingleString(h5_path) || h5_path == "")
        stop(wmsg("'h5_path' must be a a non-empty string ",
                  "specifying the path to the HDF5 file ",
                  "where to write the assays of the ", class(x), " object"))

    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    x@assays <- .write_h5_assays(x@assays, h5_path, chunkdim, level, verbose)
    .serialize_HDF5SummarizedExperiment(x, rds_path, verbose)
    invisible(x)
}

### Does a lot of checking (via .restore_absolute_assay2h5_links()) on
### the assays of the SummarizedExperiment object found in 'rds_path' and
### fails with an informative error message if they don't look as expected.
.read_HDF5SummarizedExperiment <- function(rds_path)
{
    .load_SummarizedExperiment_package()

    if (!file.exists(rds_path))
        stop(wmsg("file not found: ", rds_path))
    if (dir.exists(rds_path))
        stop(wmsg("'", rds_path, "' is a directory, not a file"))

    ans <- updateObject(readRDS(rds_path), check=FALSE)
    if (!is(ans, "SummarizedExperiment"))
        stop(wmsg("the object serialized in \"", rds_path, "\" is not ",
                  "a SummarizedExperiment object or derivative"))
    dir <- dirname(rds_path)
    ans@assays <- .restore_absolute_assay2h5_links(ans@assays, dir)
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
        stop(wmsg("cannot create directory \"", dir, "\""))
}

.replace_dir <- function(dir, replace)
{
    if (!replace)
        stop(wmsg("Directory \"", dir, "\" already exists. ",
                  "Use 'replace=TRUE' to replace it. ",
                  "Its content will be lost!"))
    if (unlink(dir, recursive=TRUE) != 0L)
        stop(wmsg("failed to delete directory \"", dir, "\""))
    if (!suppressWarnings(dir.create(dir)))
        stop(wmsg("cannot create directory \"", dir, "\""))
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
        stop(wmsg("failed to delete file \"", rds_path, "\""))
    if (unlink(h5_path, recursive=TRUE) != 0L)
        stop(wmsg("failed to delete file \"", h5_path, "\""))
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
        stop(wmsg("'x' must be a SummarizedExperiment object"))

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
    rds_path <- file.path(dir, paste0(prefix, .SE_RDS_BASENAME))
    h5_path <- file.path(dir, paste0(prefix, .ASSAYS_H5_BASENAME))
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

### Does a lot of checking (via .read_HDF5SummarizedExperiment()) and fails
### graciously if the content of 'dir' doesn't look as expected.
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

    rds_path <- file.path(dir, paste0(prefix, .SE_RDS_BASENAME))
    ans <- try(.read_HDF5SummarizedExperiment(rds_path), silent=TRUE)
    if (inherits(ans, "try-error"))
        .stop_if_bad_dir(dir, prefix)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### quickResaveHDF5SummarizedExperiment()
###

.stop_if_cannot_quick_resave <- function()
    stop(wmsg("cannot quick-resave a SummarizedExperiment ",
              "object that was not previously saved with ",
              "saveHDF5SummarizedExperiment()"))

.map_h5_path_to_rds_path <- function(h5_path)
{
    dir <- dirname(h5_path)
    h5_basename <- basename(h5_path)
    h5_suffix <- substr(h5_basename,
                        nchar(h5_basename) - nchar(.ASSAYS_H5_BASENAME) + 1L,
                        nchar(h5_basename))
    if (h5_suffix != .ASSAYS_H5_BASENAME)
        .stop_if_cannot_quick_resave()
    prefix <- substr(h5_basename,
                     1L, nchar(h5_basename) - nchar(.ASSAYS_H5_BASENAME))
    rds_basename <- paste0(prefix, .SE_RDS_BASENAME)
    file.path(dir, rds_basename)
}

### Check that the assays of the SummarizedExperiment object found
### in 'rds_path' are HDF5-based and point to 'h5_path'.
.check_HDF5SummarizedExperiment_h5_file <- function(rds_path, h5_path)
{
    se <- try(.read_HDF5SummarizedExperiment(rds_path), silent=TRUE)
    if (inherits(se, "try-error"))
        .stop_if_cannot_quick_resave()
    se_h5_path <-  try(.get_unique_assay2h5_links(se@assays), silent=TRUE)
    if (inherits(se_h5_path, "try-error") ||
        length(se_h5_path) != 1L ||
        se_h5_path != h5_path)
        .stop_if_cannot_quick_resave()
}

### 'x' must have been previously saved with saveHDF5SummarizedExperiment()
### and possibly modified since then.
### A quick-resave preserves the current HDF5 file and datasets and
### re-serializes the SummarizedExperiment object without realizing the
### delayed operations possibly carried by the assays.
quickResaveHDF5SummarizedExperiment <- function(x, verbose=FALSE)
{
    .load_SummarizedExperiment_package()

    if (!is(x, "SummarizedExperiment"))
        stop(wmsg("'x' must be a SummarizedExperiment object"))

    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))

    h5_path <-  try(.get_unique_assay2h5_links(x@assays), silent=TRUE)
    if (inherits(h5_path, "try-error") || length(h5_path) != 1L)
        .stop_if_cannot_quick_resave()
    if (verbose)
        message("All assay data already in HDF5 file:\n  ", h5_path)

    rds_path <- .map_h5_path_to_rds_path(h5_path)
    .check_HDF5SummarizedExperiment_h5_file(rds_path, h5_path)
    .serialize_HDF5SummarizedExperiment(x, rds_path, verbose)
    invisible(x)
}

