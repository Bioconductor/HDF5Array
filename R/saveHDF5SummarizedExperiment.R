### =========================================================================
### Save/load an HDF5-based SummarizedExperiment object
### -------------------------------------------------------------------------


.create_dir <- function(dir, replace)
{
    if (dir.exists(dir)) {
        if (!replace)
            stop(wmsg("Directory \"", dir, "\" already exists. ",
                      "Use 'replace=TRUE' to replace it. ",
                      "Its content will be lost!"))
        if (unlink(dir, recursive=TRUE) != 0L)
            stop("failed to delete directory \"", dir, "\"")
    } else if (file.exists(dir)) {
        stop(wmsg("\"", dir, "\" already exists and is a file, ",
                  "not a directory"))
    }
    if (!suppressWarnings(dir.create(dir)))
        stop("cannot create directory \"", dir, "\"")
}

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

### Save all the assays in HDF5 format, including in-memory assays.
### Delayed assays with delayed operations on them are realized while they
### are written to disk..
saveHDF5SummarizedExperiment <- function(x, dir="my_h5_se", replace=FALSE,
                                         chunkdim=NULL, level=NULL,
                                         verbose=FALSE)
{
    library(SummarizedExperiment)

    if (!is(x, "SummarizedExperiment"))
        stop("'x' must be a SummarizedExperiment object")
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(x),
                  " object (the directory will be created)"))
    if (!isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")

    .create_dir(dir, replace)

    h5_path <- file.path(dir, "assays.h5")
    x@assays <- .write_h5_assays(x@assays, h5_path, chunkdim, level, verbose)

    rds_path <- file.path(dir, "se.rds")
    ans <- x
    x@assays <- .shorten_h5_paths(x@assays)
    saveRDS(x, file=rds_path)

    invisible(ans)
}

.THE_EXPECTED_STUFF <- c(
    "an HDF5-based SummarizedExperiment object previously ",
    "saved with saveHDF5SummarizedExperiment()"
)

.stop_if_bad_dir <- function(dir)
    stop(wmsg("directory \"", dir, "\" does not seem to contain ",
              .THE_EXPECTED_STUFF))

### Does a lot of checking and tries to fail graciously if the content
### of 'dir' doesn't look as expected.
loadHDF5SummarizedExperiment <- function(dir="my_h5_se")
{
    library(SummarizedExperiment)
    assays <- SummarizedExperiment::assays
    assay <- SummarizedExperiment::assay
    `assay<-` <- SummarizedExperiment::`assay<-`

    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    h5_path <- file.path(dir, "assays.h5")
    rds_path <- file.path(dir, "se.rds")
    if (!file.exists(h5_path) || !file.exists(rds_path))
        .stop_if_bad_dir(dir)
    h5_content <- try(h5ls(h5_path), silent=TRUE)
    if (inherits(h5_content, "try-error"))
        .stop_if_bad_dir(dir)
    h5_datasets <- h5_content[ , "name"]
    ans <- updateObject(readRDS(rds_path), check=FALSE)
    if (!is(ans, "SummarizedExperiment"))
        .stop_if_bad_dir(dir)
    for (i in seq_along(assays(ans))) {
        a <- assay(ans, i, withDimnames=FALSE)
        if (!is(a, "HDF5Array") || !identical(a@seed@filepath, "assays.h5") ||
            !(a@seed@name %in% h5_datasets))
            .stop_if_bad_dir(dir)
        a@seed@filepath <-
            file_path_as_absolute(file.path(dir, a@seed@filepath))
        assay(ans, i, withDimnames=FALSE) <- a
    }
    ans
}

