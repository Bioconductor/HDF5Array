.find_timings_file <- function(machine_name)
{
    suppressPackageStartupMessages(library(S4Vectors))
    suppressPackageStartupMessages(library(HDF5Array))
    stopifnot(isSingleString(machine_name))
    machine_path <- system.file(package="HDF5Array",
                                "scripts", "timings_db", machine_name)
    if (machine_path == "")
        stop(wmsg("no '", machine_name, "' folder in timings db"))
    file_path <- file.path(machine_path, "timings.dcf")
    if (file.exists(file_path))
        return(file_path)
    pattern <- "^timings.*\\.dcf$"
    file_paths <- list.files(machine_path, pattern=pattern, full.names=TRUE)
    if (length(file_paths) == 0L)
        stop(wmsg("no timings files found in '", machine_path, "'"))
    sort(file_paths, decreasing=TRUE)[[1L]]
}

.get_time <- function(timings, ncells, format, step, block_size)
{
    stopifnot(is.matrix(timings), is.character(timings))
    ok1 <- timings[ , "ncells"] == ncells & timings[ , "format"] == format
    if (step == "norm") {
        ok2 <- timings[ , "norm_block_size"] == block_size
        t <- timings[ok1 & ok2, "norm_time"]
    } else {
        ok2 <- timings[ , "pca_block_size"] == block_size
        t <- timings[ok1 & ok2, "pca_time"]
    }
    if (length(t) == 0L)
        return("N/A")
    if (length(t) != 1L)
        stop(wmsg("no time (or more than one time) found for ",
                  "ncells=", ncells, ", format=\"", format, "\", ",
                  "step=\"", step, "\", and block_size=", block_size))
    round(as.numeric(t))  # rounding to the second (closest)
}

.NGENES_BEFORE_NORM <- 27998
.NGENES_AFTER_NORM <- 1000

make_timings_table <- function(machine_name,
                               step=c("norm", "pca"),
                               block_sizes=c(40, 100, 250))
{
    file_path <- .find_timings_file(machine_name)
    step <- match.arg(step)
    stopifnot(is.numeric(block_sizes))

    timings <- read.dcf(file_path)  # character matrix
    EXPECTED_COLS <- c("ncells", "format",
                       "norm_block_size", "norm_time",
                       "pca_block_size", "pca_time")
    stopifnot(setequal(colnames(timings), EXPECTED_COLS))

    cat('<table style="margin-left: 0pt; text-align: center; font-size: smaller;">\n')

    ## Header lines

    TH_STYLE <- c("background: #CCC", "border: 1pt solid #888", "padding: 3pt")
    TH_STYLE <- paste0("style='", paste(TH_STYLE, collapse="; "), "'")

    cat('  <tr>\n')
    cat('    <th></th>\n')
    colspan <- 1L + 2L * length(block_sizes)
    cat(sprintf('    <th %s colspan="%d">\n', TH_STYLE, colspan))
    cat('      sparse<br/>(TENxMatrix)\n')
    cat('    </th>\n')
    cat(sprintf('    <th %s colspan="%d">\n', TH_STYLE, colspan))
    cat('      dense<br/>(HDF5Matrix)\n')
    cat('    </th>\n')
    cat('  </tr>\n')

    cat('  <tr>\n')
    cat(sprintf('    <th %s rowspan="2">\n', TH_STYLE))
    cat('      object<br />dimensions<br />(genes&nbsp;x&nbsp;cells)\n')
    cat('    </th>\n')
    cat(sprintf('    <th %s rowspan="2">\n', TH_STYLE))
    cat('      object name\n')
    cat('    </th>\n')
    for (j in seq_along(block_sizes)) {
        cat(sprintf('    <th %s colspan="2">\n', TH_STYLE))
        cat(sprintf('      block&nbsp;size<br />= %s&nbsp;Mb\n',
                    block_sizes[[j]]))
        cat('    </th>\n')
    }
    cat(sprintf('    <th %s rowspan="2">\n', TH_STYLE))
    cat('      object name\n')
    cat('    </th>\n')
    for (j in seq_along(block_sizes)) {
        cat(sprintf('    <th %s colspan="2">\n', TH_STYLE))
        cat(sprintf('      block&nbsp;size<br />= %s&nbsp;Mb\n',
                    block_sizes[[j]]))
        cat('    </th>\n')
    }
    cat('  </tr>\n')

    cat('  <tr>\n')
    for (j in seq_len(2L * length(block_sizes))) {
        cat(sprintf('    <th %s>time<br />in<br />seconds</th>\n', TH_STYLE))
        cat(sprintf('    <th %s>max.<br />mem.<br />used</th>\n', TH_STYLE))
    }
    cat('  </tr>\n')

    ## Data lines

    TD_STYLE <- c("border: 1pt solid #888", "padding: 3pt")
    TD_STYLE <- paste0("style='", paste(TD_STYLE, collapse="; "), "'")

    ngenes <- if (step == "norm") .NGENES_BEFORE_NORM else .NGENES_AFTER_NORM
    unique_ncells <- sort(as.integer(unique(timings[ , "ncells"])))
    for (i in seq_along(unique_ncells)) {
        cat('  <tr>\n')
        ncells <- unique_ncells[[i]]
        cat(sprintf('    <td %s>%d x %d</td>\n', TD_STYLE, ngenes, ncells))
        object_name <- sprintf("sparse%d", i)
        if (step == "pca")
            object_name <- paste0(object_name, "n")
        cat(sprintf('    <td %s><code>%s</code></td>\n', TD_STYLE, object_name))
        for (j in seq_along(block_sizes)) {
            block_size <- block_sizes[[j]]
            t <- .get_time(timings, ncells, "sparse", step, block_size)
            cat(sprintf('    <td %s>%s</td>\n', TD_STYLE, t))
            cat(sprintf('    <td %s></td>\n', TD_STYLE))
        }
        object_name <- sprintf("dense%d", i)
        if (step == "pca")
            object_name <- paste0(object_name, "n")
        cat(sprintf('    <td %s><code>%s</code></td>\n', TD_STYLE, object_name))
        for (j in seq_along(block_sizes)) {
            block_size <- block_sizes[[j]]
            t <- .get_time(timings, ncells, "dense", step, block_size)
            cat(sprintf('    <td %s>%s</td>\n', TD_STYLE, t))
            cat(sprintf('    <td %s></td>\n', TD_STYLE))
        }
        cat('  </tr>\n')
    }
    cat('</table>\n')
}

