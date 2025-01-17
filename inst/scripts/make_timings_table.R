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
        return(NA_integer_)
    if (length(t) != 1L)
        stop(wmsg("no time (or more than one time) found for ",
                  "ncells=", ncells, ", format=\"", format, "\", ",
                  "step=\"", step, "\", and block_size=", block_size))
    as.integer(as.numeric(t) + 0.5)  # rounding to the second (closest)
}

.make_header_lines <- function(timings, block_sizes=c(40, 100, 250))
{
    th_style <- c(.BASE_STYLE, "background: #CCC")
    th_style <- paste0("style='", paste(th_style, collapse="; "), "'")

    cat('  <tr>\n')
    cat('    <th></th>\n')
    colspan <- 1L + 2L * length(block_sizes)
    cat(sprintf('    <th %s colspan="%d">\n', th_style, colspan))
    cat('      sparse<br/>(TENxMatrix)\n')
    cat('    </th>\n')
    cat(sprintf('    <th %s colspan="%d">\n', th_style, colspan))
    cat('      dense<br/>(HDF5Matrix)\n')
    cat('    </th>\n')
    cat('  </tr>\n')

    cat('  <tr>\n')
    cat(sprintf('    <th %s rowspan="2">\n', th_style))
    cat('      object<br />dimensions<br />(genes&nbsp;x&nbsp;cells)\n')
    cat('    </th>\n')
    cat(sprintf('    <th %s rowspan="2">\n', th_style))
    cat('      object name\n')
    cat('    </th>\n')
    for (j in seq_along(block_sizes)) {
        cat(sprintf('    <th %s colspan="2">\n', th_style))
        cat(sprintf('      block&nbsp;size<br />= %s&nbsp;Mb\n',
                    block_sizes[[j]]))
        cat('    </th>\n')
    }
    cat(sprintf('    <th %s rowspan="2">\n', th_style))
    cat('      object name\n')
    cat('    </th>\n')
    for (j in seq_along(block_sizes)) {
        cat(sprintf('    <th %s colspan="2">\n', th_style))
        cat(sprintf('      block&nbsp;size<br />= %s&nbsp;Mb\n',
                    block_sizes[[j]]))
        cat('    </th>\n')
    }
    cat('  </tr>\n')

    cat('  <tr>\n')
    for (j in seq_len(2L * length(block_sizes))) {
        cat(sprintf('    <th %s>time<br />in<br />seconds</th>\n', th_style))
        cat(sprintf('    <th %s>max.<br />mem.<br />used</th>\n', th_style))
    }
    cat('  </tr>\n')
}

.NGENES_BEFORE_NORM <- 27998
.NGENES_AFTER_NORM <- 1000
.BASE_STYLE <- c("border: 1pt solid #888", "padding: 3pt")

## Produces 2 * length(times) td elements.
.make_td_group <- function(times)
{
    stopifnot(is.integer(times))
    base_style <- paste0("style='", paste(.BASE_STYLE, collapse="; "), "'")
    green_style <- c(.BASE_STYLE, "background: #EFE")
    green_style <- paste0("style='", paste(green_style, collapse="; "), "'")
    min_time <- suppressWarnings(min(times, na.rm=TRUE))
    for (j in seq_along(times)) {
        t <- times[[j]]
        style <- if (!is.na(t) && t == min_time) green_style else base_style
        cat(sprintf('    <td %s>%d</td>\n', style, t))
        cat(sprintf('    <td %s></td>\n', base_style))
    }
}

## Produces a tr element with 3 + 4 * length(block_sizes) td elements in it.
.make_data_line <- function(timings, step=c("norm", "pca"),
                            ncells, dataset_rank, block_sizes=c(40, 100, 250))
{
    step <- match.arg(step)
    ngenes <- if (step == "norm") .NGENES_BEFORE_NORM else .NGENES_AFTER_NORM

    base_style <- paste0("style='", paste(.BASE_STYLE, collapse="; "), "'")

    cat('  <tr>\n')
    cat(sprintf('    <td %s>%d&nbsp;x&nbsp;%d</td>\n',
                base_style, ngenes, ncells))

    ## Results for sparse objects.
    object_name <- sprintf("sparse%d", dataset_rank)
    if (step == "pca")
        object_name <- paste0(object_name, "n")
    cat(sprintf('    <td %s><code>%s</code></td>\n', base_style, object_name))
    times <- vapply(block_sizes,
        function(block_size)
            .get_time(timings, ncells, "sparse", step, block_size),
        integer(1), USE.NAMES=FALSE)
    .make_td_group(times)

    ## Results for dense objects.
    object_name <- sprintf("dense%d", dataset_rank)
    if (step == "pca")
        object_name <- paste0(object_name, "n")
    cat(sprintf('    <td %s><code>%s</code></td>\n', base_style, object_name))
    times <- vapply(block_sizes,
        function(block_size)
            .get_time(timings, ncells, "dense", step, block_size),
        integer(1), USE.NAMES=FALSE)
    .make_td_group(times)

    cat('  </tr>\n')
}

.make_data_lines <- function(timings, step=c("norm", "pca"),
                             block_sizes=c(40, 100, 250))
{
    step <- match.arg(step)
    unique_ncells <- sort(as.integer(unique(timings[ , "ncells"])))
    for (i in seq_along(unique_ncells)) {
        ncells <- unique_ncells[[i]]
        .make_data_line(timings, step=step,
                        ncells=ncells, dataset_rank=i, block_sizes=block_sizes)
    }
}

### Generates an HTML table with 3 + 4 * length(block_sizes) columns.
.make_table <- function(timings, block_sizes=c(40, 100, 250))
{
    table_ncols <- 3L + 4L * length(block_sizes)

    TABLE_STYLE <- c("margin-left: 0pt",
                     "text-align: center",
                     "font-size: smaller")
    table_style <- paste0("style='", paste(TABLE_STYLE, collapse="; "), "'")
    cat(sprintf('<table %s>\n', table_style))

    .make_header_lines(timings, block_sizes=block_sizes)

    th_style <- c(.BASE_STYLE, "background: #EEE")
    th_style <- paste0("style='", paste(th_style, collapse="; "), "'")

    cat(sprintf('<tr><th %s colspan="%d">Normalization</th></tr>\n',
                th_style, table_ncols))
    .make_data_lines(timings, "norm", block_sizes=block_sizes)

    cat(sprintf('<tr><th %s colspan="%d">PCA</th></tr>\n',
                th_style, table_ncols))
    .make_data_lines(timings, "pca", block_sizes=block_sizes)

    cat('</table>\n')
}

make_timings_table <- function(machine_name,
                               block_sizes=c(40, 100, 250))
{
    file_path <- .find_timings_file(machine_name)
    stopifnot(is.numeric(block_sizes))

    timings <- read.dcf(file_path)  # character matrix
    EXPECTED_COLS <- c("ncells", "format",
                       "norm_block_size", "norm_time",
                       "pca_block_size", "pca_time")
    stopifnot(setequal(colnames(timings), EXPECTED_COLS))
    .make_table(timings, block_sizes=block_sizes)
}

