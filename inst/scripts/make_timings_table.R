### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .build_html_table()
###

.EXPECTED_TIMINGS_COLS <- c("ncells", "num_var_genes", "format",
                            "norm_block_size", "norm_time",
                            "realize_block_size", "realize_time",
                            "pca_block_size", "pca_time")

.check_and_add_missing_cols <- function(timings)
{
    stopifnot(is.matrix(timings))
    missing_cols <- setdiff(.EXPECTED_TIMINGS_COLS, colnames(timings))
    if (length(missing_cols) != 0L) {
        ## Add missing cols (filled with NAs).
        m <- matrix(NA_character_,
                    nrow=nrow(timings), ncol=length(missing_cols),
                    dimnames=list(NULL, missing_cols))
        timings <- cbind(timings, m)
    }
    timings <- timings[ , .EXPECTED_TIMINGS_COLS, drop=FALSE]
    na_idx <- which(is.na(timings[ , "num_var_genes"]))
    if (length(na_idx) != 0L)
        timings[na_idx, "num_var_genes"] <- 1000
    timings
}

.get_time <- function(timings, ncells, num_var_genes, format, step, block_size)
{
    stopifnot(is.matrix(timings), is.character(timings),
              isSingleNumber(ncells), isSingleNumber(num_var_genes),
              isSingleString(format), isSingleString(step),
              isSingleNumber(block_size))
    ok1 <- timings[ , "ncells"] == ncells &
           timings[ , "num_var_genes"] == num_var_genes &
           timings[ , "format"] == format
    block_size_colname <- paste0(step, "_block_size")
    ok2 <- timings[ , block_size_colname] == block_size
    rowidx <- which(ok1 & ok2)
    if (length(rowidx) == 0L)
        return(NA_integer_)
    if (length(rowidx) != 1L)
        stop(wmsg("no time (or more than one time) found for ",
                  "ncells=", ncells, ", num_var_genes=", num_var_genes, ", ",
                  "format=\"", format, "\", step=\"", step, "\", ",
                  "and block_size=", block_size))
    time_colname <- paste0(step, "_time")
    t <- suppressWarnings(as.numeric(timings[rowidx, time_colname]))
    as.integer(t + 0.5)  # rounding to the closest integer
}

.make_header_lines <- function(timings, block_sizes=c(40L, 100L, 250L))
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
.BASE_STYLE <- c("border: 1pt solid #888", "padding: 2pt")
.MIN_STYLE <- c(.BASE_STYLE, "background: #EFE")
.NA_STYLE <- c(.BASE_STYLE, "color: #D00")

## Produces 2 * length(times) td elements.
.make_td_group <- function(times)
{
    stopifnot(is.integer(times))
    base_style <- paste0("style='", paste(.BASE_STYLE, collapse="; "), "'")
    min_style <- paste0("style='", paste(.MIN_STYLE, collapse="; "), "'")
    na_style <- paste0("style='", paste(.NA_STYLE, collapse="; "), "'")

    min_time <- suppressWarnings(min(times, na.rm=TRUE))
    for (j in seq_along(times)) {
        t <- times[[j]]
        if (is.na(t)) {
            style <- na_style
        } else if (t == min_time) {
            style <- min_style
        } else {
            style <- base_style
        }
        cat(sprintf('    <td %s>%d</td>\n', style, t))
        cat(sprintf('    <td %s></td>\n', base_style))
    }
}

## Produces a tr element with 3 + 4 * length(block_sizes) td elements in it.
.make_data_line <- function(timings, step=c("norm", "realize", "pca"),
                            ncells, num_var_genes, dataset_rank,
                            block_sizes=c(40L, 100L, 250L))
{
    step <- match.arg(step)
    ngenes <- if (step == "norm") .NGENES_BEFORE_NORM else num_var_genes

    base_style <- paste0("style='", paste(.BASE_STYLE, collapse="; "), "'")

    cat('  <tr>\n')
    cat(sprintf('    <td %s>%d&nbsp;x&nbsp;%d</td>\n',
                base_style, ngenes, ncells))

    ## Results for sparse objects.
    object_name <- sprintf("sparse%d", dataset_rank)
    if (step != "norm")
        object_name <- paste0(object_name, "n")
    cat(sprintf('    <td %s><code>%s</code></td>\n', base_style, object_name))
    times <- vapply(block_sizes,
        function(block_size)
            .get_time(timings, ncells, num_var_genes, "sparse",
                      step, block_size),
        integer(1), USE.NAMES=FALSE)
    .make_td_group(times)

    ## Results for dense objects.
    object_name <- sprintf("dense%d", dataset_rank)
    if (step != "norm")
        object_name <- paste0(object_name, "n")
    cat(sprintf('    <td %s><code>%s</code></td>\n', base_style, object_name))
    times <- vapply(block_sizes,
        function(block_size)
            .get_time(timings, ncells, num_var_genes, "dense",
                      step, block_size),
        integer(1), USE.NAMES=FALSE)
    .make_td_group(times)

    cat('  </tr>\n')
}

.make_data_lines <- function(timings, num_var_genes=1000L,
                             step=c("norm", "realize", "pca"),
                             block_sizes=c(40L, 100L, 250L))
{
    step <- match.arg(step)
    unique_ncells <- sort(as.integer(unique(timings[ , "ncells"])))
    for (i in seq_along(unique_ncells)) {
        ncells <- unique_ncells[[i]]
        .make_data_line(timings, step=step,
                        ncells=ncells, num_var_genes=num_var_genes,
                        dataset_rank=i, block_sizes=block_sizes)
    }
}

### Generates an HTML table with 3 + 4 * length(block_sizes) columns.
.build_html_table <- function(timings, num_var_genes=1000L,
                              block_sizes=c(40L, 100L, 250L))
{
    timings <- .check_and_add_missing_cols(timings)

    table_ncols <- 3L + 4L * length(block_sizes)

    TABLE_STYLE <- c("margin-left: 0pt",
                     "text-align: center",
                     "font-size: smaller")
    table_style <- paste0("style='", paste(TABLE_STYLE, collapse="; "), "'")
    cat(sprintf('<table %s>\n', table_style))

    .make_header_lines(timings, block_sizes=block_sizes)

    th_style <- c(.BASE_STYLE, "background: #EEE")
    th_style <- paste0("style='", paste(th_style, collapse="; "), "'")

    cat(sprintf('<tr><th %s colspan="%d">1. Normalization</th></tr>\n',
                th_style, table_ncols))
    .make_data_lines(timings, num_var_genes=num_var_genes,
                     step="norm", block_sizes=block_sizes)

    cat(sprintf('<tr><th %s colspan="%d">2. On-disk realization of the normalized datasets</th></tr>\n',
                th_style, table_ncols))
    .make_data_lines(timings, num_var_genes=num_var_genes,
                     step="realize", block_sizes=block_sizes)

    cat(sprintf('<tr><th %s colspan="%d">3. PCA</th></tr>\n',
                th_style, table_ncols))
    .make_data_lines(timings, num_var_genes=num_var_genes,
                     step="pca", block_sizes=block_sizes)

    cat('</table>\n')
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_timings_table()
###

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

make_timings_table <- function(machine_name, num_var_genes=1000L,
                               block_sizes=c(40L, 100L, 250L))
{
    file_path <- .find_timings_file(machine_name)
    stopifnot(is.numeric(block_sizes))

    timings <- read.dcf(file_path)  # character matrix
    .build_html_table(timings, num_var_genes=num_var_genes,
                      block_sizes=block_sizes)
}

