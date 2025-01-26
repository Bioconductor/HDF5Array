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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### deparse_html_tree()
###
### Generate HTML from a nested list representation of the HTML document.
###
### HTML element: named ordinary list with 1 to 4 components:
###   1. tag:     single string
###   2. attribs: named character or numeric vector
###   3. style:   unnamed character vector
###   4. content: can be either
###      - a character vector: interpreted as text (including unparsed html);
###      - a named list: must represent an HTML element;
###      - an unnamed list: represents mix content where each list
###        element must be either a character vector or an HTML element.
### Only the first element (tag) is mandatory.
### Example:
###   td_elt <- list(tag="td", style="padding: 2pt", content=c("hi", "there"))
###   tr_elt <- list(tag="tr", content=list(td_elt, td_elt, td_elt))
###   table_elt <- list(tag="table", style="background: grey", content=tr_elt)
### Note that 'table_elt' is a tree structure similar to the Document Object
### Model (DOM) representation, but with a simple representation based on
### nested lists.

### Returns a single string.
.deparse_elt_attribs <- function(attribs)
{
    if (!(is.character(attribs) || is.numeric(attribs)))
        stop(wmsg("'attribs' must be a named character or numeric vector"))
    attribs_names <- names(attribs)
    if (is.null(attribs_names))
        stop(wmsg("'attribs' must be a named character or numeric vector"))
    attribs <- paste0(attribs_names, "=\"", attribs, "\"")
    paste(attribs, collapse=" ")
}

### Returns a single string.
.deparse_elt_style <- function(style)
{
    if (!is.character(style))
        stop(wmsg("'style' must be a character vector"))
    if (!is.null(names(style)))
        stop(wmsg("'style' cannot have names"))
    paste0("style=\"", paste(style, collapse="; "), "\"")
}

### Returns a character vector.
.deparse_elt_content <- function(content)
{
    if (is.character(content))
        return(content)
    if (!is.list(content))
        stop(wmsg("'content' must be either a character vector or a list"))
    if (!is.null(names(content)))
        return(.deparse_elt(content))
    unlist(lapply(content, .deparse_elt_content))
}

### Returns a character vector.
.deparse_elt <- function(elt)
{
    stopifnot(is.list(elt))
    elt_names <- names(elt)
    VALID_NAMES <- c("tag", "attribs", "style", "content")
    invalid_names <- setdiff(elt_names, VALID_NAMES)
    if (length(invalid_names) != 0L) {
        in1string <- paste(invalid_names, collapse=", ")
        stop(wmsg("invalid names on HTML element: ", in1string))
    }
    tag <- elt$tag
    if (is.null(tag))
        stop(wmsg("'tag' missing on HTML element"))
    if (!isSingleString(tag) || tag == "")
        stop(wmsg("'tag' must be a single string"))
    attribs <- elt$attribs
    if (!is.null(attribs))
        attribs <- .deparse_elt_attribs(attribs)
    style <- elt$style
    if (!is.null(style))
        style <- .deparse_elt_style(style)
    content <- elt$content
    if (!is.null(content))
        content <- paste0("  ", .deparse_elt_content(content))
    tag_html <- paste0("<", tag)
    if (!is.null(attribs))
        tag_html <- paste(tag_html, attribs)
    if (!is.null(style))
        tag_html <- paste(tag_html, style)
    tag_html <- paste0(tag_html, ">")
    c(tag_html, content, paste0("</", tag, ">"))
}

deparse_html_tree <- function(html_tree) .deparse_elt_content(html_tree)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .build_html_table()
###

.NGENES_BEFORE_NORM <- 27998
.TABLE_STYLE <- c("margin-left: 0pt",
                  "text-align: center",
                  "font-size: smaller")
.BASE_STYLE <- c("border: 1pt solid #888", "padding: 2pt")
.TH_STYLE <- c(.BASE_STYLE, "background: #CCC")
.MIN_STYLE <- c(.BASE_STYLE, "background: #EFE")
.NA_STYLE <- c(.BASE_STYLE, "color: #D00")

.make_header_tr_elts <- function(block_sizes=c(40L, 100L, 250L))
{
    ## 1st <tr> element.
    th11_elt <- list(tag="th")
    colspan <- 1L + 2L * length(block_sizes)
    th12_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.TH_STYLE,
                     content="sparse<br/>(TENxMatrix)")
    th13_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.TH_STYLE,
                     content="dense<br/>(HDF5Matrix)")
    tr1_elt <- list(tag="tr", content=list(th11_elt, th12_elt, th13_elt))

    ## 2nd <tr> element.
    th2od_elt <- list(tag="th",
                      attribs=c(rowspan=2),
                      style=.TH_STYLE,
                      content=c("object<br />",
                                "dimensions<br />",
                                "(genes&nbsp;x&nbsp;cells)"))
    th2on_elt <- list(tag="th",
                      attribs=c(rowspan=2),
                      style=.TH_STYLE,
                      content="object name")
    th2bs_elts <- lapply(block_sizes,
        function(block_size) {
            content <- sprintf("block&nbsp;size<br />= %s&nbsp;Mb", block_size)
            list(tag="th",
                 attribs=c(colspan=2),
                 style=.TH_STYLE,
                 content=content)
        })
    content <- list(th2od_elt, th2on_elt, th2bs_elts, th2on_elt, th2bs_elts)
    tr2_elt <- list(tag="tr", content=content)

    ## 3rd <tr> element.
    content <- lapply(seq_len(2L * length(block_sizes)),
        function(j) {
            content <- "time<br />in<br />seconds"
            th31_elt <- list(tag="th", style=.TH_STYLE, content=content)
            content <- "max.<br />mem.<br />used"
            th32_elt <- list(tag="th", style=.TH_STYLE, content=content)
            list(th31_elt, th32_elt)
        })
    tr3_elt <- list(tag="tr", content=content)

    list(tr1_elt, tr2_elt, tr3_elt)
}

## Produces 2 * length(times) <td> elements.
.make_td_group <- function(times)
{
    stopifnot(is.integer(times))
    min_time <- suppressWarnings(min(times, na.rm=TRUE))
    lapply(times,
        function(t) {
            if (is.na(t)) {
                style <- .NA_STYLE
            } else if (t == min_time) {
                style <- .MIN_STYLE
            } else {
                style <- .BASE_STYLE
            }
            td1_elt <- list(tag="td", style=style, content=as.character(t))
            td2_elt <- list(tag="td", style=.BASE_STYLE)
            list(td1_elt, td2_elt)
        })
}

## Produces a <tr> element with 3 + 2 * (n1 + n2) <td> elements in it,
## where n1 = length(sparse_times) and n2 = length(dense_times).
.make_data_line <- function(sparse_times, dense_times,
                            step=c("norm", "realize", "pca"),
                            ncells, num_var_genes, dataset_rank)
{
    stopifnot(is.integer(sparse_times), is.integer(dense_times))
    step <- match.arg(step)
    ngenes <- if (step == "norm") .NGENES_BEFORE_NORM else num_var_genes
    content <- sprintf("%d&nbsp;x&nbsp;%d", ngenes, ncells)
    td0_elt <- list(tag="td", style=.BASE_STYLE, content=content)

    ## Results for sparse objects.
    object_name <- sprintf("sparse%d", dataset_rank)
    if (step != "norm")
        object_name <- paste0(object_name, "n")
    code_elt <- list(tag="code", content=object_name)
    td1_elt <- list(tag="td", style=.BASE_STYLE, content=code_elt)
    td_group1 <- .make_td_group(sparse_times)

    ## Results for dense objects.
    object_name <- sprintf("dense%d", dataset_rank)
    if (step != "norm")
        object_name <- paste0(object_name, "n")
    code_elt <- list(tag="code", content=object_name)
    td2_elt <- list(tag="td", style=.BASE_STYLE, content=code_elt)
    td_group2 <- .make_td_group(dense_times)

    content <- list(td0_elt, td1_elt, td_group1, td2_elt, td_group2)
    list(tag="tr", content=content)
}

.make_data_lines <- function(timings, num_var_genes=1000L,
                             step=c("norm", "realize", "pca"),
                             block_sizes=c(40L, 100L, 250L))
{
    step <- match.arg(step)
    unique_ncells <- sort(as.integer(unique(timings[ , "ncells"])))
    lapply(seq_along(unique_ncells),
        function(i) {
            ncells <- unique_ncells[[i]]
            sparse_times <- vapply(block_sizes,
                function(block_size) {
                    .get_time(timings, ncells, num_var_genes, "sparse",
                              step, block_size)
                }, integer(1), USE.NAMES=FALSE)
            dense_times <- vapply(block_sizes,
                function(block_size) {
                    .get_time(timings, ncells, num_var_genes, "dense",
                              step, block_size)
                }, integer(1), USE.NAMES=FALSE)
            .make_data_line(sparse_times, dense_times, step=step,
                            ncells=ncells, num_var_genes=num_var_genes,
                            dataset_rank=i)
        })
}

### Generates an HTML table with 3 + 4 * length(block_sizes) columns.
.build_html_table <- function(timings, num_var_genes=1000L,
                              block_sizes=c(40L, 100L, 250L))
{
    timings <- .check_and_add_missing_cols(timings)
    table_ncols <- 3L + 4L * length(block_sizes)
    header_tr_elts <- .make_header_tr_elts(block_sizes=block_sizes)
    th_style <- c(.BASE_STYLE, "background: #EEE")

    ## 1. Normalization.
    content <- c("1.&nbsp;Normalization (&amp;&nbsp;selection&nbsp;of&nbsp;",
                 num_var_genes, "&nbsp;most&nbsp;variable&nbsp;genes)")
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content=content)
    norm_tr_elts <- .make_data_lines(timings, num_var_genes=num_var_genes,
                                     step="norm", block_sizes=block_sizes)
    norm_tr_elts <- list(list(tag="tr", content=th_elt), norm_tr_elts)

    ## 2. Realization.
    content <- "2.&nbsp;On-disk realization of the normalized datasets"
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content=content)
    realize_tr_elts <- .make_data_lines(timings, num_var_genes=num_var_genes,
                                        step="realize", block_sizes=block_sizes)
    realize_tr_elts <- list(list(tag="tr", content=th_elt), realize_tr_elts)

    ## 3. PCA.
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content="3.&nbsp;PCA")
    pca_tr_elts <- .make_data_lines(timings, num_var_genes=num_var_genes,
                                    step="pca", block_sizes=block_sizes)
    pca_tr_elts <- list(list(tag="tr", content=th_elt), pca_tr_elts)

    content <- list(header_tr_elts, norm_tr_elts, realize_tr_elts, pca_tr_elts)
    list(tag="table",
         style=.TABLE_STYLE,
	 content=content)
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
                               block_sizes=c(40L, 100L, 250L),
                               file="")
{
    file_path <- .find_timings_file(machine_name)
    stopifnot(is.numeric(block_sizes))

    timings <- read.dcf(file_path)  # character matrix
    html_table <- .build_html_table(timings, num_var_genes=num_var_genes,
                                    block_sizes=block_sizes)
    cat(deparse_html_tree(html_table), sep="\n", file=file)
}

