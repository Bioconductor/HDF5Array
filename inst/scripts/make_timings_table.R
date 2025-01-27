.EXPECTED_TIMINGS_COLS <- c("ncells", "num_var_genes", "format",
                            "norm_block_size", "norm_time",
                            "realize_block_size", "realize_time",
                            "pca_block_size", "pca_time")

.VALID_FORMATS <- c("sparse", "dense")
.VALID_STEPS <- c("norm", "realize", "pca")

.check_and_add_missing_timings_cols <- function(timings)
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

### Returns a single integer or NA_integer_.
.get_time <- function(timings, ncells, num_var_genes, format, block_size, step)
{
    stopifnot(is.matrix(timings), is.character(timings),
              isSingleString(ncells), isSingleString(num_var_genes),
              isSingleString(format), isSingleString(step),
              isSingleString(block_size))
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

### Returns a 5D integer array.
.fold_timings_matrix_into_5D_array <- function(timings)
{
    timings <- .check_and_add_missing_timings_cols(timings)
    stopifnot(all(timings[ , "format"] %in% .VALID_FORMATS))
    block_size_colnames <- paste0(.VALID_STEPS, "_block_size")
    unique_block_sizes <- as.integer(timings[ , block_size_colnames])
    unique_block_sizes <- sort(unique(unique_block_sizes))
    unique_num_var_genes <- as.integer(timings[ , "num_var_genes"])
    unique_num_var_genes <- sort(unique(unique_num_var_genes))
    unique_ncells <- as.integer(timings[ , "ncells"])
    unique_ncells <- sort(unique(unique_ncells))
    ans_dimnames <- list(step=.VALID_STEPS,
                         block_size=as.character(unique_block_sizes),
                         format=.VALID_FORMATS,
                         num_var_genes=as.character(unique_num_var_genes),
                         ncells=as.character(unique_ncells))
    ans_dim <- lengths(ans_dimnames)
    ans <- array(NA_integer_, dim=ans_dim, dimnames=ans_dimnames)
    for (ncells in dimnames(ans)[[5L]]) {
        for (num_var_genes in dimnames(ans)[[4L]]) {
            for (format in dimnames(ans)[[3L]]) {
                for (block_size in dimnames(ans)[[2L]]) {
                    for (step in dimnames(ans)[[1L]]) {
                      t <- .get_time(timings, ncells, num_var_genes,
                                              format, block_size, step)
                      ans[step, block_size, format, num_var_genes, ncells] <- t
                    }
                }
            }
        }
    }
    ans
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
### .make_td_group()
###

.BASE_STYLE <- c("border: 1pt solid #888", "padding: 2pt")
.TH_STYLE <- c(.BASE_STYLE, "background: #CCC")

.make_td_style <- function(t, min_time, base_style=NULL)
{
    style <- if (is.null(base_style)) .BASE_STYLE else base_style
    if (is.na(t))
        return(c(style, "color: #D00"))
    if (t != min_time)
        return(style)
    if (is.null(base_style)) {
        xstyle <- "background: #EFE"
    } else {
        xstyle <- "font-weight: bold"
    }
    c(style, xstyle)
}

### Produces 2 * length(times) <td> elements.
.make_td_group <- function(times, base_style=NULL)
{
    stopifnot(is.integer(times))
    min_time <- suppressWarnings(min(times, na.rm=TRUE))
    lapply(unname(times),
        function(t) {
            style <- .make_td_style(t, min_time, base_style=base_style)
            td1_elt <- list(tag="td", style=style, content=as.character(t))
	    style <- if (is.null(base_style)) .BASE_STYLE else base_style
            td2_elt <- list(tag="td", style=style)
            list(td1_elt, td2_elt)
        })
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_table_OLD()
###

.NGENES_BEFORE_NORM <- 27998
.TABLE_STYLE <- c("margin-left: 0pt",
                  "text-align: center",
                  "font-size: smaller")

.make_table_header_OLD <- function(block_sizes)
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
    th2bs_elts <- lapply(unname(block_sizes),
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

### Produces a <tr> element with 3 + 2 * (n1 + n2) <td> elements in it,
### where n1 = length(sparse_times) and n2 = length(dense_times).
.make_data_line_OLD <- function(sparse_times, dense_times,
                                step, ncells, num_var_genes, dataset_rank)
{
    stopifnot(is.integer(sparse_times), is.integer(dense_times))
    ngenes <- if (step == "norm") .NGENES_BEFORE_NORM else num_var_genes
    content <- sprintf("%s&nbsp;x&nbsp;%s", ngenes, ncells)
    td0_elt <- list(tag="td", style=.BASE_STYLE, content=content)

    ## Results for sparse objects.
    object_name <- sprintf("sparse%d", dataset_rank)
    if (step != "norm")
        object_name <- paste0(object_name, "n")
    content <- sprintf("<code>%s</code>", object_name)
    td1_elt <- list(tag="td", style=.BASE_STYLE, content=content)
    td_group1 <- .make_td_group(sparse_times)

    ## Results for dense objects.
    object_name <- sprintf("dense%d", dataset_rank)
    if (step != "norm")
        object_name <- paste0(object_name, "n")
    content <- sprintf("<code>%s</code>", object_name)
    td2_elt <- list(tag="td", style=.BASE_STYLE, content=content)
    td_group2 <- .make_td_group(dense_times)

    content <- list(td0_elt, td1_elt, td_group1, td2_elt, td_group2)
    list(tag="tr", content=content)
}

.make_data_lines_OLD <- function(timings, step, num_var_genes="1000")
{
    stopifnot(isSingleString(step), step %in% .VALID_STEPS)
    unique_ncells <- dimnames(timings)$ncells
    lapply(seq_along(unique_ncells),
        function(i) {
            ncells <- unique_ncells[[i]]
            sparse_times <- timings[step, , "sparse", num_var_genes, ncells]
            dense_times  <- timings[step, , "dense" , num_var_genes, ncells]
            .make_data_line_OLD(sparse_times, dense_times,
                                step, ncells, num_var_genes, dataset_rank=i)
        })
}

### Generates an HTML table with 3 + 4 * length(unique_block_sizes) columns,
### where 'unique_block_sizes' is 'dimnames(timings)$block_size'.
.make_table_OLD <- function(timings, num_var_genes="1000")
{
    stopifnot(length(dim(timings)) == 5L, isSingleString(num_var_genes))
    unique_block_sizes <- dimnames(timings)$block_size
    header_tr_elts <- .make_table_header_OLD(unique_block_sizes)
    table_ncols <- 3L + 4L * length(unique_block_sizes)
    th_style <- c(.BASE_STYLE, "background: #EEE")

    ## 1. Normalization.
    content <- c("1.&nbsp;Normalization (&amp;&nbsp;selection&nbsp;of&nbsp;",
                 num_var_genes, "&nbsp;most&nbsp;variable&nbsp;genes)")
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content=content)
    norm_tr_elts <- .make_data_lines_OLD(timings, "norm",
                                         num_var_genes=num_var_genes)
    norm_tr_elts <- list(list(tag="tr", content=th_elt), norm_tr_elts)

    ## 2. Realization.
    content <- "2.&nbsp;On-disk realization of the normalized datasets"
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content=content)
    realize_tr_elts <- .make_data_lines_OLD(timings, "realize",
                                            num_var_genes=num_var_genes)
    realize_tr_elts <- list(list(tag="tr", content=th_elt), realize_tr_elts)

    ## 3. PCA.
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content="3.&nbsp;PCA")
    pca_tr_elts <- .make_data_lines_OLD(timings, "pca",
                                        num_var_genes=num_var_genes)
    pca_tr_elts <- list(list(tag="tr", content=th_elt), pca_tr_elts)

    content <- list(header_tr_elts, norm_tr_elts, realize_tr_elts, pca_tr_elts)
    list(tag="table",
         style=.TABLE_STYLE,
         content=content)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_table()
###

.NORM_TH_STYLE <- c(.BASE_STYLE, "background: #CDC")
.NORM_TD_STYLE <- c(.BASE_STYLE, "background: #EFE")
.REALIZE_TH_STYLE <- c(.BASE_STYLE, "background: #CDD")
.REALIZE_TD_STYLE <- c(.BASE_STYLE, "background: #EFF")
.PCA_TH_STYLE <- c(.BASE_STYLE, "background: #DCC")
.PCA_TD_STYLE <- c(.BASE_STYLE, "background: #FEE")

### Produces 2 <tr> elements with 4 + 6*n <td> elements in each, where
### n = length(block_sizes).
.make_top_header <- function(block_sizes)
{
    ## 1st <tr> element.
    content <- "dimensions<br />of test<br />dataset"
    th1a_elt <- list(tag="th",
                     attribs=c(rowspan=2, colspan=2),
                     style=.TH_STYLE,
                     content=content)
    content <- "dimensions<br />of<br />normalized<br />object"
    th1b_elt <- list(tag="th",
                     attribs=c(rowspan=2, colspan=2),
                     style=.TH_STYLE,
                     content=content)
    th1bs_elts <- lapply(unname(block_sizes),
        function(block_size) {
            content <- sprintf("block&nbsp;size<br />= %s&nbsp;Mb", block_size)
            list(tag="th",
                 attribs=c(colspan=2),
                 style=.TH_STYLE,
                 content=content)
        })
    content <- list(th1a_elt, th1bs_elts, th1b_elt, th1bs_elts, th1bs_elts)
    tr1_elt <- list(tag="tr", content=content)

    ## 2nd <tr> element.
    content <- lapply(seq_len(3L * length(block_sizes)),
        function(j) {
            content <- "time<br />in<br />seconds"
            th21_elt <- list(tag="th", style=.TH_STYLE, content=content)
            content <- "max.<br />mem.<br />used"
            th22_elt <- list(tag="th", style=.TH_STYLE, content=content)
            list(th21_elt, th22_elt)
        })
    tr2_elt <- list(tag="tr", content=content)

    list(tr1_elt, tr2_elt)
}

### Produces a <tr> element with 4 + 6 * num_block_sizes <td> elements in it.
.make_step_header <- function(num_block_sizes, num_var_genes)
{
    content <- "nrow<br />(#&nbsp;genes)"
    th1a_elt <- list(tag="th", style=.TH_STYLE, content=content)
    content <- "nrow<br />(#&nbsp;sel.<br />genes)"
    th1b_elt <- list(tag="th", style=.TH_STYLE, content=content)
    content <- "ncol<br />(#&nbsp;cells)"
    th2_elt <- list(tag="th", style=.TH_STYLE, content=content)

    colspan <- 2L * num_block_sizes
    content <- c("1.&nbsp;NORMALIZATION<br />",
		 "(&amp;&nbsp;selection&nbsp;of&nbsp;", num_var_genes,
                 "&nbsp;most&nbsp;variable&nbsp;genes)")
    N_th_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.NORM_TH_STYLE,
                     content=content)
    R_th_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.REALIZE_TH_STYLE,
                     content="2.&nbsp;ON-DISK&nbsp;REALIZATION")
    P_th_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.PCA_TH_STYLE,
                     content="3.&nbsp;PCA")
    content <- list(th1a_elt, th2_elt, N_th_elt,
                    th1b_elt, th2_elt, R_th_elt, P_th_elt)
    list(tag="tr", content=content)
}

### Produces a <tr> element with 4 + 2 * (n1 + n2 + n3) <td> elements in it,
### where n1 = length(Ntimes), n2 = length(Rtimes), and n3 = length(Ptimes).
.make_data_line <- function(Ntimes, Rtimes, Ptimes, ncells, num_var_genes)
{
    stopifnot(is.integer(Ntimes), is.integer(Rtimes), is.integer(Ptimes))

    td1a_elt <- list(tag="td",
                     style=.BASE_STYLE,
                     content=as.character(.NGENES_BEFORE_NORM))
    td1b_elt <- list(tag="td",
                     style=.BASE_STYLE,
                     content=as.character(num_var_genes))
    td2_elt <- list(tag="td",
                    style=.BASE_STYLE,
                    content=as.character(ncells))

    ## Normalization results.
    td_groupN <- .make_td_group(Ntimes, base_style=.NORM_TD_STYLE)

    ## Realization results.
    td_groupR <- .make_td_group(Rtimes, base_style=.REALIZE_TD_STYLE)

    ## PCA results.
    td_groupP <- .make_td_group(Ptimes, base_style=.PCA_TD_STYLE)

    content <- list(td1a_elt, td2_elt, td_groupN,
                    td1b_elt, td2_elt, td_groupR, td_groupP)
    list(tag="tr", content=content)
}

.make_data_lines <- function(timings, format, num_var_genes)
{
    stopifnot(isSingleString(format), format %in% .VALID_FORMATS)
    unique_ncells <- dimnames(timings)$ncells
    lapply(unique_ncells,
        function(ncells) {
            Ntimes <- timings["norm",    , format, num_var_genes, ncells]
            Rtimes <- timings["realize", , format, num_var_genes, ncells]
            Ptimes <- timings["pca",     , format, num_var_genes, ncells]
            .make_data_line(Ntimes, Rtimes, Ptimes, ncells, num_var_genes)
        })
}

.make_table_section <- function(timings, num_block_sizes, num_var_genes)
{
    stopifnot(isSingleString(num_var_genes))
    step_tr_elt <- .make_step_header(num_block_sizes, num_var_genes)
    table_ncols <- 4L + 6L * num_block_sizes
    th_style <- c(.BASE_STYLE, "background: #EEE")

    ## 1. Sparse representation.
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content="Sparse representation (TENxMatrix object)")
    sparse_tr_elts <- .make_data_lines(timings, "sparse", num_var_genes)
    sparse_tr_elts <- list(list(tag="tr", content=th_elt), sparse_tr_elts)

    ## 2. Dense representation.
    th_elt <- list(tag="th",
                   attribs=c(colspan=table_ncols),
                   style=th_style,
                   content="Dense representation (HDF5Matrix object)")
    dense_tr_elts <- .make_data_lines(timings, "dense", num_var_genes)
    dense_tr_elts <- list(list(tag="tr", content=th_elt), dense_tr_elts)

    list(step_tr_elt, sparse_tr_elts, dense_tr_elts)
}

.make_table <- function(timings)
{
    stopifnot(length(dim(timings)) == 5L)
    unique_block_sizes <- dimnames(timings)$block_size
    num_block_sizes <- length(unique_block_sizes)
    top_header <- .make_top_header(unique_block_sizes)
    section1 <- .make_table_section(timings, num_block_sizes,
                                    num_var_genes="1000")
    section2 <- .make_table_section(timings, num_block_sizes,
                                    num_var_genes="2000")
    content <- list(top_header, section1, section2)
    list(tag="table",
         style=.TABLE_STYLE,
         content=content)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_timings_table_OLD()
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

make_timings_table_OLD <- function(machine_name, num_var_genes="1000", file="")
{
    stopifnot(isSingleString(num_var_genes))
    file_path <- .find_timings_file(machine_name)
    timings <- read.dcf(file_path)  # character matrix
    timings <- .fold_timings_matrix_into_5D_array(timings)
    table_elt <- .make_table_OLD(timings, num_var_genes=num_var_genes)
    cat(deparse_html_tree(table_elt), sep="\n", file=file)
}

make_timings_table <- function(machine_name, file="")
{
    file_path <- .find_timings_file(machine_name)
    timings <- read.dcf(file_path)  # character matrix
    timings <- .fold_timings_matrix_into_5D_array(timings)
    table_elt <- .make_table(timings)
    cat(deparse_html_tree(table_elt), sep="\n", file=file)
}

