.prefixes <- c("_block_size", "_time", "_max_vsz", "_max_rss")
.TIMINGS_DB_COLS <- c("ncells", "num_var_genes", "format",
                      paste0("norm", .prefixes),
                      paste0("realize", .prefixes),
                      paste0("pca", .prefixes))

.VALID_FORMATS <- c("s", "D", "Ds")
.VALID_STEPS <- c("norm", "realize", "pca")

.check_and_add_missing_timings_db_cols <- function(timings_db)
{
    stopifnot(is.matrix(timings_db), is.character(timings_db))

    ## Only for compatibility with old timings db files.
    format <- timings_db[ , "format"]
    format[format %in% "sparse"] <- "s"
    format[format %in% "dense"]  <- "D"
    timings_db[ , "format"] <- format

    missing_cols <- setdiff(.TIMINGS_DB_COLS, colnames(timings_db))
    if (length(missing_cols) != 0L) {
        ## Add missing cols (filled with NAs).
        m <- matrix(NA_character_,
                    nrow=nrow(timings_db), ncol=length(missing_cols),
                    dimnames=list(NULL, missing_cols))
        timings_db <- cbind(timings_db, m)
    }
    timings_db <- timings_db[ , .TIMINGS_DB_COLS, drop=FALSE]
    na_idx <- which(is.na(timings_db[ , "num_var_genes"]))
    if (length(na_idx) != 0L)
        timings_db[na_idx, "num_var_genes"] <- 1000
    timings_db
}

### Returns a single integer or NA_integer_.
.get_val_from_timings_db <-
    function(timings_db, varname=c("time", "max_vsz", "max_rss"),
             ncells, num_var_genes, format, block_size, step)
{
    stopifnot(is.matrix(timings_db), is.character(timings_db),
              isSingleString(ncells), isSingleString(num_var_genes),
              isSingleString(format), isSingleString(step),
              isSingleString(block_size))
    varname <- match.arg(varname)
    val_colname <- paste0(step, "_", varname)
    ok1 <- timings_db[ , "ncells"] == ncells &
           timings_db[ , "num_var_genes"] == num_var_genes &
           timings_db[ , "format"] == format
    block_size_colname <- paste0(step, "_block_size")
    ok2 <- timings_db[ , block_size_colname] == block_size
    rowidx <- which(ok1 & ok2)
    if (length(rowidx) == 0L)
        return(NA_integer_)
    if (length(rowidx) != 1L)
        stop(wmsg("no (or more than one) \"", val_colname, "\" value found for",
                  "ncells=", ncells, ", num_var_genes=", num_var_genes, ", ",
                  "format=\"", format, "\", step=\"", step, "\", ",
                  "and block_size=", block_size))
    val <- suppressWarnings(as.numeric(timings_db[rowidx, val_colname]))
    as.integer(val + 0.5)  # rounding to the closest integer
}

### Returns a 5D integer array.
.extract_var_from_timings_db <-
    function(timings_db, varname=c("time", "max_vsz", "max_rss"))
{
    stopifnot(is.matrix(timings_db),
              identical(colnames(timings_db), .TIMINGS_DB_COLS))
    varname <- match.arg(varname)
    stopifnot(all(timings_db[ , "format"] %in% .VALID_FORMATS))
    block_size_colnames <- paste0(.VALID_STEPS, "_block_size")
    unique_block_sizes <- as.integer(timings_db[ , block_size_colnames])
    unique_block_sizes <- sort(unique(unique_block_sizes))
    unique_num_var_genes <- as.integer(timings_db[ , "num_var_genes"])
    unique_num_var_genes <- sort(unique(unique_num_var_genes))
    unique_ncells <- as.integer(timings_db[ , "ncells"])
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
                val <- .get_val_from_timings_db(timings_db, varname,
                                                ncells, num_var_genes,
                                                format, block_size, step)
                ans[step, block_size, format, num_var_genes, ncells] <- val
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
    closing_tag <- paste0("</", tag, ">")
    tag <- paste0("<", tag)
    attribs <- elt$attribs
    if (!is.null(attribs)) {
        attribs <- .deparse_elt_attribs(attribs)
        tag <- paste(tag, attribs)
    }
    style <- elt$style
    if (!is.null(style)) {
        style <- .deparse_elt_style(style)
        tag <- paste(tag, style)
    }
    tag <- paste0(tag, ">")
    content <- elt$content
    if (is.null(content))
        return(paste0(tag, closing_tag))
    content <- paste0("  ", .deparse_elt_content(content))
    c(tag, content, closing_tag)
}

deparse_html_tree <- function(html_tree) .deparse_elt_content(html_tree)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_td_group()
###

.MEM_THRESHOLD <- 4
.LIGHT_RED <- "#D77"  # to display memory usage that is NA or > .MEM_THRESHOLD
.BASE_STYLE <- c("border: 1pt solid #BBB", "padding: 2pt")

.make_time_td_style <- function(t, min_time, base_style=NULL)
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

.make_mem_td_style <- function(m, base_style=NULL)
{
    style <- if (is.null(base_style)) .BASE_STYLE else base_style
    #style <- c(style, "font-style: italic")
    ## Use light red if NA or > .MEM_THRESHOLD, otherwise light grey.
    color <- if (is.na(m) || m > .MEM_THRESHOLD) .LIGHT_RED else "#777"
    c(style, paste0("color: ", color))
}

### Produces 2 * length(times) <td> elements.
.make_td_group <- function(times, mem, base_style=NULL, draw_box=FALSE)
{
    stopifnot(is.integer(times), is.integer(mem),
              length(times) == length(mem))
    min_time <- suppressWarnings(min(times, na.rm=TRUE))
    lapply(seq_along(times),
        function(i) {
            t <- times[[i]]
            m <- mem[[i]] / 1024  # from Mb to Gb
            style <- .make_time_td_style(t, min_time, base_style=base_style)
            content <- as.character(t)
            if (draw_box && !is.na(t) && t == min_time) {
                span_style <- "border: 1pt solid black"
                content <- sprintf("<span style=\"%s\">&nbsp;%s&nbsp;</span>",
                                   span_style, content)
            }
            td1_elt <- list(tag="td", style=style, content=content)
            style <- .make_mem_td_style(m, base_style=base_style)
            content <- sprintf("%.1f", m)  # max. mem. used in Gb
            if (!is.na(m)) {
                Gb <- "Gb"
                if (m <= .MEM_THRESHOLD)
                    Gb <- sprintf("<span style=\"color: %s\">%s</span>",
                                  "#AAA", Gb)
                content <- paste0(content, Gb)
            }
            td2_elt <- list(tag="td", style=style, content=content)
            list(td1_elt, td2_elt)
        })
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .make_table()
###

.NGENES_BEFORE_NORM <- 27998
.TABLE_STYLE <- c("border-spacing: 0px",
                  "border-collapse: collapse",
                  "margin-left: 0pt",
                  "text-align: center",
                  "font-size: smaller")

.make_hline <- function(colspan, height="0pt", color="#BBB")
{
    style <- c(paste0("border: 1pt solid ", color), "padding: 0px",
               paste("height:", height), paste("background:", color))
    td_elt <- list(tag="td", attribs=c(colspan=colspan), style=style)
    list(tag="tr", content=td_elt)
}

.make_footnote <- function(colspan, title=NULL)
{
    style <- "font-style: italic"
    formats <- sprintf("<span style=\"font-weight: bold\">[%s]</span>",
                       .VALID_FORMATS)
    fmt_explained <- sprintf("%s&nbsp;%s", formats,
                             c("TENxMatrix&nbsp;(sparse)",
                               "HDF5Matrix&nbsp;(dense)",
                               "HDF5Matrix&nbsp;as&nbsp;sparse"))
    #fmt_explained <- paste0(paste(fmt_explained, collapse="; "), ".")
    fmt_explained <- paste0(paste(fmt_explained, collapse=" &mdash; "), ".")

    ## Replace "four" with whatever is the new number of block sizes
    ## if we ever happen to change that.
    content <- c(
        "Formats:&nbsp;", fmt_explained, "<br />",
        "For each operation, the best time across the ",
        "four different block sizes is displayed in ",
        "<span style=\"font-weight: bold\">bold</span>.<br />",
        "In addition, if it's also the best time across the three formats ",
        "(", formats[[1L]], ",", formats[[2L]], ", and", formats[[3L]], "), ",
        "then we <span style=\"font-weight: bold; border: 1pt solid black\">",
        "&nbsp;box&nbsp;</span> it ",
        "(only for Normalization and PCA).<br />",
        "The \"max. mem. used\" is the max RSS (Resident Set Size) ",
        "value obtained by running <code>ps u -p &lt;PID&gt;</code> ",
        "every second while performing a given operation.<br />",
        "\"max. mem. used\" values > ", .MEM_THRESHOLD, "Gb are displayed ",
        "in <span style=\"color: ", .LIGHT_RED, "\">light red</span>.")
    if (!is.null(title)) {
        title <- sprintf("<span style=\"font-weight: bold\">%s</span><br />",
                         title)
        content <- c(title, content)
    }
    td_elt <- list(tag="td",
                   attribs=c(colspan=colspan),
                   style=style,
                   content=content)
    list(tag="tr", content=td_elt)
}

.TH_BASE_STYLE <- c(.BASE_STYLE, "color: #555")
.TH_STYLE <- c(.TH_BASE_STYLE, "background: #CCC")
.TH_LIGHTER_STYLE <- c(.TH_BASE_STYLE, "background: #E6E6E6")

.NORM_TH_STYLE <- c(.TH_BASE_STYLE, "background: #C7CFC7")
.NORM_TH_LIGHTER_STYLE <- c(.TH_BASE_STYLE, "background: #E7EFE7")
.NORM_TD_STYLE <- c(.BASE_STYLE, "background: #F7FFF7")
.NORM_TD_DENSE_STYLE <- c(.BASE_STYLE, "background: #F0F8F0")

.REALIZE_TH_STYLE <- c(.TH_BASE_STYLE, "background: #C7CFCF")
.REALIZE_TH_LIGHTER_STYLE <- c(.TH_BASE_STYLE, "background: #E7EFEF")
.REALIZE_TD_STYLE <- c(.BASE_STYLE, "background: #F7FFFF")
.REALIZE_TD_DENSE_STYLE <- c(.BASE_STYLE, "background: #F0F8F8")

.PCA_TH_STYLE <- c(.TH_BASE_STYLE, "background: #CFC7C7")
.PCA_TH_LIGHTER_STYLE <- c(.TH_BASE_STYLE, "background: #EFE7E7")
.PCA_TD_STYLE <- c(.BASE_STYLE, "background: #FFF7F7")
.PCA_TD_DENSE_STYLE <- c(.BASE_STYLE, "background: #F8F0F0")

### Produces 2 <tr> elements that span 3 + 6 * n columns each, where
### n = length(block_sizes).
.make_top_header <- function(block_sizes)
{
    ## 1st <tr> element.
    content <- "F<br />o<br />r<br />m<br />a<br />t"
    th1a_elt <- list(tag="th",
                     attribs=c(rowspan=2),
                     style=c(.TH_STYLE, "font-size: smaller"),
                     content=content)
    content <- "Test&nbsp;Dataset"
    th1b_elt <- list(tag="th",
                     style=.TH_STYLE,
                     content=content)
    content <- "Normalized<br />Test&nbsp;Dataset"
    th1c_elt <- list(tag="th",
                     style=.TH_STYLE,
                     content=content)
    make_th1_elts <- function(block_sizes, style) {
        lapply(unname(block_sizes),
            function(block_size) {
                content <- sprintf("block&nbsp;size<br />= %s&nbsp;Mb",
                                   block_size)
                list(tag="th",
                     attribs=c(colspan=2),
                     style=style,
                     content=content)
            })
    }
    N_th1_elts <- make_th1_elts(block_sizes, .NORM_TH_STYLE)
    R_th1_elts <- make_th1_elts(block_sizes, .REALIZE_TH_STYLE)
    P_th1_elts <- make_th1_elts(block_sizes, .PCA_TH_STYLE)
    content <- list(th1a_elt, th1b_elt, N_th1_elts,
                    th1c_elt, R_th1_elts, P_th1_elts)
    tr1_elt <- list(tag="tr", content=content)

    ## 2nd <tr> element.
    content <- c("nrow&nbsp;x&nbsp;ncol<br />",
                 "(#&nbsp;genes&nbsp;x&nbsp;#&nbsp;cells)")
    th2a_elt <- list(tag="th",
                     style=.TH_STYLE,
                     content=content)
    content <- c("nrow&nbsp;x&nbsp;ncol<br />",
                 "(#&nbsp;sel.&nbsp;genes<br />x&nbsp;#&nbsp;cells)")
    th2c_elt <- list(tag="th",
                     style=.TH_STYLE,
                     content=content)
    make_th2_elts <- function(block_sizes, style) {
        lapply(seq_along(block_sizes),
            function(j) {
                content <- "time<br />in<br />sec."
                th21_elt <- list(tag="th", style=style, content=content)
                content <- "max.<br />mem.<br />used"
                #style <- c(style, "font-style: italic", "color: #777")
                style <- c(style, "color: #777")
                th22_elt <- list(tag="th", style=style, content=content)
                list(th21_elt, th22_elt)
            })
    }
    N_th2_elts <- make_th2_elts(block_sizes, .NORM_TH_STYLE)
    R_th2_elts <- make_th2_elts(block_sizes, .REALIZE_TH_STYLE)
    P_th2_elts <- make_th2_elts(block_sizes, .PCA_TH_STYLE)
    content <- c(list(th2a_elt), N_th2_elts,
                 list(th2c_elt), R_th2_elts, P_th2_elts)
    tr2_elt <- list(tag="tr", style="font-size: smaller", content=content)

    list(tr1_elt, tr2_elt)
}

### Produces a <tr> element that spans 3 + 6 * num_block_sizes columns.
.make_steps_header <- function(num_block_sizes, num_var_genes)
{
    th1_elt <- list(tag="th", style=.TH_LIGHTER_STYLE)
    #th2_elt <- list(tag="th", style=.TH_LIGHTER_STYLE, content="format")
    th2_elt <- list(tag="th", style=.TH_LIGHTER_STYLE)

    colspan <- 2L * num_block_sizes
    content <- c("1.&nbsp;NORMALIZATION<br />",
                 "&amp;&nbsp;selection&nbsp;of&nbsp;", num_var_genes,
                 "&nbsp;most&nbsp;variable&nbsp;genes")
    N_th_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.NORM_TH_LIGHTER_STYLE,
                     content=content)
    content <- c("2.&nbsp;ON-DISK&nbsp;REALIZATION<br />",
                 "of the normalized dataset")
    R_th_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.REALIZE_TH_LIGHTER_STYLE,
                     content=content)
    content <- "3.&nbsp;PCA<br />of the normalized dataset"
    P_th_elt <- list(tag="th",
                     attribs=c(colspan=colspan),
                     style=.PCA_TH_LIGHTER_STYLE,
                     content=content)
    content <- list(th1_elt, th2_elt, N_th_elt, th1_elt, R_th_elt, P_th_elt)
    list(tag="tr", content=content)
}

### Produces a <tr> element that spans 3 + 2 * (n1 + n2 + n3) columns,
### where n1 = length(Ntimes), n2 = length(Rtimes), and n3 = length(Ptimes).
.make_data_line <- function(ncells, format, num_var_genes,
                            Ntimes, Nbox, Nmem,
                            Rtimes, Rbox, Rmem,
                            Ptimes, Pbox, Pmem)
{
    stopifnot(isSingleString(format),
              is.integer(Ntimes), is.integer(Rtimes), is.integer(Ptimes),
              isTRUEorFALSE(Nbox), isTRUEorFALSE(Rbox), isTRUEorFALSE(Pbox),
              is.integer(Nmem), is.integer(Rmem), is.integer(Pmem),
              length(Ntimes) == length(Nmem),
              length(Rtimes) == length(Rmem),
              length(Ptimes) == length(Pmem))

    style <- c(.BASE_STYLE, "font-weight: bold", "color: #888")
    if (format != "s")
        style <- c(style, "background: #F8F8F8")
    td0_elt <- list(tag="td",
                    style=style,
                    content=paste0("[", format, "]"))

    content <- sprintf("<span style=\"%s\">%s&nbsp;x&nbsp;</span>%s",
                       "color: #888", .NGENES_BEFORE_NORM, ncells)
    td1_elt <- list(tag="td",
                    attribs=c(rowspan=3),
                    style=.BASE_STYLE,
                    content=content)

    content <- sprintf("%s<span style=\"%s\">&nbsp;x&nbsp;</span>%s",
                       num_var_genes, "color: #888", ncells)
    td2_elt <- list(tag="td",
                    attribs=c(rowspan=3),
                    style=.BASE_STYLE,
                    content=content)

    ## Normalization results.
    base_style <-
        if (format != "s") .NORM_TD_DENSE_STYLE else .NORM_TD_STYLE
    td_groupN <- .make_td_group(Ntimes, Nmem,
                                base_style=base_style, draw_box=Nbox)

    ## Realization results.
    base_style <-
        if (format != "s") .REALIZE_TD_DENSE_STYLE else .REALIZE_TD_STYLE
    td_groupR <- .make_td_group(Rtimes, Rmem,
                                base_style=base_style, draw_box=Rbox)

    ## PCA results.
    base_style <-
        if (format != "s") .PCA_TD_DENSE_STYLE else .PCA_TD_STYLE
    td_groupP <- .make_td_group(Ptimes, Pmem,
                                base_style=base_style, draw_box=Pbox)

    if (format == "s") {
        content <- list(td0_elt,
                        td1_elt, td_groupN,
                        td2_elt, td_groupR, td_groupP)
    } else {
        content <- list(td0_elt,
                                 td_groupN,
                                 td_groupR, td_groupP)
    }
    list(tag="tr", content=content)
}

### Produces 3 <tr> elements, one for each format in .VALID_FORMATS.
.make_data_line_triplet <- function(times, memused, ncells, num_var_genes)
{
    ## Ntimes, Nmemused: 3-col matrices with colnames .VALID_FORMATS and
    ## 1 row per block size.
    Ntimes <- times["norm", , , num_var_genes, ncells]
    stopifnot(identical(colnames(Ntimes), .VALID_FORMATS))
    Nmemused <- memused["norm", , , num_var_genes, ncells]
    stopifnot(identical(colnames(Nmemused), .VALID_FORMATS))
    Ntimes_min <- suppressWarnings(min(Ntimes, na.rm=TRUE))
    ## Logical vector of length 3 with names .VALID_FORMATS on it.
    Nbox <- vapply(colnames(Ntimes),
        function(j)
            suppressWarnings(min(Ntimes[ , j],  na.rm=TRUE)) == Ntimes_min,
        logical(1))

    ## Rtimes, Rmemused: 3-col matrices with colnames .VALID_FORMATS and
    ## 1 row per block size.
    Rtimes <- times["realize", , , num_var_genes, ncells]
    Rmemused <- memused["realize", , , num_var_genes, ncells]
    Rtimes_min <- suppressWarnings(min(Rtimes, na.rm=TRUE))
    ## Logical vector of length 3 with names .VALID_FORMATS on it.
    Rbox <- vapply(colnames(Rtimes),
        function(j)
            suppressWarnings(min(Rtimes[ , j],  na.rm=TRUE)) == Rtimes_min,
        logical(1))
    ## Disable boxing of the best realization time for now (too many boxes!
    ## which is distracting and we don't really care about realization anyway).
    Rbox[] <- FALSE

    ## Ptimes, Pmemused: 3-col matrices with colnames .VALID_FORMATS and
    ## 1 row per block size.
    Ptimes <- times["pca", , , num_var_genes, ncells]
    Pmemused <- memused["pca", , , num_var_genes, ncells]
    Ptimes_min <- suppressWarnings(min(Ptimes, na.rm=TRUE))
    ## Logical vector of length 3 with names .VALID_FORMATS on it.
    Pbox <- vapply(colnames(Ptimes),
        function(j)
            suppressWarnings(min(Ptimes[ , j],  na.rm=TRUE)) == Ptimes_min,
        logical(1))

    lapply(.VALID_FORMATS,
        function(format)
            .make_data_line(ncells, format, num_var_genes,
                    Ntimes[ , format], Nbox[[format]], Nmemused[ , format],
                    Rtimes[ , format], Rbox[[format]], Rmemused[ , format],
                    Ptimes[ , format], Pbox[[format]], Pmemused[ , format])
    )
}

.make_table_section <- function(times, memused,
                                num_block_sizes, num_var_genes,
                                hline=NULL)
{
    stopifnot(isSingleString(num_var_genes))
    steps_header <- .make_steps_header(num_block_sizes, num_var_genes)
    unique_ncells <- dimnames(times)$ncells
    tr_elts <- lapply(unique_ncells,
        function(ncells) {
            line_triplet <- .make_data_line_triplet(times, memused,
                                                    ncells, num_var_genes)
            if (is.null(hline))
                return(line_triplet)
            c(list(hline), line_triplet)
        })
    section <- list(steps_header, tr_elts)
    if (is.null(hline))
        return(section)
    c(list(hline), section)
}

### times, memused: 5D integer arrays of same dimensions and dimnames.
.make_table <- function(times, memused, title=NULL)
{
    stopifnot(length(dim(times)) == 5L,
              identical(dim(times), dim(memused)),
              identical(dimnames(times), dimnames(memused)))
    unique_block_sizes <- dimnames(times)$block_size
    num_block_sizes <- length(unique_block_sizes)
    top_header <- .make_top_header(unique_block_sizes)
    hline <- .make_hline(3L+6L*num_block_sizes)
    section1 <- .make_table_section(times, memused, num_block_sizes,
                                    num_var_genes="1000", hline=hline)
    section2 <- .make_table_section(times, memused, num_block_sizes,
                                    num_var_genes="2000", hline=hline)
    footnote <- .make_footnote(3L+6L*num_block_sizes, title=title)
    content <- list(top_header, section1, section2, hline, footnote)
    list(tag="table",
         style=.TABLE_STYLE,
         content=content)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### make_timings_table()
###

.find_timings_db_file <- function(machine_name)
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
        stop(wmsg("no timings db files found in '", machine_path, "'"))
    sort(file_paths, decreasing=TRUE)[[1L]]
}

make_timings_table <- function(machine_name, title=NULL, file="")
{
    timings_db_file <- .find_timings_db_file(machine_name)
    timings_db <- read.dcf(timings_db_file)  # character matrix
    timings_db <- .check_and_add_missing_timings_db_cols(timings_db)
    times <- .extract_var_from_timings_db(timings_db, varname="time")
    ## We choose to populate the "max. mem. used" table columns with
    ## the "max_rss" values, not the "max_vsz" values, because the VSZ
    ## as reported by 'ps u -p <PID>' seems meaningless on macOS.
    memused <- .extract_var_from_timings_db(timings_db, varname="max_rss")
    table_elt <- .make_table(times, memused, title)
    cat(deparse_html_tree(table_elt), sep="\n", file=file)
}

