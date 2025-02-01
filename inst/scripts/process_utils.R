### Only supported on Linux and Mac.

### interval: sleep time in ps infinite loop.
### Returns the PID of the loop process.
start_log_process_info <- function(pid, logfile, interval=1)
{
    oldwd <- getwd()
    setwd(system.file(package="HDF5Array", "scripts", mustWork=TRUE))
    on.exit(setwd(oldwd))
    cmd <- "./start_log_process_info.sh"
    system2(cmd, args=c(pid, logfile, interval), stdout=TRUE)
}

stop_log_process_info <- function(loop_pid)
{
    ## We ignore the result.
    res <- suppressWarnings(
        system2("kill", args=loop_pid, stdout=TRUE, stderr=TRUE)
    )
}

### Returns a 11-col matrix.
.import_logfile_as_matrix <- function(logfile)
{
    ## 11 fields expected on both supported platforms, in the same order, but
    ## with subtle differences in some of the names.
    LINUX_FIELDS <- c("USER", "PID", "%CPU", "%MEM", "VSZ", "RSS",
                      "TTY", "STAT", "START", "TIME", "COMMAND")
    MAC_FIELDS <- c("USER", "PID", "%CPU", "%MEM", "VSZ", "RSS",
                    "TT", "STAT", "STARTED", "TIME", "COMMAND")
    stopifnot(length(LINUX_FIELDS) == length(MAC_FIELDS))
    lines <- readLines(logfile)
    fragments <- strsplit(lines, " +")
    is_Linux_header <-
        vapply(fragments, function(frags) identical(frags, LINUX_FIELDS),
               logical(1))
    is_Mac_header <-
        vapply(fragments, function(frags) identical(frags, MAC_FIELDS),
               logical(1))
    is_header <- is_Linux_header | is_Mac_header
    data <- fragments[!is_header]
    data <- lapply(data, head, n=length(LINUX_FIELDS))
    if (length(data) == 0L)
        stop(wmsg("File '", logfile, "' not in 'ps u' format"))
    ## Maybe last line got truncated in which case we drop it.
    if (length(data[[length(data)]]) < length(LINUX_FIELDS))
        data <- data[-length(data)]
    if (length(data) == 0L)
        stop(wmsg("File '", logfile, "' not in 'ps u' format"))
    matrix(unlist(data), nrow=length(data), byrow=TRUE)
}

### Returns max memory used in a named integer vector made of 2 elements:
### the max VSZ (Virtual Memory Size) and max RSS (Resident Set Size), both
### reported **in Mb**. The names on the vector are "max_vsz" and "max_rss".
extract_max_mem_used <- function(logfile, pid)
{
    data <- .import_logfile_as_matrix(logfile)
    PID <- data[ , 2L]
    if (!all(PID == pid))
        stop(wmsg("File '", logfile, "' does not contain 'ps u' ",
                  "output for expected process (pid ", pid, ")"))
    VSZ <- as.integer(data[ , 5L])
    RSS <- as.integer(data[ , 6L])
    ans <- c(max_vsz=max(VSZ) , max_rss=max(RSS)) / 1024  # in Mb
    setNames(as.integer(ans + 0.5), names(ans))
}

