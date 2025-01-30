# To run this R script:
#
#   Rscript normalize_and_PCA.R <ncells> <num_var_genes> <format> \
#                     <norm_block_size> <realize_block_size> <pca_block_size>
#
# To run it in "batch mode":
#
#   Rscript normalize_and_PCA.R 12500 1000 s \
#                     250 40 100 >normalize_and_PCA.log 2>&1 &
#

suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(RSpectra))
process_utils_path <- system.file(package="HDF5Array",
                                  "scripts", "process_utils.R", mustWork=TRUE)
source(process_utils_path)
pid <- Sys.getpid()
process_info_log <- tempfile()

## Retrieve and check script arguments.

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 6L)
ncells <- as.integer(args[[1L]])
num_var_genes <- as.integer(args[[2L]])
format <- args[[3L]]
norm_block_size <- as.integer(args[[4L]])     # block size in Mb (normalization)
realize_block_size <- as.integer(args[[5L]])  # block size in Mb (realization)
pca_block_size <- as.integer(args[[6L]])      # block size in Mb (PCA)

stopifnot(isSingleInteger(ncells), ncells > 0L,
          isSingleInteger(num_var_genes), num_var_genes > 0L,
          format %in% c("s", "D", "Ds"),
          isSingleInteger(norm_block_size), norm_block_size > 0L,
          isSingleInteger(realize_block_size), realize_block_size > 0L,
          isSingleInteger(pca_block_size), pca_block_size > 0L)

cat("ncells = ", ncells, "\n", sep="")
cat("num_var_genes = ", num_var_genes, "\n", sep="")
cat("format = ", format, "\n", sep="")
cat("norm_block_size = ", norm_block_size, "\n", sep="")
cat("realize_block_size = ", realize_block_size, "\n", sep="")
cat("pca_block_size = ", pca_block_size, "\n", sep="")
cat("\n")

## Prepare dataset.

hub <- ExperimentHub()
brain_s_path <- suppressMessages(hub[["EH1039"]])
brain_s <- TENxMatrix(brain_s_path, group="mm10")
stopifnot(is_sparse(brain_s),
          identical(chunkdim(brain_s), c(27998L, 1L)))

if (format == "s") {
    full_brain <- brain_s
} else {
    as_sparse <- format == "Ds"
    brain_D_path <- suppressMessages(hub[["EH1040"]])
    full_brain <- HDF5Array(brain_D_path, name="counts", as.sparse=as_sparse)
    stopifnot(identical(chunkdim(full_brain), c(100L, 100L)),
              is_sparse(full_brain) == as_sparse)
    ## The dense HDF5 file does not contain the dimnames of the matrix
    ## so we manually add them:
    dimnames(full_brain) <- dimnames(brain_s)
}
stopifnot(identical(dim(full_brain), c(27998L, 1306127L)))
dataset <- full_brain[ , seq_len(ncells)]

## Define functions simple_normalize() and simple_PCA().

simple_normalize <- function(mat, num_var_genes=1000)
{
    stopifnot(length(dim(mat)) == 2, !is.null(rownames(mat)))
    mat <- mat[rowSums(mat) > 0, ]
    mat <- t(t(mat) * 10000 / colSums(mat))
    row_vars <- rowVars(mat)
    rv_order <- order(row_vars, decreasing=TRUE)
    variable_idx <- head(rv_order, n=num_var_genes)
    mat <- log1p(mat[variable_idx, ])
    mat / rowSds(mat)
}

simple_PCA <- function(mat, k=25)
{
    stopifnot(length(dim(mat)) == 2)
    row_means <- rowMeans(mat)
    Ax <- function(x, args)
        (as.numeric(mat %*% x) - row_means * sum(x))
    Atx <- function(x, args)
        (as.numeric(x %*% mat) - as.vector(row_means %*% x))
    RSpectra::svds(Ax, Atrans=Atx, k=k, dim=dim(mat))
}

## Normalization.

cat("Running normalization ...\n")
DelayedArray::setAutoBlockSize(norm_block_size * 1e6)
loop_pid <- start_log_process_info(pid, process_info_log)
on.exit(stop_log_process_info(loop_pid))
timing <- system.time(normalized <- simple_normalize(dataset, num_var_genes=num_var_genes))
stop_log_process_info(loop_pid)
norm_max_mem_used <- extract_max_mem_used(process_info_log, pid)
gc()
norm_time <- timing[["elapsed"]]
cat("---> normalization completed in ", norm_time, " s.\n\n", sep="")

## On-disk realization of normalized dataset.

cat("On-disk realization of normalized dataset ...\n")
DelayedArray::setAutoBlockSize(realize_block_size * 1e6)
normalized_path <- tempfile()
loop_pid <- start_log_process_info(pid, process_info_log)
on.exit(stop_log_process_info(loop_pid))
if (format == "s") {
    timing <- system.time(
        normalized <- writeTENxMatrix(normalized, normalized_path,
                                      group="matrix", level=0)
    )
} else {
    timing <- system.time(
        normalized <- writeHDF5Array(normalized, normalized_path,
                                     name="normalized_counts", level=0)
    )
}
stop_log_process_info(loop_pid)
realize_max_mem_used <- extract_max_mem_used(process_info_log, pid)
gc()
realize_time <- timing[["elapsed"]]
cat("---> realization completed in ", realize_time, " s.\n\n", sep="")

## PCA.

cat("Running PCA ...\n")
DelayedArray::setAutoBlockSize(pca_block_size * 1e6)
loop_pid <- start_log_process_info(pid, process_info_log)
on.exit(stop_log_process_info(loop_pid))
timing <- system.time(pca <- simple_PCA(normalized))
stop_log_process_info(loop_pid)
pca_max_mem_used <- extract_max_mem_used(process_info_log, pid)
gc()
pca_time <- timing[["elapsed"]]
cat("---> PCA completed in ", pca_time, " s.\n\n", sep="")

cat("ncells: ", ncells, "\n",
    "num_var_genes: ", num_var_genes, "\n",
    "format: ", format, "\n",
    "norm_block_size: ", norm_block_size, "\n",
    "norm_time: ", norm_time, "\n",
    "norm_max_vsz: ", norm_max_mem_used[["max_vsz"]], "\n",
    "norm_max_rss: ", norm_max_mem_used[["max_rss"]], "\n",
    "realize_block_size: ", realize_block_size, "\n",
    "realize_time: ", realize_time, "\n",
    "realize_max_vsz: ", realize_max_mem_used[["max_vsz"]], "\n",
    "realize_max_rss: ", realize_max_mem_used[["max_rss"]], "\n",
    "pca_block_size: ", pca_block_size, "\n",
    "pca_time: ", pca_time, "\n",
    "pca_max_vsz: ", pca_max_mem_used[["max_vsz"]], "\n",
    "pca_max_rss: ", pca_max_mem_used[["max_rss"]], "\n",
    "\n", sep="", file="timings.dcf", append=TRUE)

