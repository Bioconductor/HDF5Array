# To run this R script:
#
#   Rscript normalize_and_PCA.R <ncells> <format> \
#                               <norm_block_size> <pca_block_size>
#
# To run it in "batch mode":
#
#   Rscript normalize_and_PCA.R 12500 sparse \
#                               250 100 >normalize_and_PCA.log 2>&1 &
#

suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(RSpectra))

## Retrieve and check script arguments.

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 4L)
ncells <- as.integer(args[[1L]])
format <- args[[2L]]
norm_block_size <- as.integer(args[[3L]])  # block size in Mb (normalization)
pca_block_size <- as.integer(args[[4L]])   # block size in Mb (PCA)

stopifnot(isSingleInteger(ncells), ncells > 0L,
          format %in% c("sparse", "dense"),
          isSingleInteger(norm_block_size), norm_block_size > 0L,
          isSingleInteger(pca_block_size), pca_block_size > 0L)

cat("ncells = ", ncells, "\n", sep="")
cat("format = ", format, "\n", sep="")
cat("norm_block_size = ", norm_block_size, "\n", sep="")
cat("pca_block_size = ", pca_block_size, "\n", sep="")
cat("\n")

## Prepare dataset.

hub <- ExperimentHub()
h5_path <- suppressMessages(hub[["EH1039"]])
full_dataset <- TENxMatrix(h5_path, group="mm10")
stopifnot(is_sparse(full_dataset),
          identical(chunkdim(full_dataset), c(27998L, 1L)))

if (format == "dense") {
    full_sparse_dataset <- full_dataset
    h5_path <- suppressMessages(hub[["EH1040"]])
    full_dataset <- HDF5Array(h5_path, name="counts")
    stopifnot(!is_sparse(full_dataset),
              identical(chunkdim(full_dataset), c(100L, 100L)))
    ## The dense HDF5 file does not contain the dimnames of the matrix
    ## so we manually add them:
    dimnames(full_dataset) <- dimnames(full_sparse_dataset)
}
stopifnot(identical(dim(full_dataset), c(27998L, 1306127L)))
dataset <- full_dataset[ , seq_len(ncells)]

## Define functions simple_normalize() and simple_PCA().

simple_normalize <- function(mat, num_variable_genes=1000)
{
    stopifnot(length(dim(mat)) == 2, !is.null(rownames(mat)))
    mat <- mat[rowSums(mat) > 0, ]
    mat <- t(t(mat) * 10000 / colSums(mat))
    row_vars <- rowVars(mat)
    rv_order <- order(row_vars, decreasing=TRUE)
    variable_idx <- head(rv_order, n=num_variable_genes)
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
timing <- system.time(normalized <- simple_normalize(dataset))
gc()
norm_time <- timing[["elapsed"]]
cat("---> normalization completed in ", norm_time, " s.\n\n", sep="")
normalized_path <- tempfile()
if (format == "sparse") {
    normalized <- writeTENxMatrix(normalized, normalized_path,
                                  group="matrix", level=0)
} else {
    normalized <- writeHDF5Array(normalized, normalized_path,
                                 name="normalized_counts", level=0)
}
gcm <- gc()

## PCA.

cat("Running PCA ...\n")
DelayedArray::setAutoBlockSize(pca_block_size * 1e6)
if (format == "sparse") {
    normalized <- TENxMatrix(normalized_path)
} else {
    normalized <- HDF5Array(normalized_path, name="normalized_counts")
}
timing <- system.time(pca <- simple_PCA(normalized))
gc()
pca_time <- timing[["elapsed"]]
cat("---> PCA completed in ", pca_time, " s.\n\n", sep="")

cat("ncells: ", ncells, "\n",
    "format: ", format, "\n",
    "norm_block_size: ", norm_block_size, "\n",
    "norm_time: ", norm_time, "\n",
    "pca_block_size: ", pca_block_size, "\n",
    "pca_time: ", pca_time, "\n",
    "\n", sep="", file="timings.dcf", append=TRUE)

