\name{TENxMatrix-class}
\docType{class}

\alias{class:TENxMatrix}
\alias{TENxMatrix-class}
\alias{TENxMatrix}

\alias{DelayedArray,TENxMatrixSeed-method}

\alias{sparsity,TENxMatrix-method}
\alias{read_sparse_block,TENxMatrix-method}
\alias{extractNonzeroDataByCol,TENxMatrix-method}

\alias{coerce,TENxMatrix,dgCMatrix-method}
\alias{coerce,TENxMatrix,sparseMatrix-method}

\title{10x Genomics datasets as DelayedMatrix objects}

\description{
  A 10x Genomics dataset like the "1.3 Million Brain Cell Dataset" is an
  HDF5 sparse matrix stored in CSR/CSC/Yale format ("Compressed Sparse Row").

  The TENxMatrix class is a \link[DelayedArray]{DelayedMatrix} subclass
  for representing and operating on this kind of dataset.

  All the operations available for \link[DelayedArray]{DelayedMatrix}
  objects work on TENxMatrix objects.
}

\usage{
## Constructor function:
TENxMatrix(filepath, group="matrix")
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the HDF5 file where the 10x Genomics
    dataset is located.
  }
  \item{group}{
    The name of the group in the HDF5 file containing the 10x Genomics data.
  }
}

\value{
  \code{TENxMatrix()} returns a TENxMatrix object.
}

\details{
  In addition to all the methods defined for
  \link[DelayedArray]{DelayedMatrix} objects, TENxMatrix objects
  support the following specialized methods: \code{sparsity()} and
  \code{extractNonzeroDataByCol()}. See \code{?\link{H5SparseMatrixSeed}}
  for more information about what these methods do.
}

\note{
  If your dataset uses the HDF5 sparse matrix representation from
  10x Genomics, use the \code{TENxMatrix()} constructor documented here.

  But if your dataset uses the conventional (a.k.a. dense) HDF5
  representation, use the \code{\link{HDF5Array}()} constructor instead.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing conventional (a.k.a.
          dense) HDF5 datasets as \link[DelayedArray]{DelayedArray} objects.

    \item \link[DelayedArray]{DelayedMatrix} objects in the \pkg{DelayedArray}
          package.

    \item \code{\link{writeTENxMatrix}} for writing a matrix-like object
          as an HDF5-based sparse matrix.

    \item The \code{\link[TENxBrainData]{TENxBrainData}} dataset (in the
          \pkg{TENxBrainData} package).

    \item \code{\link[parallel]{detectCores}} from the \pkg{parallel}
          package.

    \item \code{\link[DelayedArray]{setAutoBPPARAM}} and
          \code{\link[DelayedArray]{setAutoBlockSize}} in the
          \pkg{DelayedArray} package.

    \item \code{\link[DelayedArray]{colAutoGrid}} and
          \code{\link[DelayedArray]{blockApply}} in the
          \pkg{DelayedArray} package.

    \item The \link{TENxMatrixSeed} helper class.

    \item \code{\link{h5ls}} to list the content of an HDF5 file.

    \item \link[IRanges]{NumericList} and \link[IRanges]{IntegerList}
          objects in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## THE "1.3 Million Brain Cell Dataset" AS A DelayedMatrix OBJECT
## ---------------------------------------------------------------------

## The 1.3 Million Brain Cell Dataset from 10x Genomics is available
## via ExperimentHub:

library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]

## 'fname' is an HDF5 file. Use h5ls() to list its content:
h5ls(fname)

## The 1.3 Million Brain Cell Dataset is represented by the "mm10"
## group. We point the TENxMatrix() constructor to this group to
## create a TENxMatrix object representing the dataset:
oneM <- TENxMatrix(fname, group="mm10")
oneM

is(oneM, "DelayedMatrix")  # TRUE
seed(oneM)
path(oneM)
sparsity(oneM)

## Some examples of delayed operations:
oneM != 0
oneM^2

## ---------------------------------------------------------------------
## SOME EXAMPLES OF ROW/COL SUMMARIZATION
## ---------------------------------------------------------------------

## In order to reduce computation times, we'll use only the first
## 25000 columns of the 1.3 Million Brain Cell Dataset:
oneM25k <- oneM[ , 1:25000]

## Row/col summarization methods like rowSums() use a block-processing
## mechanism behind the scene that can be controlled via global
## settings. 2 important settings that can have a strong impact on
## performance are the automatic number of workers and automatic block
## size, controlled by setAutoBPPARAM() and setAutoBlockSize()
## respectively.
library(BiocParallel)
if (.Platform$OS.type != "windows") {
    ## On a modern Linux laptop with 8 cores (as reported by
    ## parallel::detectCores()) and 16 Gb of RAM, reasonably good
    ## performance is achieved by setting the automatic number of workers
    ## to 5 or 6 and the automatic block size between 300 Mb and 400 Mb:
    workers <- 5
    block_size <- 3e8  # 300 Mb
    setAutoBPPARAM(MulticoreParam(workers))
} else {
    ## MulticoreParam() is not supported on Windows so we use SnowParam()
    ## on this platform. Also we reduce the block size to 200 Mb on
    ## 32-bit Windows to avoid memory allocation problems (they tend to
    ## be common there because a process cannot use more than 3 Gb of
    ## memory).
    workers <- 4
    setAutoBPPARAM(SnowParam(workers))
    block_size <- if (.Platform$r_arch == "i386") 2e8 else 3e8
}
setAutoBlockSize(block_size)

## We're ready to compute the library sizes, number of genes expressed
## per cell, and average expression across cells:
system.time(lib_sizes <- colSums(oneM25k))
system.time(n_exprs <- colSums(oneM25k != 0))
system.time(ave_exprs <- rowMeans(oneM25k))

## Note that the 3 computations above load the data in oneM25k 3 times
## in memory. This can be avoided by computing the 3 summarizations in
## a single pass with blockApply(). First we define the function that
## we're going to apply to each block of data:
FUN <- function(block)
  list(colSums(block), colSums(block != 0), rowSums(block))

## Then we call blockApply() to apply FUN() to each block. The blocks
## are defined by the grid passed to the 'grid' argument. In this case
## we supply a grid made with colAutoGrid() to generate blocks of full
## columns (see ?colAutoGrid for more information):
system.time({
  block_results <- blockApply(oneM25k, FUN, grid=colAutoGrid(oneM25k),
                              verbose=TRUE)
})

## 'block_results' is a list with 1 list element per block in
## colAutoGrid(oneM25k). Each list element is the result that was
## obtained by applying FUN() on the block so is itself a list of
## length 3.
## Let's combine the results:
lib_sizes2 <- unlist(lapply(block_results, `[[`, 1L))
n_exprs2 <- unlist(lapply(block_results, `[[`, 2L))
block_rowsums <- unlist(lapply(block_results, `[[`, 3L), use.names=FALSE)
tot_exprs <- rowSums(matrix(block_rowsums, nrow=nrow(oneM25k)))
ave_exprs2 <- setNames(tot_exprs / ncol(oneM25k), rownames(oneM25k))

## Sanity checks:
stopifnot(all.equal(lib_sizes, lib_sizes2))
stopifnot(all.equal(n_exprs, n_exprs2))
stopifnot(all.equal(ave_exprs, ave_exprs2))

## Turn off parallel evaluation and reset automatic block size to factory
## settings:
setAutoBPPARAM()
setAutoBlockSize()

## ---------------------------------------------------------------------
## extractNonzeroDataByCol()
## ---------------------------------------------------------------------

## extractNonzeroDataByCol() provides a convenient and very efficient
## way to extract the nonzero data in a compact form:
nonzeros <- extractNonzeroDataByCol(oneM, 1:25000)  # takes < 5 sec.

## The data is returned as an IntegerList object with one list element
## per column and no row indices associated to the values in the object.
## Furthermore, the values within a given list element can be returned
## in any order:
nonzeros

names(nonzeros) <- colnames(oneM25k)

## This can be used to compute some simple summaries like the library
## sizes and the number of genes expressed per cell. For these use
## cases, it is a lot more efficient than using colSums(oneM25k) and
## colSums(oneM25k != 0):
lib_sizes3 <- sum(nonzeros)
n_exprs3 <- lengths(nonzeros)

## Sanity checks:
stopifnot(all.equal(lib_sizes, lib_sizes3))
stopifnot(all.equal(n_exprs, n_exprs3))
}
\keyword{classes}
\keyword{methods}
