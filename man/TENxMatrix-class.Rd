\name{TENxMatrix-class}
\docType{class}

\alias{class:TENxMatrix}
\alias{TENxMatrix-class}
\alias{TENxMatrix}

\alias{DelayedArray,TENxMatrixSeed-method}

% sparsity

\alias{sparsity}
\alias{extractNonzeroDataByCol}

\alias{sparsity,TENxMatrix-method}
\alias{read_sparse_block,TENxMatrix-method}
\alias{extractNonzeroDataByCol,TENxMatrix-method}

\alias{coerce,TENxMatrix,dgCMatrix-method}
\alias{coerce,TENxMatrix,sparseMatrix-method}

\title{10x Genomics datasets as DelayedMatrix objects}

\description{
  The 1.3 Million Brain Cell Dataset and other datasets published by
  10x Genomics use an HDF5-based sparse matrix representation instead
  of the conventional (i.e. dense) HDF5 representation.

  The TENxMatrix class is a \link[DelayedArray]{DelayedMatrix} subclass
  for representing an HDF5-based sparse matrix like one used by
  10x Genomics for the 1.3 Million Brain Cell Dataset.

  All the operations available for \link[DelayedArray]{DelayedMatrix}
  objects work on TENxMatrix objects.
}

\usage{
## Constructor functions:
TENxMatrix(filepath, group="mm10")

## sparsity() and a convenient data extractor:
sparsity(x)
extractNonzeroDataByCol(x, j)
}

\arguments{
  \item{filepath}{
    The path (as a single string) to the HDF5 file where the 10x Genomics
    dataset is located.
  }
  \item{group}{
    The name of the group in the HDF5 file containing the 10x Genomics data.
  }
  \item{x}{
    A TENxMatrix (or \link{TENxMatrixSeed}) object.
  }
  \item{j}{
    An integer vector containing valid column indices.
  }
}

\value{
  \code{TENxMatrix}: A TENxMatrix object.

  \code{sparsity}: The number of zero-valued matrix elements in the
  object divided by its total number of elements (a.k.a. its length).

  \code{extractNonzeroDataByCol}: A \link[IRanges]{NumericList} or
  \link[IRanges]{IntegerList} object \emph{parallel} to \code{j} i.e.
  with one list element per column index in \code{j}. The row indices
  of the values are not returned. Furthermore, the values within a given
  list element can be returned in any order. In particular you should not
  assume that they are ordered by ascending row index.
}

\note{
  If your dataset uses the HDF5-based sparse matrix representation from
  10x Genomics, use the \code{TENxMatrix()} constructor.

  If your dataset uses the conventional (i.e. dense) HDF5 representation,
  use the \code{\link{HDF5Array}()} constructor.
}

\seealso{
  \itemize{
    \item \link{HDF5Array} objects for representing conventional (i.e.
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

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

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

## The structure of this HDF5 file can be seen using the h5ls() command
## from the rhdf5 package:
library(rhdf5)
h5ls(fname)

## The 1.3 Million Brain Cell Dataset is represented by the "mm10"
## group. We point the TENxMatrix() constructor to this group to
## create a TENxMatrix object representing the dataset:
oneM <- TENxMatrix(fname, "mm10")
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
nonzeroes <- extractNonzeroDataByCol(oneM, 1:25000)  # takes < 5 sec.

## The data is returned as an IntegerList object with one list element
## per column and no row indices associated to the values in the object.
## Furthermore, the values within a given list element can be returned
## in any order:
nonzeroes

names(nonzeroes) <- colnames(oneM25k)

## This can be used to compute some simple summaries like the library
## sizes and the number of genes expressed per cell. For these use
## cases, it is a lot more efficient than using colSums(oneM25k) and
## colSums(oneM25k != 0):
lib_sizes3 <- sum(nonzeroes)
n_exprs3 <- lengths(nonzeroes)

## Sanity checks:
stopifnot(all.equal(lib_sizes, lib_sizes3))
stopifnot(all.equal(n_exprs, n_exprs3))
}
\keyword{classes}
\keyword{methods}
