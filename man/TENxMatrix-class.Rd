\name{TENxMatrix-class}
\docType{class}

% TENxMatrixSeed

\alias{class:TENxMatrixSeed}
\alias{TENxMatrixSeed-class}
\alias{TENxMatrixSeed}

\alias{path,TENxMatrixSeed-method}
\alias{path<-,TENxMatrixSeed-method}
\alias{dim,TENxMatrixSeed-method}
\alias{dimnames,TENxMatrixSeed-method}
\alias{extract_array,TENxMatrixSeed-method}
\alias{chunkdim,TENxMatrixSeed-method}

\alias{show,TENxMatrixSeed-method}

% TENxMatrix

\alias{class:TENxMatrix}
\alias{TENxMatrix-class}
\alias{TENxMatrix}

\alias{DelayedArray,TENxMatrixSeed-method}

% sparsity() and extractNonZeroDataByCol()

\alias{sparsity}
\alias{sparsity,TENxMatrixSeed-method}
\alias{sparsity,TENxMatrix-method}

\alias{extractNonZeroDataByCol}
\alias{extractNonZeroDataByCol,TENxMatrixSeed-method}
\alias{extractNonZeroDataByCol,TENxMatrix-method}

\title{10x Genomics datasets as DelayedArray objects}

\description{
  The 1.3 Million Brain Cell Dataset and other datasets published by
  10x Genomics use an HDF5-based sparse matrix representation instead
  of the conventional (i.e. dense) HDF5 representation.

  We provide 2 classes for representing a 10x Genomics dataset as a
  matrix-like object in R:
  \itemize{
    \item TENxMatrix: A high-level class that extends
          \link[DelayedArray]{DelayedArray}.
          All the operations available for \link[DelayedArray]{DelayedArray}
          objects work on TENxMatrix objects.

    \item TENxMatrixSeed: A low-level class for pointing to a 10x Genomics
          dataset.
          No operation can be performed directly on a TENxMatrixSeed object.
          It first needs to be wrapped in a \link[DelayedArray]{DelayedArray}
          object.
          A TENxMatrix object is just a TENxMatrixSeed object wrapped in
          a \link[DelayedArray]{DelayedArray} object.
  }
}

\usage{
## Constructor functions:
TENxMatrix(filepath, group="mm10")
TENxMatrixSeed(filepath, group="mm10")

## sparsity() and a convenient data extractor:
sparsity(x)
extractNonZeroDataByCol(x, j)
}

\arguments{
  \item{filepath}{
    The path (as a single character string) to the HDF5 file where the
    10x Genomics dataset is located.
  }
  \item{group}{
    The name of the group in the HDF5 file containing the 10x Genomics data.
  }
  \item{x}{
    A TENxMatrix or TENxMatrixSeed object.
  }
  \item{j}{
    An integer vector containing valid column indices.
  }
}

\value{
  \code{TENxMatrix}: A TENxMatrix object.

  \code{TENxMatrixSeed}: A TENxMatrixSeed object.

  \code{sparsity}: The number of zero-valued matrix elements in the
  object divided by its total number of elements (a.k.a. its length).

  \code{extractNonZeroDataByCol}: A \link[IRanges]{NumericList} or
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

    \item \link[DelayedArray]{DelayedArray} objects in the \pkg{DelayedArray}
          package.

    \item \code{\link{writeTENxMatrix}} for writing a matrix-like object
          as an HDF5-based sparse matrix.

    \item \code{\link[parallel]{detectCores}} from the \pkg{parallel}
          package.

    \item \code{\link[BiocParallel]{MulticoreParam}} and
          \code{\link[BiocParallel]{register}} from the \pkg{BiocParallel}
          package.

    \item \code{\link[DelayedArray]{setDefaultBlockSize}} in the
          \pkg{DelayedArray} package.

    \item \code{\link[DelayedArray]{colGrid}} and
          \code{\link[DelayedArray]{blockApply}} in the
          \pkg{DelayedArray} package.

    \item \code{\link[rhdf5]{h5ls}} in the \pkg{rhdf5} package.

    \item The \pkg{\link[rhdf5]{rhdf5}} package on top of which TENxMatrix
          and TENxMatrixSeed objects are implemented.

    \item \link[IRanges]{NumericList} and \link[IRanges]{IntegerList}
          objects in the \pkg{IRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## THE "1.3 Million Brain Cell Dataset" AS A DelayedArray OBJECT
## ---------------------------------------------------------------------
## The 1.3 Million Brain Cell Dataset from 10x Genomics is available via
## ExperimentHub:
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "TENxBrainData")
fname <- hub[["EH1039"]]

## The structure of the file can be seen using the h5ls() command from
## the rhdf5 package:
library(rhdf5)
h5ls(fname)

## The 1.3 Million Brain Cell Dataset is represented by the "mm10"
## group. We point the TENxMatrix() constructor to this group to
## create a TENxMatrix object representing the dataset:
oneM <- TENxMatrix(fname, "mm10")
oneM
sparsity(oneM)

is(oneM, "DelayedArray")  # TRUE

## Some examples of delayed operations:
oneM != 0
oneM^2

## ---------------------------------------------------------------------
## SOME EXAMPLES OF ROW/COL SUMMARIZATION
## ---------------------------------------------------------------------
## In order to reduce computation times, we'll use only the first
## 50000 columns of the 1.3 Million Brain Cell Dataset:
oneM50k <- oneM[ , 1:50000]

## Row/col summarization methods like rowSums() use a block-processing
## mechanism behind the scene that can be controlled via global
## settings. 2 important settings that can have a strong impact on
## performance are the number of workers and the default block size.
## On a modern Linux laptop with 8 core (as reported by
## parallel::detectCores()) and 16 Gb of RAM, reasonably good
## performance is achieved by setting the number of workers to 6
## and the default block size to 500 Mb:
register(MulticoreParam(workers=6))
setDefaultBlockSize(5e8)
DelayedArray:::set_verbose_block_processing(TRUE)

## We're ready to compute the library sizes, number of genes expressed
## per cell, and average expression across cells:
system.time(lib_sizes <- colSums(oneM50k))
system.time(n_exprs <- colSums(oneM50k != 0))
system.time(ave_exprs <- rowMeans(oneM50k))

## Note that the 3 computations above load the data in oneM50k 3 times
## in memory. This can be avoided by computing the 3 summarizations in
## a single pass with blockApply(). First we define the function that
## we're going to apply to each block of data:
FUN <- function(block)
  list(colSums(block), colSums(block != 0), rowSums(block))

## Then we call blockApply() to apply FUN() to each block. The blocks
## are defined by the grid passed to the 'grid' argument. In this case
## we supply a grid made with colGrid() to generate blocks of full
## columns (see ?colGrid for more information):
system.time({
  block_results <- blockApply(oneM50k, FUN, grid=colGrid(oneM50k))
})

## 'block_results' is a list with 1 list element per block in
## colGrid(oneM50k). Each list element is the result that was obtained
## by applying FUN() on the block so is itself a list of length 3.
## Let's combine the results:
lib_sizes2 <- unlist(lapply(block_results, `[[`, 1L))
n_exprs2 <- unlist(lapply(block_results, `[[`, 2L))
block_rowsums <- unlist(lapply(block_results, `[[`, 3L), use.names=FALSE)
tot_exprs <- rowSums(matrix(block_rowsums, nrow=nrow(oneM50k)))
ave_exprs2 <- setNames(tot_exprs / ncol(oneM50k), rownames(oneM50k))

## Sanity checks:
stopifnot(all.equal(lib_sizes, lib_sizes2))
stopifnot(all.equal(n_exprs, n_exprs2))
stopifnot(all.equal(ave_exprs, ave_exprs2))

## Reset default block size to factory settings:
setDefaultBlockSize()
DelayedArray:::set_verbose_block_processing(FALSE)

## ---------------------------------------------------------------------
## extractNonZeroDataByCol()
## ---------------------------------------------------------------------
## extractNonZeroDataByCol() provides a convenient and very efficient
## way to extract the non-zero data in a compact form:
non_zeroes <- extractNonZeroDataByCol(oneM, 1:50000)  # takes < 5 sec.

## The data is returned as an IntegerList object with one list element
## per column and no row indices associated to the values in the object.
## Furthermore, the values within a given list element can be returned
## in any order:
non_zeroes

names(non_zeroes) <- colnames(oneM50k)

## This can be used to compute some simple summaries like the library
## sizes and the number of genes expressed per cell. For these use
## cases, it is a lot more efficient than using colSums(oneM50k) and
## colSums(oneM50k != 0):
lib_sizes3 <- sum(non_zeroes)
n_exprs3 <- lengths(non_zeroes)

## Sanity checks:
stopifnot(all.equal(lib_sizes, lib_sizes3))
stopifnot(all.equal(n_exprs, n_exprs3))
}
\keyword{methods}
\keyword{classes}
