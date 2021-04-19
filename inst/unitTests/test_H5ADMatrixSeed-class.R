library(SingleCellExperiment)
library(zellkonverter)

.make_TEST_SCE1 <- function(assays.as=NULL)
{
    cells <- letters[1:10]
    genes <- LETTERS[1:20]
    ncells <- length(cells)
    ngenes <- length(genes)
    u <- matrix(rpois(ngenes * ncells, 5),
                ncol=ncells,
                dimnames=list(genes, cells))
    v <- log2(u + 1)    
    if (!is.null(assays.as)) {
        u <- as(u, assays.as)
        v <- as(v, assays.as)
    }
    pca <- matrix(runif(ncells*5), ncells)
    tsne <- matrix(rnorm(ncells*2), ncells)
    SingleCellExperiment(assays=list(counts=u, logcounts=v),
                         reducedDims=SimpleList(PCA=pca, tSNE=tsne))
}

test_Dense_H5ADMatrixSeed <- function()
{
    ## zellkonverter::writeH5AD() relies on basilisk which isn't supported
    ## on 32-bit Windows.
    if (.Platform$OS.type == "windows" && .Platform$r_arch == "i386")
        return()

    sce <- .make_TEST_SCE1()
    zellkonverter::writeH5AD(sce, "test.h5ad", X_name="counts")
    on.exit(unlink("test.h5ad"))

    seed <- H5ADMatrixSeed("test.h5ad")
    checkTrue(is(seed, "Dense_H5ADMatrixSeed"))
    checkIdentical(dim(sce), dim(seed))
    checkIdentical(dimnames(sce), dimnames(seed))
    checkEquals(assay(sce, "counts"), as.array(seed))
    checkTrue(!is_sparse(seed))

    seed <- H5ADMatrixSeed("test.h5ad", layer="logcounts")
    checkTrue(is(seed, "Dense_H5ADMatrixSeed"))
    checkIdentical(dim(sce), dim(seed))
    checkIdentical(dimnames(sce), dimnames(seed))
    checkEquals(assay(sce, "logcounts"), as.array(seed))
    checkTrue(!is_sparse(seed))
}

test_CSR_H5ADMatrixSeed <- function()
{
    ## zellkonverter::writeH5AD() relies on basilisk which isn't supported
    ## on 32-bit Windows.
    if (.Platform$OS.type == "windows" && .Platform$r_arch == "i386")
        return()

    sce <- .make_TEST_SCE1("dgCMatrix")
    zellkonverter::writeH5AD(sce, "test.h5ad", X_name="counts")
    on.exit(unlink("test.h5ad"))

    seed <- H5ADMatrixSeed("test.h5ad")
    checkTrue(is(seed, "CSR_H5ADMatrixSeed"))
    checkIdentical(dim(sce), dim(seed))
    checkIdentical(dimnames(sce), dimnames(seed))
    checkEquals(as.array(assay(sce, "counts")), as.array(seed))
    checkTrue(is_sparse(seed))

    seed <- H5ADMatrixSeed("test.h5ad", layer="logcounts")
    checkTrue(is(seed, "CSR_H5ADMatrixSeed"))
    checkIdentical(dim(sce), dim(seed))
    checkIdentical(dimnames(sce), dimnames(seed))
    checkEquals(as.array(assay(sce, "logcounts")), as.array(seed))
    checkTrue(is_sparse(seed))
}

