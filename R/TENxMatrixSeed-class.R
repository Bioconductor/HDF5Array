### =========================================================================
### TENxMatrixSeed objects
### -------------------------------------------------------------------------


setClass("TENxMatrixSeed", contains="H5SparseMatrixSeed")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

### Return the rownames of the matrix.
.get_tenx_genes <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/genes")))
        return(NULL)
    read_h5sparse_component(filepath, group, "genes")
}

### Currently unused.
.get_tenx_gene_names <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/gene_names")))
        return(NULL)
    read_h5sparse_component(filepath, group, "gene_names")
}

### Return the colnames of the matrix.
.get_tenx_barcodes <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/barcodes")))
        return(NULL)
    read_h5sparse_component(filepath, group, "barcodes")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

TENxMatrixSeed <- function(filepath, group="mm10")
{
    seed0 <- H5SparseMatrixSeed(filepath, group)

    ## dimnames
    rownames <- .get_tenx_genes(seed0@filepath, seed0@group)
    stopifnot(is.null(rownames) || length(rownames) == seed0@dim[[1L]])
    colnames <- .get_tenx_barcodes(seed0@filepath, seed0@group)
    stopifnot(is.null(colnames) || length(colnames) == seed0@dim[[2L]])
    dimnames <- list(rownames, colnames)

    new2("TENxMatrixSeed", seed0, dimnames=dimnames)
}

