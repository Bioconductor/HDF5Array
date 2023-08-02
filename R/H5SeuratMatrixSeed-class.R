### =========================================================================
### H5SeuratMatrixSeed objects
### -------------------------------------------------------------------------

setClass("H5SeuratMatrixSeed", contains="CSC_H5SparseMatrixSeed")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

.read_h5seurat_component <- function(filepath, name){
    if (h5exists(filepath, name))
        return(as.vector(h5mread(filepath, name)))
    NULL
}

.load_tenx_rownames <- function(filepath, group)
{
    name <- .find_rownames_dataset(filepath, group)
    # For h5seurat data
    # Find /path/to/group/../
    if (is.null(name))
        name <- .find_rownames_dataset(filepath, gsub('/[^/]+$','',group))
    if (is.null(name))
        return(NULL)
    read_h5sparse_component(filepath, group, name)
}

#.find_rownames_dataset <- function(filepath, group)
#{
#    name <- "genes"
#    if (h5exists(filepath, paste(group, name, sep="/")))
#        return(name)
#    name <- "features/id"
#    if (h5exists(filepath, paste(group, name, sep="/")))
#        return(name)
#    NULL
#}

### Return the rownames of the matrix.
.load_tenx_rownames <- function(filepath, group)
{
    name <- .find_rownames_dataset(filepath, group)
    if (is.null(name)){
        # Support h5seurat
        h5SeuratName <- file.path(gsub('/[^/]+$','',group), 'features')
        return(.read_h5seurat_component(filepath, h5SeuratName))}
    read_h5sparse_component(filepath, group, name)
}

### Return the colnames of the matrix.
.load_tenx_barcodes <- function(filepath, group)
{
    if (!h5exists(filepath, paste0(group, "/barcodes"))){
        return(.read_h5seurat_component(filepath, 'cell.names'))}
    read_h5sparse_component(filepath, group, "barcodes")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

H5SeuratMatrixSeed <- function(filepath, group="matrix")
{
    seed0 <- H5SparseMatrixSeed(filepath, group)

    ## dimnames
    rownames <- .load_tenx_rownames(seed0@filepath, seed0@group)
    stopifnot(is.null(rownames) || length(rownames) == seed0@dim[[1L]])
    colnames <- .load_tenx_barcodes(seed0@filepath, seed0@group)
    stopifnot(is.null(colnames) || length(colnames) == seed0@dim[[2L]])
    dimnames <- list(rownames, colnames)

    new2("H5SeuratMatrixSeed", seed0, dimnames=dimnames)
}

