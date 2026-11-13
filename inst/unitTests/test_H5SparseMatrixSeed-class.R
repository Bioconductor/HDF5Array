.write_h5_sparse_matrix <-
    function(filepath, group, shape, data, indices, indptr, subdata=NULL)
{
    rhdf5::h5createGroup(filepath, group)
    rhdf5::h5write(shape, filepath, paste0(group, "/shape"))
    name <- "/data"
    if (!is.null(subdata)) {
        rhdf5::h5createGroup(filepath, paste0(group, name))
        name <- paste0(name, "/", subdata)
    }
    rhdf5::h5write(data, filepath, paste0(group, name))
    rhdf5::h5write(indices, filepath, paste0(group, "/indices"))
    rhdf5::h5write(indptr, filepath, paste0(group, "/indptr"))
}

test_H5SparseMatrixSeed <- function()
{
    set.seed(2009)
    svt <- poissonSparseMatrix(12, 20, density=0.4) * c(0.1, -0.1)
    tsvt <- t(svt)
    m <- as.matrix(svt)
    dgc <- as(svt, "dgCMatrix")

    filepath <- tempfile(fileext=".h5")
    rhdf5::h5createFile(filepath)
    shape <- dim(dgc)
    data <- dgc@x
    indices <- dgc@i
    indptr <- dgc@p

    ## --- Create the H5SparseMatrixSeed objects ---

    ## Store nonzero data in "<group>/data" dataset (the default):
    .write_h5_sparse_matrix(filepath, "SM1", shape, data, indices, indptr)
    seed1 <- H5SparseMatrixSeed(filepath, "SM1")
    checkIdentical(group(seed1), "/SM1")

    ## Store nonzero data in "<group>/data/counts" dataset:
    .write_h5_sparse_matrix(filepath, "SM2", shape, data, indices, indptr,
                            subdata="counts")
    checkException(H5SparseMatrixSeed(filepath, "SM2"))
    seed2 <- H5SparseMatrixSeed(filepath, "SM2", subdata="counts")
    checkIdentical(group(seed2), "/SM2")

    ## Store nonzero data in "<group>/data/flags" dataset.
    ## WARNING: This dataset will actually contain zeros. Is this allowed?
    ## This kind of breaks 'nzcount(seed3)'. More generally, how will it
    ## impact the proper operations of H5SparseMatrixSeed 'seed3'?
    rhdf5::h5write(as.raw(data > 0), filepath, "SM2/data/flags")
    seed3 <- H5SparseMatrixSeed(filepath, "SM2", subdata="flags")
    checkIdentical(group(seed3), "/SM2")

    ## --- Test the H5SparseMatrixSeed objects ---

    index_list <- list(list(NULL, NULL),
                       list(5:8, NULL),
                       list(NULL, 2:15),
                       list(2:9, 4:18),
                       list(integer(0), 4:18),
                       list(12L, integer(0)),
                       list(c(11:12, 1L), c(4:18, 20:15)))

    ## Leaving 'seed3' out of this for now.
    #for (seed in list(seed1, seed2, seed3)) {
    for (seed in list(seed1, seed2)) {
        checkTrue(is(seed, "CSC_H5SparseMatrixSeed"))
        checkIdentical(path(seed), filepath)
        checkIdentical(dim(seed), dim(svt))
        checkIdentical(chunkdim(seed), c(nrow(svt), 1L))
        checkIdentical(type(seed), type(svt))
        checkTrue(is_sparse(seed))
        checkIdentical(nzcount(seed), nzcount(svt))
        for (index in index_list) {
            checkIdentical(extract_sparse_array(seed, index),
                           extract_sparse_array(svt, index))
            checkIdentical(extract_array(seed, index),
                           extract_array(svt, index))
        }
        checkIdentical(as.matrix(seed), as.matrix(svt))

        tseed <- t(seed)
        checkTrue(is(tseed, "CSR_H5SparseMatrixSeed"))
        checkIdentical(path(tseed), filepath)
        checkIdentical(group(tseed), group(seed))
        checkIdentical(dim(tseed), rev(dim(svt)))
        checkIdentical(chunkdim(tseed), c(1L, nrow(svt)))
        checkIdentical(type(tseed), type(svt))
        checkTrue(is_sparse(tseed))
        checkIdentical(nzcount(tseed), nzcount(svt))
        for (index in index_list) {
            tindex <- rev(index)
            checkIdentical(extract_sparse_array(tseed, tindex),
                           extract_sparse_array(tsvt, tindex))
            checkIdentical(extract_array(tseed, tindex),
                           extract_array(tsvt, tindex))
        }
        checkIdentical(as.matrix(tseed), as.matrix(tsvt))
    }
}

