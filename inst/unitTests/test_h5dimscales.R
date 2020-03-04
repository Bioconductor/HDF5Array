test_h5setdimscales_h5getdimscales <- function()
{
    h5setdimscales <- HDF5Array:::h5setdimscales
    h5getdimscales <- HDF5Array:::h5getdimscales

    h5file <- tempfile(fileext=".h5")
    h5createFile(h5file)
    h5createGroup(h5file, "stuff")
    h5createGroup(h5file, "more stuff")
    writeHDF5Array(array(1:24, 4:1), h5file, "A")
    writeHDF5Array(matrix(11:16, ncol=3), h5file, "B")
    writeHDF5Array(matrix(letters, ncol=2), h5file, "stuff/C")
    writeHDF5Array(matrix(runif(10), ncol=5), h5file, "more stuff/D")
    writeHDF5Array(matrix(1:4, ncol=2), h5file, "more stuff/E")
    writeHDF5Array(matrix(1:8, ncol=2), h5file, "E")

    target <- rep(NA_character_, 4)
    current <- h5getdimscales(h5file, "A")
    checkIdentical(target, current)

    h5setdimscales(h5file, "A", c("B", "stuff/C", NA, "more stuff/D"))
    target <- c("/B", "/stuff/C", NA, "/more stuff/D")
    current <- h5getdimscales(h5file, "A")
    checkIdentical(target, current)

    target <- rep(NA_character_, 4)
    current <- h5getdimscales(h5file, "A", "foo")
    checkIdentical(target, current)

    h5setdimscales(h5file, "A", c(NA, "E", "more stuff/E", "E"), "foo")
    target <- c(NA, "/E", "/more stuff/E", "/E")
    current <- h5getdimscales(h5file, "A", "foo")
    checkIdentical(target, current)
}

test_h5setdimlabels_h5getdimlabels <- function()
{
    h5setdimlabels <- HDF5Array:::h5setdimlabels
    h5getdimlabels <- HDF5Array:::h5getdimlabels

    h5file <- tempfile(fileext=".h5")
    h5createFile(h5file)
    h5createGroup(h5file, "stuff")
    h5createGroup(h5file, "more stuff")
    writeHDF5Array(array(1:24, 4:1), h5file, "stuff/A")

    checkIdentical(NULL, h5getdimlabels(h5file, "stuff/A"))

    labels <- c("dim1", NA, NA, "dim4")
    h5setdimlabels(h5file, "stuff/A", labels)
    target <- c("dim1", "", "", "dim4")
    current <- h5getdimlabels(h5file, "stuff/A")
    checkIdentical(target, current)

    labels <- c(NA, "dim2", NA, "XXXX")
    h5setdimlabels(h5file, "stuff/A", labels)
    target <- c("dim1", "dim2", "", "XXXX")
    current <- h5getdimlabels(h5file, "stuff/A")
    checkIdentical(target, current)

    labels <- c("", "", NA, "")
    h5setdimlabels(h5file, "stuff/A", labels)
    checkIdentical(NULL, h5getdimlabels(h5file, "stuff/A"))
}

