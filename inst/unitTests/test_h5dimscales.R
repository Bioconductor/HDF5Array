test_h5getdimscales_h5setdimscales <- function()
{
    h5getdimscales <- HDF5Array:::h5getdimscales
    h5setdimscales <- HDF5Array:::h5setdimscales

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

    dimscales <- c("B", "stuff/C", NA, "more stuff/D")
    target <- c(TRUE, TRUE, FALSE, TRUE)
    current <- h5setdimscales(h5file, "A", dimscales, dry.run=TRUE)
    checkIdentical(target, current)
    current <- h5setdimscales(h5file, "A", dimscales)
    checkIdentical(target, current)
    checkIdentical(logical(4), h5setdimscales(h5file, "A", dimscales))

    target <- paste0("/", dimscales); target[is.na(dimscales)] <- NA
    current <- h5getdimscales(h5file, "A")
    checkIdentical(target, current)

    target <- rep(NA_character_, 4)
    current <- h5getdimscales(h5file, "A", scalename="foo")
    checkIdentical(target, current)

    dimscales <- c(NA, "E", "more stuff/E", "E")
    target <- c(FALSE, TRUE, TRUE, TRUE)
    current <- h5setdimscales(h5file, "A", dimscales, scalename="foo",
                              dry.run=TRUE)
    checkIdentical(target, current)
    current <- h5setdimscales(h5file, "A", dimscales, scalename="foo")
    checkIdentical(target, current)
    checkIdentical(logical(4), h5setdimscales(h5file, "A", dimscales,
                                              scalename="foo"))

    target <- paste0("/", dimscales); target[is.na(dimscales)] <- NA
    current <- h5getdimscales(h5file, "A", scalename="foo")
    checkIdentical(target, current)

    for (name in paste0("Adim", 1:4))
        writeHDF5Array(matrix(1:2), h5file, name)
    for (scalename in c(NA, "foo", "bar"))
        dimscales0 <- h5getdimscales(h5file, "A", scalename=scalename)
        for (scale1 in c(NA, "Adim1", "B"))
          for (scale2 in c(NA, "Adim2", "E"))
            for (scale3 in c(NA, "Adim3", "bogus"))
              for (scale4 in c(NA, "Adim4", "bogus")) {
                dimscales <- c(scale1, scale2, scale3, scale4)
                dimscales2 <- paste0("/", dimscales)
                ok <- is.na(dimscales) |
                      is.na(dimscales0) |
                      dimscales2 == dimscales0
                if (scale1 %in% "B" && !(scalename %in% NA))
                    ok[[1]] <- FALSE
                if (scale2 %in% "E" && !(scalename %in% "foo"))
                    ok[[2]] <- FALSE
                if (!all(ok) || "bogus" %in% dimscales) {
                    checkException(h5setdimscales(h5file, "A", dimscales,
                                                  scalename=scalename),
                                   silent=TRUE)
                    next
                }
                target <- !is.na(dimscales)
                current <- h5setdimscales(h5file, "A", dimscales,
                                          scalename=scalename, dry.run=TRUE)
                checkIdentical(target, current)
              }
}

test_h5getdimlabels_h5setdimlabels <- function()
{
    h5getdimlabels <- HDF5Array:::h5getdimlabels
    h5setdimlabels <- HDF5Array:::h5setdimlabels

    h5file <- tempfile(fileext=".h5")
    h5createFile(h5file)
    h5createGroup(h5file, "stuff")
    h5createGroup(h5file, "more stuff")
    writeHDF5Array(array(1:24, 4:1), h5file, "stuff/A")

    checkIdentical(NULL, h5getdimlabels(h5file, "stuff/A"))

    dimlabels <- c("dim1", NA, NA, "dim4")
    h5setdimlabels(h5file, "stuff/A", dimlabels)
    target <- c("dim1", "", "", "dim4")
    current <- h5getdimlabels(h5file, "stuff/A")
    checkIdentical(target, current)

    dimlabels <- c(NA, "dim2", NA, "XXXX")
    h5setdimlabels(h5file, "stuff/A", dimlabels)
    target <- c("dim1", "dim2", "", "XXXX")
    current <- h5getdimlabels(h5file, "stuff/A")
    checkIdentical(target, current)

    dimlabels <- c("", "", NA, "")
    h5setdimlabels(h5file, "stuff/A", dimlabels)
    checkIdentical(NULL, h5getdimlabels(h5file, "stuff/A"))
}

