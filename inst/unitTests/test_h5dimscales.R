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

    dsnames <- c("B", "stuff/C", NA, "more stuff/D")
    target <- c(TRUE, TRUE, FALSE, TRUE)
    current <- h5setdimscales(h5file, "A", dsnames, dry.run=TRUE)
    checkIdentical(target, current)
    current <- h5setdimscales(h5file, "A", dsnames)
    checkIdentical(target, current)
    checkIdentical(logical(4), h5setdimscales(h5file, "A", dsnames))

    target <- paste0("/", dsnames); target[is.na(dsnames)] <- NA
    current <- h5getdimscales(h5file, "A")
    checkIdentical(target, current)

    target <- rep(NA_character_, 4)
    current <- h5getdimscales(h5file, "A", "foo")
    checkIdentical(target, current)

    dsnames <- c(NA, "E", "more stuff/E", "E")
    target <- c(FALSE, TRUE, TRUE, TRUE)
    current <- h5setdimscales(h5file, "A", dsnames, "foo", dry.run=TRUE)
    checkIdentical(target, current)
    current <- h5setdimscales(h5file, "A", dsnames, "foo")
    checkIdentical(target, current)
    checkIdentical(logical(4), h5setdimscales(h5file, "A", dsnames, "foo"))

    target <- paste0("/", dsnames); target[is.na(dsnames)] <- NA
    current <- h5getdimscales(h5file, "A", "foo")
    checkIdentical(target, current)

    for (name in paste0("Adim", 1:4))
        writeHDF5Array(matrix(1:2), h5file, name)
    for (scalename in c(NA, "foo", "bar"))
        scales <- h5getdimscales(h5file, "A", scalename)
        for (dsname1 in c(NA, "Adim1", "B"))
          for (dsname2 in c(NA, "Adim2", "E"))
            for (dsname3 in c(NA, "Adim3", "bogus"))
              for (dsname4 in c(NA, "Adim4", "bogus")) {
                dsnames <- c(dsname1, dsname2, dsname3, dsname4)
                scales2 <- paste0("/", dsnames)
                ok <- scales2 == scales | is.na(scales2) | is.na(scales)
                if (dsname1 %in% "B" && !(scalename %in% NA))
                    ok[[1]] <- FALSE
                if (dsname2 %in% "E" && !(scalename %in% "foo"))
                    ok[[2]] <- FALSE
                if (!all(ok) || "bogus" %in% dsnames) {
                    checkException(h5setdimscales(h5file, "A", dsnames,
                                                  scalename),
                                   silent=TRUE)
                    next
                }
                target <- !is.na(dsnames)
                current <- h5setdimscales(h5file, "A", dsnames, scalename,
                                          dry.run=TRUE)
                checkIdentical(target, current)
              }
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

