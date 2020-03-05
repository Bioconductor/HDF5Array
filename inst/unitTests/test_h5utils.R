test_h5writeDimnames_h5readDimnames <- function()
{
    h5writeDimnames <- HDF5Array:::h5writeDimnames
    h5readDimnames <- HDF5Array:::h5readDimnames

    h5file <- tempfile(fileext=".h5")
    h5createFile(h5file)
    h5createGroup(h5file, "stuff")
    Aname <- "stuff/A"
    A <- writeHDF5Array(array(1:24, 4:1), h5file, Aname)
    Bname <- "stuff/B"
    B <- writeHDF5Array(array(101:124, 2:4), h5file, Bname)
    Cname <- "stuff/C"
    C <- writeHDF5Array(matrix(201:206, ncol=2), h5file, Cname)

    ## Write dimnames for 'A'.
    Adimnames <- list(letters[1:4], NULL, 11:12, NULL)
    Adimnames_group <- paste(dirname(Aname),
                             paste0(basename(Aname), "dimnames"),
                             sep="/")
    Adsnames <- paste(Adimnames_group, sprintf("%03d", seq_along(dim(A))),
                      sep="/")
    h5createGroup(h5file, Adimnames_group)
    h5writeDimnames(h5file, Aname, Adimnames, Adsnames)
    target <- lapply(Adimnames,
        function(dn) if (is.null(dn)) NULL else as.character(dn)
    )
    current <- h5readDimnames(h5file, Aname)
    checkIdentical(target, current)

    ## Write dimnames (with dimlabels) for 'B'.
    Bdimnames <- list(X=letters[1:2], Y=NULL, LETTERS[1:4])
    h5writeDimnames(h5file, Bname, Bdimnames, c("X", "Y", "Z"))
    current <- h5readDimnames(h5file, Bname)
    checkIdentical(Bdimnames, current)

    ## Write dimnames for 'C'.
    Cdimnames <- list(NULL, NULL)

    ## Does not actually write anything to the HDF5 file.
    h5writeDimnames(h5file, Cname, Cdimnames, c("C1", "C2"))
    current <- h5readDimnames(h5file, Cname)
    checkIdentical(NULL, current)

    names(Cdimnames) <- c("", "")
    ## Does not actually write anything to the HDF5 file.
    h5writeDimnames(h5file, Cname, Cdimnames, c("C1", "C2"))
    current <- h5readDimnames(h5file, Cname)
    checkIdentical(NULL, current)

    names(Cdimnames)[[1]] <- "x"
    h5writeDimnames(h5file, Cname, Cdimnames, c("C1", "C2"))
    current <- h5readDimnames(h5file, Cname)
    checkIdentical(Cdimnames, current)
}

