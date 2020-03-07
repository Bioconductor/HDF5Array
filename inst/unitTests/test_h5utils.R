test_write_read_h5dimnames <- function()
{
    write_h5dimnames <- HDF5Array:::write_h5dimnames
    read_h5dimnames <- HDF5Array:::read_h5dimnames

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
    write_h5dimnames(Adimnames, h5file, Aname)
    current <- read_h5dimnames(h5file, Aname)
    checkIdentical(Adimnames, current)

    ## Write dimnames (with dimlabels) for 'B'.
    Bdimnames <- list(x=letters[1:2], y=NULL, LETTERS[1:4])
    dimscales <- c("X", "Y", "Z")
    write_h5dimnames(Bdimnames, h5file, Bname, dimscales=dimscales)
    current <- read_h5dimnames(h5file, Bname)
    checkIdentical(Bdimnames, current)

    ## Write dimnames for 'C'.
    Cdimnames <- list(NULL, NULL)

    ## Does not actually write anything to the HDF5 file.
    write_h5dimnames(Cdimnames, h5file, Cname)
    current <- read_h5dimnames(h5file, Cname)
    checkIdentical(NULL, current)

    names(Cdimnames) <- c("", "")
    ## Does not actually write anything to the HDF5 file.
    write_h5dimnames(Cdimnames, h5file, Cname)
    current <- read_h5dimnames(h5file, Cname)
    checkIdentical(NULL, current)

    names(Cdimnames)[[1]] <- "x"
    write_h5dimnames(Cdimnames, h5file, Cname, group="more stuff")
    current <- read_h5dimnames(h5file, Cname)
    checkIdentical(Cdimnames, current)
}

