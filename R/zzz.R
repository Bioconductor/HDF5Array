.onLoad <- function(libname, pkgname)
{
    setHDF5DumpFile()
    setHDF5DumpName()
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

