.onLoad <- function(libname, pkgname)
{
    setHDF5OutputFile()
    setHDF5OutputName()
    options(HDF5Array.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

