.onLoad <- function(libname, pkgname)
{
    setHDF5ArrayOutputFile()
    setHDF5ArrayOutputName()
    options(HDF5Array.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

