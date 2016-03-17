.onLoad <- function(libname, pkgname)
{
    setHDF5ArrayOutputFile()
    options(HDF5Array.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

