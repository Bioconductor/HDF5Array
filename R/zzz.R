.onLoad <- function(libname, pkgname)
{
    setHDF5DumpFile()
    setHDF5DumpName()
    options(HDF5Array.block.size=DEFAULT_BLOCK_SIZE)
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

