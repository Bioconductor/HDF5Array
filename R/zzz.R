.onLoad <- function(libname, pkgname)
{
    init_HDF5_dump_files_global_counter()
    init_HDF5_dump_names_global_counter()
    setHDF5DumpFile()
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

