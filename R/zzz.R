.onLoad <- function(libname, pkgname)
{
    init_HDF5_dump_files_global_counter()
    init_HDF5_dump_names_global_counter()
    setHDF5DumpDir()
    file.create(get_HDF5_dump_logfile())
    init_HDF5_dataset_creation_global_counter()
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

