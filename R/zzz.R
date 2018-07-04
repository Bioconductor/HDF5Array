.onLoad <- function(libname, pkgname)
{
    init_HDF5_dump_files_global_counter()
    init_HDF5_dump_names_global_counter()
    setHDF5DumpDir()
    setHDF5DumpChunkMaxLength()
    setHDF5DumpChunkShape()
    setHDF5DumpCompressionLevel()
    file.create(get_HDF5_dump_logfile())
    init_HDF5_dataset_creation_global_counter()
}

.test <- function()
{
    BiocGenerics:::testPackage("HDF5Array")

    ## Skip this on 32-bit Windows to avoid 'R CMD check' TIMEOUT on the
    ## Windows build machines.
    if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
        setRealizationBackend("HDF5Array")
        BiocGenerics:::testPackage("DelayedArray")
    }
}

