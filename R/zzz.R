### Placeholders, initialized in .onLoad()
.DUMP_FILES_GLOBAL_COUNTER_FILEPATH <- NULL
.DUMP_NAMES_GLOBAL_COUNTER_FILEPATH <- NULL
.DATASET_CREATION_GLOBAL_COUNTER_FILEPATH <- NULL

.onLoad <- function(libname, pkgname)
{
    .DUMP_FILES_GLOBAL_COUNTER_FILEPATH <<-
        init_HDF5_dump_files_global_counter()
    .DUMP_NAMES_GLOBAL_COUNTER_FILEPATH <<-
        init_HDF5_dump_names_global_counter()
    setHDF5DumpDir()
    setHDF5DumpChunkLength()
    setHDF5DumpChunkShape()
    setHDF5DumpCompressionLevel()
    file.create(get_HDF5_dump_logfile())
    .DATASET_CREATION_GLOBAL_COUNTER_FILEPATH <<-
        init_HDF5_dataset_creation_global_counter()
}

.onUnload <- function(libpath)
{
    file.remove(.DUMP_FILES_GLOBAL_COUNTER_FILEPATH)
    file.remove(.DUMP_NAMES_GLOBAL_COUNTER_FILEPATH)
    file.remove(.DATASET_CREATION_GLOBAL_COUNTER_FILEPATH)
}

.test <- function() BiocGenerics:::testPackage("HDF5Array")

