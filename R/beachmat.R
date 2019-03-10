setupHDF5Matrix <- function(dims, storage.mode)
# This is a convenience wrapper that, if given the matrix dimensions,
# runs a series of R functions and returns the result to the C++ code.
# The idea is to avoid multiple environment calls from C++.
{
    fname <- getHDF5DumpFile(for.use=TRUE)
    dname <- getHDF5DumpName(for.use=TRUE)
    chunk <- getHDF5DumpChunkDim(dims)
    compress <- getHDF5DumpCompressionLevel()
    appendDatasetCreationToHDF5DumpLog(fname, dname, dims, storage.mode, chunk, compress)
    list(fname=fname, dname=dname, chunk=chunk, compress=compress)
}

# Add flags to indicate that, in fact, HDF5Matrix read/write is supported.
beachmat_HDF5Matrix_integer_input <- TRUE
beachmat_HDF5Matrix_logical_input <- TRUE
beachmat_HDF5Matrix_numeric_input <- TRUE
beachmat_HDF5Matrix_character_input <- TRUE
beachmat_HDF5Matrix_integer_output <- TRUE
beachmat_HDF5Matrix_logical_output <- TRUE
beachmat_HDF5Matrix_numeric_output <- TRUE
beachmat_HDF5Matrix_character_output <- TRUE
