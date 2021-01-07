### =========================================================================
### Summarization of an HDF5 dataset
### -------------------------------------------------------------------------


h5summarize <- function(filepath, name, index=NULL, as.integer=FALSE,
                        op, na.rm=FALSE, verbose=FALSE)
{
    .Call2("C_h5summarize", filepath, name, index, as.integer,
                            op, na.rm, verbose,
                            PACKAGE="HDF5Array")
}

