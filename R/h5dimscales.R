### =========================================================================
### Low-level manipulation of HDF5 Dimension Scale datasets
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5setdimscales() / h5getdimscales()
###

### name:      The name of the dataset on which to set Dimension Scales.
### dsnames:   A character vector containing the names of the existing HDF5
###            datasets (1 per dimension in dataset 'name') to attach along
###            the dimensions of dataset 'name'. NAs are allowed and, if
###            present, nothing will get attached along the corresponding
###            dimensions.
### scalename: The name of the Dimension Scale (analog to the name of an
###            attribute in R).
h5setdimscales <- function(filepath, name, dsnames, scalename=NA_character_,
                           dry.run=FALSE)
{
    stopifnot(isSingleStringOrNA(scalename),
              is.character(dsnames),
              isTRUEorFALSE(dry.run))
    scalename <- as.character(scalename)
    .Call2("C_h5setdimscales", filepath, name, dsnames, scalename, dry.run,
                               PACKAGE="HDF5Array")
}

### Retrieve the names of the existing HDF5 datasets (1 per dimension in
### dataset 'name') currently attached along the dimensions of dataset 'name'
### for Dimension Scale 'scalename'.
h5getdimscales <- function(filepath, name, scalename=NA_character_)
{
    stopifnot(isSingleStringOrNA(scalename))
    scalename <- as.character(scalename)
    .Call2("C_h5getdimscales", filepath, name, scalename,
                               PACKAGE="HDF5Array")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set/get the "dimension labels" of an HDF5 dataset
###
### The "dimension labels" the HDF5 equivalent of the names on 'dimnames(a)'
### in R.
###

h5setdimlabels <- function(filepath, name, labels)
{
    stopifnot(is.character(labels))
    invisible(.Call2("C_h5setdimlabels", filepath, name, labels,
                                         PACKAGE="HDF5Array"))
}

h5getdimlabels <- function(filepath, name)
{
    .Call2("C_h5getdimlabels", filepath, name, PACKAGE="HDF5Array")
}

