### =========================================================================
### Low-level manipulation of HDF5 Dimension Scale datasets
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5isdimscale()
###

h5isdimscale <- function(filepath, name)
{
    .Call2("C_h5isdimscale", filepath, name, PACKAGE="HDF5Array")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5getdimscales() / h5setdimscales()
###

### Retrieve the names of the existing HDF5 datasets (1 per dimension in
### dataset 'name') currently attached along the dimensions of dataset 'name'
### for Dimension Scale 'scalename'.
h5getdimscales <- function(filepath, name, scalename=NA)
{
    stopifnot(isSingleStringOrNA(scalename))
    scalename <- as.character(scalename)
    .Call2("C_h5getdimscales", filepath, name, scalename,
                               PACKAGE="HDF5Array")
}

### name:      The name of the dataset on which to set Dimension Scales.
### dimscales: A character vector containing the names of the existing HDF5
###            datasets (1 per dimension in dataset 'name') to attach along
###            the dimensions of dataset 'name'. NAs are allowed and, if
###            present, nothing gets attached along the corresponding
###            dimensions.
### scalename: The name of the Dimension Scale (analog to the name of an
###            attribute in R).
h5setdimscales <- function(filepath, name, dimscales, scalename=NA,
                           dry.run=FALSE)
{
    stopifnot(isSingleStringOrNA(scalename),
              is.character(dimscales),
              isTRUEorFALSE(dry.run))
    scalename <- as.character(scalename)
    .Call2("C_h5setdimscales", filepath, name, dimscales, scalename, dry.run,
                               PACKAGE="HDF5Array")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get/set the "dimension labels" of an HDF5 dataset
###
### The "dimension labels" the HDF5 equivalent of the names on 'dimnames(a)'
### in R.
###

h5getdimlabels <- function(filepath, name)
{
    .Call2("C_h5getdimlabels", filepath, name, PACKAGE="HDF5Array")
}

h5setdimlabels <- function(filepath, name, dimlabels)
{
    stopifnot(is.character(dimlabels))
    invisible(.Call2("C_h5setdimlabels", filepath, name, dimlabels,
                                         PACKAGE="HDF5Array"))
}

