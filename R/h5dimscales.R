### =========================================================================
### Low-level manipulation of HDF5 Dimension Scale datasets
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### h5setdimscales() / h5getdimscales()
###

### 'name':      The name of the dataset on which to set dimension scale
###              datasets.
### 'scalename': The name of the dimension scale (analog to the name of an
###              attribute in R).
### 'dsnames':   The names of the existing HDF5 datasets (1 per dimension in
###              dataset 'name') to attach along the dimensions of dataset
###              'name'. An NA means that nothing is attached along the
###              corresponding dimension.
h5setdimscales <- function(filepath, name, scalename, dsnames)
{
    stopifnot(is.null(scalename) || isSingleString(scalename),
              is.character(dsnames))
    .Call2("C_h5setdimscales", filepath, name, scalename, dsnames,
                               PACKAGE="HDF5Array")
}

### Retrieve the names of the existing HDF5 datasets (1 per dimension in
### dataset 'name') currently attached along the dimensions of dataset 'name'
### for dimension scale 'scalename'.
h5getdimscales <- function(filepath, name, scalename)
{
    stopifnot(is.null(scalename) || isSingleString(scalename))
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
    .Call2("C_h5setdimlabels", filepath, name, labels, PACKAGE="HDF5Array")
}

h5getdimlabels <- function(filepath, name)
{
    .Call2("C_h5getdimlabels", filepath, name, PACKAGE="HDF5Array")
}

