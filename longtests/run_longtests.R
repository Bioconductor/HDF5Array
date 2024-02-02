require("HDF5Array") || stop("unable to load HDF5Array package")
setAutoRealizationBackend("HDF5Array")
## The tests in DelayedArray require the DelayedMatrixStats and genefilter
## packages so HDF5Array must have them in Suggests.
BiocGenerics:::testPackage("DelayedArray")

