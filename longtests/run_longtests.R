require("HDF5Array") || stop("unable to load HDF5Array package")
setAutoRealizationBackend("HDF5Array")
BiocGenerics:::testPackage("DelayedArray")

