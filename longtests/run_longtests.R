require("HDF5Array") || stop("unable to load HDF5Array package")
setRealizationBackend("HDF5Array")
BiocGenerics:::testPackage("DelayedArray")

