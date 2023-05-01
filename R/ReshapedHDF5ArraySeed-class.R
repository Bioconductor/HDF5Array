### =========================================================================
### ReshapedHDF5ArraySeed objects
### -------------------------------------------------------------------------


setClass("ReshapedHDF5ArraySeed",
    contains="HDF5ArraySeed",
    representation(
        reshaped_dim="integer",
        reshaped_chunkdim="integer_OR_NULL"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

### TODO


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() getter
###

### Does NOT access the file.
setMethod("dim", "ReshapedHDF5ArraySeed", function(x) x@reshaped_dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###

.extract_array_from_ReshapedHDF5ArraySeed <- function(x, index)
{
    ans_dim <- S4Arrays:::get_Nindex_lengths(index, dim(x))
    h5mread_from_reshaped(path(x), x@name, x@reshaped_dim, starts=index)
}

setMethod("extract_array", "ReshapedHDF5ArraySeed",
    .extract_array_from_ReshapedHDF5ArraySeed
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### chunkdim() getter
###

### Does NOT access the file.
setMethod("chunkdim", "ReshapedHDF5ArraySeed", function(x) x@reshaped_chunkdim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Return a ReshapedHDF5ArraySeed object with NO dimnames!
### See HDF5ArraySeed() constructor in the HDF5ArraySeed-class.R file for more
### information about the dimnames issue.
ReshapedHDF5ArraySeed <- function(filepath, name, dim, type=NA)
{
    seed <- HDF5ArraySeed(filepath, name, type=type)
    reshaped_dim <- S4Arrays:::normarg_dim(dim)
    collapse_along <- find_dims_to_collapse(reshaped_dim, seed@dim)
    if (is.null(seed@chunkdim)) {
        reshaped_chunkdim <- NULL
    } else {
        reshaped_chunkdim <- collapse_dims(seed@chunkdim, collapse_along)
        reshaped_chunkdim <- as.integer(reshaped_chunkdim)
    }
    new2("ReshapedHDF5ArraySeed", seed,
                                  reshaped_dim=reshaped_dim,
                                  reshaped_chunkdim=reshaped_chunkdim)
}

