### =========================================================================
### DelayedMatrix objects
### -------------------------------------------------------------------------


### Extending DataTable gives us a few things for free (head(), tail(),
### etc...)
setClass("DelayedMatrix",
    contains=c("DelayedArray", "DataTable"),
    representation(
        ## x@N1 and x@N2 must be 2 integers such that
        ##     1 <= x@N1 < x@N2 <= length(x@index)
        N1="integer",  # single integer
        N2="integer"   # single integer
    ),
    prototype(
        seeds=list(new("matrix")),
        index=list(integer(0), integer(0)),
        N1=1L,
        N2=2L
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.is_valid_N <- function(N, ndim)
{
    isSingleInteger(N) && N >= 1L && N <= ndim
}

.validate_DelayedMatrix <- function(x)
{
    if (!.is_valid_N(x@N1, length(x@index))
     || !.is_valid_N(x@N2, length(x@index)))
        return(wmsg("'x@N1' and 'x@N2' must be single integers ",
                    ">= 1 and <= 'length(x@index)'"))
    if (x@N1 >= x@N2)
        return("'x@N1' must be < 'x@N2'")
    array_dim <- lengths(x@index)
    if (!all(array_dim[-c(x@N1, x@N2)] == 1L))
        return(wmsg("'x@N1' and 'x@N2' are incompatible with the ",
                    "dimensions of the underlying DelayedArray object"))
    TRUE
}

setValidity2("DelayedMatrix", .validate_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###
### Defining the internal index() getter and setter is enough to make all the
### DelayedArray accessors (length, isEmpty, dim, dimnames, dimnames<-) work
### on a DelayedMatrix object.

setMethod("index", "DelayedMatrix",
    function(x) x@index[c(x@N1, x@N2)]
)
setReplaceMethod("index", "DelayedMatrix",
    function(x, value) { x@index[c(x@N1, x@N2)] <- value; x }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.from_DelayedArray_to_DelayedMatrix <- function(from)
{
    from_dim <- dim(from)
    if (length(from_dim) < 2L)
        stop(wmsg(class(from), " object with less than 2 dimensions cannot ",
                  "be coerced to a DelayedMatrix object at the moment"))
    if (length(from_dim) == 2L) {
        N1 <- 1L
        N2 <- 2L
    } else {
        idx <- which(from_dim != 1L)
        if (length(idx) > 2L)
            stop(wmsg("Array-like object with more than 2 effective ",
                      "dimensions cannot be coerced to a DelayedMatrix ",
                      "object. ", slicing_tip))
        if (length(idx) == 2L) {
            N1 <- idx[[1L]]
            N2 <- idx[[2L]]
        } else if (length(idx) == 0L) {
            N1 <- 1L
            N2 <- 2L
        } else {
            ## length(idx) == 1L
            N1 <- idx[[1L]]
            if (N1 == length(from_dim))
                stop(wmsg("A ", class(from), " object where the only ",
                          "effective dimension is its last dimension cannot ",
                          "be coerced to a DelayedMatrix object at the ",
                          "moment"))
            N2 <- N1 + 1L
        }
    }
    new2("DelayedMatrix", from, N1=N1, N2=N2)
}

setAs("DelayedArray", "DelayedMatrix", .from_DelayedArray_to_DelayedMatrix)

### array -> DelayedMatrix

.from_array_to_DelayedMatrix <- function(from)
{
    as(as(from, "DelayedArray"), "DelayedMatrix")
}

setAs("array", "DelayedMatrix", .from_array_to_DelayedMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Specialized constructor
###

### 'x' must be an array-like object with 3 dimensions.
### 'MARGIN' is the dimension to drop.
make_DelayedMatrix_from_3D_DelayedArray <- function(x, MARGIN)
{
    if (!is(x, "DelayedArray"))
        x <- as(x, "DelayedArray")

    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (x_ndim != 3L)
        stop("'x' must have 3 dimensions")
    if (!isSingleNumber(MARGIN))
        stop("'MARGIN' must be a single integer")
    if (!is.integer(MARGIN))
        MARGIN <- as.integer(MARGIN)
    if (MARGIN < 1L || MARGIN > x_ndim)
        stop("'MARGIN' must be >= 1 and <= length(dim(x))")
    if (x_dim[[MARGIN]] != 1L)
        stop("'dim(x)[[MARGIN]]' must be 1")

    if (x@is_transposed)
        MARGIN <- x_ndim + 1L - MARGIN
    tmp <- seq_along(x_dim)[-MARGIN]
    N1 <- tmp[[1L]]
    N2 <- tmp[[2L]]
    new2("DelayedMatrix", x, N1=N1, N2=N2)
}

