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
### on an DelayedMatrix object.

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
            stop(wmsg(class(from), " object with more than 2 effective ",
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

### matrix -> DelayedMatrix

.from_matrix_to_DelayedMatrix <- function(from)
{
    as(as(from, "DelayedArray"), "DelayedMatrix")
}

setAs("matrix", "DelayedMatrix", .from_matrix_to_DelayedMatrix)


