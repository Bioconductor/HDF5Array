.make_TEST_H5SparseMatrixSeeds <- function(m)
{
    tenx <- seed(writeTENxMatrix(m))
    csc <- H5SparseMatrixSeed(tenx@filepath, tenx@group)
    list(tenx=tenx, csc=csc, csr=t(csc))
}

### Number of calls to h5mread::h5mread() made by evaluating 'expr'.
.count_h5mread_calls <- function(expr)
{
    calls <- new.env()
    calls$n <- 0L
    suppressMessages(trace("h5mread", where=asNamespace("h5mread"),
                           print=FALSE,
                           tracer=substitute(assign("n", E$n + 1L, envir=E),
                                             list(E=calls))))
    on.exit(suppressMessages(untrace("h5mread",
                                     where=asNamespace("h5mread"))))
    expr
    calls$n
}

test_H5SparseMatrixSeed_type <- function()
{
    set.seed(123)
    m <- matrix(rpois(200, 0.5), nrow=10, ncol=20)
    m[ , 5:7] <- 0L
    ## integer and double data
    for (x in list(m, m + 0.5 * (m != 0))) {
        for (s in .make_TEST_H5SparseMatrixSeeds(x)) {
            checkTrue(!is.na(s@type))
            checkIdentical(type(s), typeof(x))
            checkIdentical(.count_h5mread_calls(type(s)), 0L)
            checkIdentical(.count_h5mread_calls(
                               S4Arrays:::extract_empty_array(s)), 0L)
            checkIdentical(type(S4Arrays:::extract_empty_array(s)), typeof(x))

            ## selections without nonzero data don't access the file
            xs <- if (is(s, "CSR_H5SparseMatrixSeed")) t(x) else x
            if (is(s, "CSR_H5SparseMatrixSeed")) {
                index <- list(5:7, NULL)
            } else {
                index <- list(NULL, 5:7)
            }
            checkIdentical(.count_h5mread_calls(extract_sparse_array(s,
                                                                index)), 0L)
            checkIdentical(as.array(extract_sparse_array(s, index)),
                           xs[index[[1L]] %||% seq_len(nrow(xs)),
                              index[[2L]] %||% seq_len(ncol(xs)), drop=FALSE])
            checkIdentical(unname(as.matrix(DelayedArray(s))), unname(xs))

            ## objects without the "type" slot, e.g. serialized with an
            ## older version of HDF5Array, obtain the type from the file
            old <- s
            attr(old, "type") <- NULL
            checkTrue(!.hasSlot(old, "type"))
            checkIdentical(type(old), typeof(x))
            checkIdentical(S4Arrays:::extract_empty_array(old),
                           S4Arrays:::extract_empty_array(s))
            checkIdentical(as.array(extract_sparse_array(old, index)),
                           as.array(extract_sparse_array(s, index)))
        }
    }

    ## binding seeds and setting their dimnames, which validates the
    ## resulting DelayedArray objects, doesn't access the files
    s1 <- .make_TEST_H5SparseMatrixSeeds(m[ , 1:8])$csc
    s2 <- .make_TEST_H5SparseMatrixSeeds(m[ , 9:20])$csc
    n <- .count_h5mread_calls({
        B <- cbind(DelayedArray(s1), DelayedArray(s2))
        dimnames(B) <- list(letters[1:10], LETTERS[1:20])
        type(B)
    })
    checkIdentical(n, 0L)
    checkIdentical(unname(as.matrix(B)), m)
}

test_H5SparseMatrixSeed_subdata <- function()
{
    set.seed(123)
    m <- Matrix::rsparsematrix(10, 20, 0.3)
    m@x <- round(10 * m@x)
    ## the nonzero data are stored in datasets "counts" (double) and
    ## "flags" (integer) of the HDF5 group "g/data"
    filepath <- tempfile(fileext=".h5")
    rhdf5::h5createFile(filepath)
    rhdf5::h5createGroup(filepath, "g")
    rhdf5::h5createGroup(filepath, "g/data")
    rhdf5::h5write(m@x, filepath, "g/data/counts")
    rhdf5::h5write(as.integer(m@x > 0), filepath, "g/data/flags")
    rhdf5::h5write(m@i, filepath, "g/indices")
    rhdf5::h5write(m@p, filepath, "g/indptr")
    rhdf5::h5write(dim(m), filepath, "g/shape")
    flags <- as.matrix(m) > 0
    storage.mode(flags) <- "integer"
    expected <- list(counts=as.matrix(m), flags=flags)

    for (subdata in names(expected)) {
        s <- H5SparseMatrixSeed(filepath, "g", subdata=subdata,
                                sparse.layout="CSC")
        checkIdentical(s@subdata, subdata)
        checkIdentical(type(s), typeof(expected[[subdata]]))
        checkIdentical(nzcount(s), length(m@x))
        checkIdentical(as.matrix(DelayedArray(s)), expected[[subdata]])
        checkIdentical(as.array(extract_sparse_array(s, list(NULL, 2:4))),
                       expected[[subdata]][ , 2:4])

        ## objects without the "type" slot obtain it from the same dataset
        old <- s
        attr(old, "type") <- NULL
        checkIdentical(type(old), typeof(expected[[subdata]]))
    }
}
