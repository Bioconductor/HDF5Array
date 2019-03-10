# This tests the ability of the API to properly access logical HDF5Matrices.
# library(testthat); source("setup.R"); source("test-logical.R")

hFUN <- logical_hFUN

set.seed(34567)
test_that("HDF5 logical matrix input is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_read_all(hFUN, mode="logical")
    check_read_all(hFUN, nr=5, nc=30, mode="logical")
    check_read_all(hFUN, nr=30, nc=5, mode="logical")

    check_read_slice(hFUN, mode="logical")
    check_read_slice(hFUN, nr=5, nc=30, mode="logical")
    check_read_slice(hFUN, nr=30, nc=5, mode="logical")

    check_read_varslice(hFUN, mode="logical")
    check_read_varslice(hFUN, nr=5, nc=30, mode="logical")
    check_read_varslice(hFUN, nr=30, nc=5, mode="logical")

    check_read_multi(hFUN, mode="logical")
    check_read_multi(hFUN, nr=5, nc=30, mode="logical")
    check_read_multi(hFUN, nr=30, nc=5, mode="logical")

    check_read_type(hFUN, mode="logical")
    check_read_class(hFUN(), mode="logical", "HDF5Matrix")

    check_read_errors(hFUN, mode="logical")
    check_read_all(hFUN, nr=0, nc=0, mode="logical")
    check_read_all(hFUN, nr=10, nc=0, mode="logical")
    check_read_all(hFUN, nr=0, nc=10, mode="logical")
})

set.seed(34568)
test_that("HDF5 logical matrix output is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_write_all(hFUN, mode="logical")
    check_write_all(hFUN, nr=5, nc=30, mode="logical")
    check_write_all(hFUN, nr=30, nc=5, mode="logical")

    check_write_slice(hFUN, mode="logical")
    check_write_slice(hFUN, nr=5, nc=30, mode="logical")
    check_write_slice(hFUN, nr=30, nc=5, mode="logical")

    check_write_varslice(hFUN, mode="logical")
    check_write_varslice(hFUN, nr=5, nc=30, mode="logical")
    check_write_varslice(hFUN, nr=30, nc=5, mode="logical")

    check_write_indexed(hFUN, mode="logical")
    check_write_indexed(hFUN, nr=5, nc=30, mode="logical")
    check_write_indexed(hFUN, nr=30, nc=5, mode="logical")

    check_write_type(hFUN, mode="logical")
    check_write_errors(hFUN, mode="logical")

    check_write_HDF5(hFUN, mode="logical")

    check_write_all(hFUN, nr=0, nc=0, mode="logical")
    check_write_all(hFUN, nr=10, nc=0, mode="logical")
    check_write_all(hFUN, nr=0, nc=10, mode="logical")
})
