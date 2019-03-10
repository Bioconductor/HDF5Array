# This tests the ability of the API to properly access integer HDF5Matrices.
# library(testthat); source("setup.R"); source("test-integer.R")

hFUN <- integer_hFUN

set.seed(34567)
test_that("HDF5 integer matrix input is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_read_all(hFUN, mode="integer")
    check_read_all(hFUN, nr=5, nc=30, mode="integer")
    check_read_all(hFUN, nr=30, nc=5, mode="integer")

    check_read_slice(hFUN, mode="integer")
    check_read_slice(hFUN, nr=5, nc=30, mode="integer")
    check_read_slice(hFUN, nr=30, nc=5, mode="integer")

    check_read_varslice(hFUN, mode="integer")
    check_read_varslice(hFUN, nr=5, nc=30, mode="integer")
    check_read_varslice(hFUN, nr=30, nc=5, mode="integer")

    check_read_multi(hFUN, mode="integer")
    check_read_multi(hFUN, nr=5, nc=30, mode="integer")
    check_read_multi(hFUN, nr=30, nc=5, mode="integer")

    check_read_type(hFUN, mode="integer")
    check_read_class(hFUN(), mode="integer", "HDF5Matrix")

    check_read_errors(hFUN, mode="integer")
    check_read_all(hFUN, nr=0, nc=0, mode="integer")
    check_read_all(hFUN, nr=10, nc=0, mode="integer")
    check_read_all(hFUN, nr=0, nc=10, mode="integer")
})

set.seed(34568)
test_that("HDF5 integer matrix output is okay", {
    expect_s4_class(hFUN(), "HDF5Matrix")

    check_write_all(hFUN, mode="integer")
    check_write_all(hFUN, nr=5, nc=30, mode="integer")
    check_write_all(hFUN, nr=30, nc=5, mode="integer")

    check_write_slice(hFUN, mode="integer")
    check_write_slice(hFUN, nr=5, nc=30, mode="integer")
    check_write_slice(hFUN, nr=30, nc=5, mode="integer")

    check_write_varslice(hFUN, mode="integer")
    check_write_varslice(hFUN, nr=5, nc=30, mode="integer")
    check_write_varslice(hFUN, nr=30, nc=5, mode="integer")

    check_write_indexed(hFUN, mode="integer")
    check_write_indexed(hFUN, nr=5, nc=30, mode="integer")
    check_write_indexed(hFUN, nr=30, nc=5, mode="integer")

    check_write_type(hFUN, mode="integer")
    check_write_errors(hFUN, mode="integer")

    check_write_HDF5(hFUN, mode="integer")

    check_write_all(hFUN, nr=0, nc=0, mode="integer")
    check_write_all(hFUN, nr=10, nc=0, mode="integer")
    check_write_all(hFUN, nr=0, nc=10, mode="integer")
})
