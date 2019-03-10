# Tests whether HDF5Matrix can be successfully read by beachmat for strings.
# library(testthat); source("setup.R"); source("test-character.R")

hFUN <- character_hFUN

test_that("HDF5Matrix can be read by beachmat", {
    check_read_all(hFUN, mode="character")
    check_read_all(hFUN, nr=5, nc=30, mode="character")
    check_read_all(hFUN, nr=30, nc=5, mode="character")

    check_read_slice(hFUN, mode="character")
    check_read_slice(hFUN, nr=5, nc=30, mode="character")
    check_read_slice(hFUN, nr=30, nc=5, mode="character")

    check_read_varslice(hFUN, mode="character")
    check_read_varslice(hFUN, nr=5, nc=30, mode="character")
    check_read_varslice(hFUN, nr=30, nc=5, mode="character")

    check_read_multi(hFUN, mode="character")
    check_read_multi(hFUN, nr=5, nc=30, mode="character")
    check_read_multi(hFUN, nr=30, nc=5, mode="character")

    check_read_class(hFUN(), mode="character", "HDF5Matrix")

    check_read_errors(hFUN, mode="character")
    check_read_all(hFUN, nr=0, nc=0, mode="character")
    check_read_all(hFUN, nr=10, nc=0, mode="character")
    check_read_all(hFUN, nr=0, nc=10, mode="character")
})

test_that("HDF5Matrix can be written by beachmat", {
    check_write_all(hFUN, mode="character")
    check_write_all(hFUN, nr=5, nc=30, mode="character")
    check_write_all(hFUN, nr=30, nc=5, mode="character")

    check_write_slice(hFUN, mode="character")
    check_write_slice(hFUN, nr=5, nc=30, mode="character")
    check_write_slice(hFUN, nr=30, nc=5, mode="character")

    check_write_varslice(hFUN, mode="character")
    check_write_varslice(hFUN, nr=5, nc=30, mode="character")
    check_write_varslice(hFUN, nr=30, nc=5, mode="character")

    check_write_indexed(hFUN, mode="character")
    check_write_indexed(hFUN, nr=5, nc=30, mode="character")
    check_write_indexed(hFUN, nr=30, nc=5, mode="character")

    check_write_errors(hFUN, mode="character")

    check_write_all(hFUN, nr=0, nc=0, mode="character")
    check_write_all(hFUN, nr=10, nc=0, mode="character")
    check_write_all(hFUN, nr=0, nc=10, mode="character")
})
