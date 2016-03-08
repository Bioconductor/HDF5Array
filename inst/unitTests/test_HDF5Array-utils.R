
test_ArrayBlocks_class <- function()
{
    split_array_in_blocks <- HDF5Array:::split_array_in_blocks
    unsplit_array_from_blocks <- HDF5Array:::unsplit_array_from_blocks

    a1 <- array(1:300, c(10, 3, 2, 5))
    A1 <- as(a1, "HDF5Array")
    for (max_block_len in c(1:10, 19:20, 29:30, 39:40, 59:60, 119:120)) {
        subarrays <- split_array_in_blocks(a1, max_block_len)
        a1b <- unsplit_array_from_blocks(subarrays, a1)
        checkIdentical(a1, a1b)

        subarrays <- split_array_in_blocks(A1, max_block_len)
        A1b <- unsplit_array_from_blocks(subarrays, A1)
        checkIdentical(a1b, A1b)
    }

    m1 <- a1[c(9, 3:7), 2, 2, -4]
    M1 <- as(A1[c(9, 3:7), 2, 2, -4], "HDF5Matrix")
    checkIdentical(m1, as.matrix(M1))
    for (max_block_len in seq_len(length(m1) - 1L)) {
        subarrays <- split_array_in_blocks(m1, max_block_len)
        m1b <- unsplit_array_from_blocks(subarrays, m1)
        checkIdentical(m1, m1b)

        subarrays <- split_array_in_blocks(M1, max_block_len)
        M1b <- unsplit_array_from_blocks(subarrays, M1)
        checkIdentical(m1b, M1b)
    }
}

