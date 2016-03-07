
test_ArrayBlocks_class <- function()
{
    ArrayBlocks <- HDF5Array:::ArrayBlocks
    break_array_in_blocks <- HDF5Array:::break_array_in_blocks
    rebuild_array_from_blocks <- HDF5Array:::rebuild_array_from_blocks

    a1 <- array(1:600, c(10, 3:5))
    for (max_block_len in 1:599) {
        subarrays <- break_array_in_blocks(a1, max_block_len)
        a1b <- rebuild_array_from_blocks(subarrays, a1)
        checkIdentical(a1, a1b)
    }
}

