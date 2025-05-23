library(testthat)

test_that("metabodecon_file works", {
    # Unambiguous paths
    expect_length(metabodecon_file("urine_1"), 1)
    expect_length(metabodecon_file("urine_1.dx"), 1)
    expect_length(metabodecon_file("sim/sim_01"), 1)

    # Ambiguous paths (i.e. multiple matches)
    expect_length(metabodecon_file("sim_01"), 2)
    expect_length(metabodecon_file("urine"), 2)

    # Non-existing paths (i.e. a character vector of length zero gets returned)
    expect_length(metabodecon_file("asdfasdf"), 0)
})
