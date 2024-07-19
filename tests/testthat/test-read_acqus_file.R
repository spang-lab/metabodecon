library(testthat)

test_that("read_acqus_file works", {
    urine1_dir <- pkg_file("example_datasets/bruker/urine/urine_1")
    acqus <- read_acqus_file(urine1_dir)
    expect_equal(as.numeric(acqus$SW), 20.0236144338963)
    expect_equal(as.numeric(acqus$SW_h), 12019.2307692308)
})
