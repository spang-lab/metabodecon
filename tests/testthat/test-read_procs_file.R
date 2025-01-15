library(testthat)

test_that("read_procs_file works", {
    spldir <- pkg_file("example_datasets/bruker/urine/urine_1")
    procs <- read_procs_file(spldir)
    expect_equal(length(procs), 125)
    expect_equal(procs$TITLE, "Parameter file, TOPSPIN\t\tVersion 3.1")
    expect_equal(procs$JCAMPDX, 5)
    expect_equal(procs$BYTORDP, 0)
    expect_equal(procs$NC_proc, -2)
    expect_equal(procs$DTYPP, 0)
    expect_equal(procs$SI, 131072)
})
