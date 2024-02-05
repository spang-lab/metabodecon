test_that("Deconvolution of 1 jcampdx works", {
    skip_if(Sys.getenv("SKIP_SLOW_TESTS") == "TRUE", "Skipped because SKIP_SLOW_TESTS=TRUE")

    x <- with(
        testdir = "MetaboDecon1D/1", answers = c("y", "y"),
        output = "out.txt", message = "err.txt", plots = "plots.pdf",
        inputs = c(urine.dx = "jcampdx/urine/urine.dx"),
        expr = {
            set.seed(123)
            MetaboDecon1D(filepath = ".", filename = "urine.dx", file_format = "jcampdx", number_iterations = 1)
        }
    )

    checksum_expected <- c(
        err.txt = 518, out.txt = 18, plots.pdf = 321364,
        urine.dx = 1192696,
        urine.dx_approximated_spectrum.txt = 2581877,
        urine.dx_parameters.txt = 72100
    )
    expect_identical(checksum(x$testdir), checksum_expected)
})
