# Keep this test in a seperate file, because it's very slow. This is useful
# because testthat runs test files in parallel, i.e. by distributing slow tests
# across different files we ensure they can be also distributed across different
# cores.
library(testthat)

test_that("1. load_jcampdx_spectrum", {
    x <- evalwith(
        testdir = "load_jcampdx_spectrum/1",
        inputs = c(urine_1.dx = "jcampdx/urine/urine_1.dx"),
        expr = {
            spectrum_data <- load_jcampdx_spectrum("urine_1.dx")
        }
    )
    expect_identical(str2(spectrum_data), c("List of 7",
        " $ x                : num [1:131072] 131 131 131 131 131 ...",
        " $ y                : num [1:131072] 0.001265 0.001003 0.000105 -0.000937 -0.001062 ...",
        " $ x_ppm            : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        " $ length           : int 131072",
        " $ ppm_range        : num 20",
        " $ ppm_highest_value: num 14.8",
        " $ ppm_lowest_value : num -5.22"
    ))
    expect_equal(dir(x$testdir), c("urine_1.dx"))
})
