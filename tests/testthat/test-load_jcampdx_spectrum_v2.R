# Keep this test in a seperate file, because it's very slow. This is useful
# because testthat runs test files in parallel, i.e. by distributing slow tests
# across different files we ensure they can be also distributed across different
# cores.
library(testthat)

test_that("1. load_jcampdx_spectrum_v1", {
    x <- with(
        testdir = "load_jcampdx_spectrum_v1/1",
        inputs = c(urine_1.dx = "jcampdx/urine/urine_1.dx"),
        expr = {
            spectrum_data <- load_jcampdx_spectrum_v1("urine_1.dx")
        }
    )
    expect_identical(str2(spectrum_data), c(
        "List of 19",
        " $ ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        " $ dp               : num [1:131072] 131071 131070 131069 131068 131067 ...",
        " $ sdp              : num [1:131072] 131 131 131 131 131 ...",
        " $ ss               : int [1:131072] 1265 1003 105 -937 -1062 -93 701 515 106 110 ...",
        " $ y              : num [1:131072] 0.001265 0.001003 0.000105 -0.000937 -0.001062 ...",
        " $ n                : int 131072",
        " $ sfx              : num 1000",
        " $ sfy              : num 1e+06",
        " $ ppm_min          : num -5.22",
        " $ ppm_max          : num 14.8",
        " $ ppm_range        : num 20",
        " $ ppm_step         : num 0.000153",
        " $ ppm_nstep        : num 0.000153",
        " $ length           : int 131072",
        " $ x_ppm            : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...",
        " $ x                : num [1:131072] 131 131 131 131 131 ...",
        " $ y                : num [1:131072] 0.001265 0.001003 0.000105 -0.000937 -0.001062 ...",
        " $ ppm_highest_value: num 14.8",
        " $ ppm_lowest_value : num -5.22"
    ))
    expect_identical(lapply(spectrum_data, list(
        ppm = c(14.80254, 14.80238723078, 14.80223446156, 14.80208169234, 14.80192892312, 14.8017761539),
        dp = c(131071, 131070, 131069, 131068, 131067, 131066),
        sdp = c(131.071, 131.07, 131.069, 131.068, 131.067, 131.066),
        ss = c(1265L, 1003L, 105L, -937L, -1062L, -93L),
        y = c(0.001265, 0.001003, 0.000105, -0.000937, -0.001062, -9.3e-05),
        n = 131072L,
        sfx = 1000,
        sfy = 1e+06,
        ppm_min = -5.2210744338963,
        ppm_max = 14.80254,
        ppm_range = 20.0236144338963,
        ppm_step = 0.000152769219994479,
        ppm_nstep = 0.000152768054457827,
        length = 131072L,
        x_ppm = c(14.80254, 14.80238723078, 14.80223446156, 14.80208169234, 14.80192892312, 14.8017761539),
        x = c(131.071, 131.07, 131.069, 131.068, 131.067, 131.066),
        y = c(0.001265, 0.001003, 0.000105, -0.000937, -0.001062, -9.3e-05),
        ppm_highest_value = 14.80254,
        ppm_lowest_value = -5.2210744338963
    )))
    expect_equal(dir(x$testdir), c("urine_1.dx"))
})
