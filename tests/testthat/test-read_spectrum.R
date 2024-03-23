expected_str_urine_1 <- c( # for str(x$rv, digits.d=12, vec.len=1)
    "List of 13",
    " $ y_raw    : num [1:131072] 1265 1003 ...",
    " $ y_scaled : num [1:131072] 0.001265 0.001003 ...",
    " $ n        : int 131072",
    " $ sfx      : num 1000",
    " $ sfy      : num 1e+06",
    " $ dp       : num [1:131072] 131071 131070 ...",
    " $ sdp      : num [1:131072] 131.071 131.07 ...",
    " $ ppm      : num [1:131072] 14.80254 ...",
    " $ ppm_min  : num -5.2210744339",
    " $ ppm_max  : num 14.80254",
    " $ ppm_range: num 20.0236144339",
    " $ ppm_step : num 0.000152769219994",
    " $ ppm_nstep: num 0.000152768054458"
)

test_that("read_spectrum works for bruker/urine_1", {
    x <- evalwith(
        testdir = "read_spectrum_urine1",
        inputs = "bruker/urine/urine_1",
        expr = read_spectrum("urine_1")
    )
    expect_identical(str2(x$rv, digits.d = 12, vec.len = 1), expected_str_urine_1)
})

test_that("read_spectrum works for jcampdx/urine_1.dx", {
    x <- evalwith(
        testdir = "read_spectrum_urine1dx",
        inputs = "jcampdx/urine/urine_1.dx",
        expr = read_spectrum("urine_1.dx", type = "jcampdx")
    )
    expect_identical(str2(x$rv, digits.d = 12, vec.len = 1), expected_str_urine_1)
})
