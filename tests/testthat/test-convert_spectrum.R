test_that("convert_spectrum works", {
    urine_1 <- pkg_file("example_datasets/bruker/urine/urine_1")
    X <- read_spectrum(urine_1, raw = TRUE)
    Z <- convert_spectrum(X, sfx = 1e3, sfy = 1e6)
    Y <- MD1D()$rv$debuglist
    expect_equal(Z$y_raw, Y$data_read$spectrum_y_raw)
    expect_equal(Z$y_scaled, Y$data_read$spectrum_y)
    expect_identical(str2(Z, digits.d = 12, vec.len = 1), c(
        "List of 14",
        " $ y_raw    : int [1:131072] 1265 1003 ...",
        " $ y_scaled : num [1:131072] 0.001265 0.001003 ...",
        " $ n        : int 131072",
        " $ sfx      : num 1000",
        " $ sfy      : num 1e+06",
        " $ dp       : num [1:131072] 131071 131070 ...",
        " $ sdp      : num [1:131072] 131.071 131.07 ...",
        " $ ppm      : num [1:131072] 14.80254 ...",
        " $ fq       : num [1:131072] 600243921.684 ...",
        " $ ppm_min  : num -5.2210744339",
        " $ ppm_max  : num 14.80254",
        " $ ppm_range: num 20.0236144339",
        " $ ppm_step : num 0.000152769219994",
        " $ ppm_nstep: num 0.000152768054458"
    ))
})
