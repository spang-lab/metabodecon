test_that("convert_spectrum works", {
    sim <- pkg_file("example_datasets/bruker/sim")
    sim_1 <- pkg_file("example_datasets/bruker/sim/sim_01")
    X <- read_spectrum(sim_1, raw = TRUE)
    Z <- convert_spectrum(X, sfx = 1e3, sfy = 1e6)
    expect_str(Z, c(
        "List of 14",
        " $ y_raw    : int [1:1309] 10151 10664 10496 9604 9056 9283 9401 8950 8822 9639 ...",
        " $ y_scaled : num [1:1309] 0.01015 0.01066 0.0105 0.0096 0.00906 ...",
        " $ n        : int 1309",
        " $ sfx      : num 1000",
        " $ sfy      : num 1e+06",
        " $ dp       : num [1:1309] 1308 1307 1306 1305 1304 ...",
        " $ sdp      : num [1:1309] 1.31 1.31 1.31 1.31 1.3 ...",
        " $ ppm      : num [1:1309] 3.6 3.6 3.6 3.6 3.6 ...",
        " $ hz       : num [1:1309] 6e+08 6e+08 6e+08 6e+08 6e+08 ...",
        " $ ppm_min  : num 3.4",
        " $ ppm_max  : num 3.6",
        " $ ppm_range: num 0.2",
        " $ ppm_step : num 0.000153",
        " $ ppm_nstep: num 0.000153"
    ))
})

test_that("convert_spectrum produces same output as MetaboDecon1D did", {
    sim <- pkg_file("example_datasets/bruker/sim")
    sim_1 <- pkg_file("example_datasets/bruker/sim/sim_01")
    X <- read_spectrum(sim_1, raw = TRUE)
    Z <- convert_spectrum(X, sfx = 1e3, sfy = 1e6)
    Y <- md1d("sim_01", nfit = 1, cache = FALSE)$rv$debuglist
    expect_equal(Z$y_raw, Y$data_read$spectrum_y_raw)
    expect_equal(Z$y_scaled, Y$data_read$spectrum_y)
})
