test_that("as_gspec works", {
    sim <- pkg_file("example_datasets/bruker/sim")
    sim_1 <- pkg_file("example_datasets/bruker/sim/sim_01")
    X <- read_spectrum(sim_1, raw = TRUE)
    Z <- as_gspec(X, sf = c(1e3, 1e6))
    expect_true(all(gspec_members %in% names(Z)))
})

test_that("as_gspec produces same output as MetaboDecon1D did", {
    sim <- pkg_file("example_datasets/bruker/sim")
    sim_1 <- pkg_file("example_datasets/bruker/sim/sim_01")
    X <- read_spectrum(sim_1, raw = TRUE)
    Z <- as_gspec(X, sf = c(1e3, 1e6))
    Y <- md1d("sim_01", nfit = 1, cache = FALSE)$rv$debuglist
    expect_equal(Z$y_raw, Y$data_read$spectrum_y_raw)
    expect_equal(Z$y_scaled, Y$data_read$spectrum_y)
})
