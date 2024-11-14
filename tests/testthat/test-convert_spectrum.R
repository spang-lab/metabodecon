test_that("as_ispec works", {
    sim <- pkg_file("example_datasets/bruker/sim")
    sim_1 <- pkg_file("example_datasets/bruker/sim/sim_01")
    x <- read_spectrum(sim_1, raw = TRUE)
    x <- as_ispec(x, sf = c(1e3, 1e6))
    expect_true(all(ispec_members %in% names(x)))
})

test_that("as_ispec produces same output as MetaboDecon1D did", {
    sim <- pkg_file("example_datasets/bruker/sim")
    sim_1 <- pkg_file("example_datasets/bruker/sim/sim_01")
    x <- read_spectrum(sim_1, raw = TRUE)
    x <- as_ispec(x, sf = c(1e3, 1e6))
    y <- MetaboDecon1D_silent_sim(sim, "sim_01", debug = TRUE)
    expect_equal(x$y_raw, y$debuglist$data$spectrum_y_raw)
    expect_equal(x$y_scaled, y$debuglist$data$spectrum_y)
})
