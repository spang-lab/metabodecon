make_decon1 <- function(show = TRUE) {
    cs <- seq(1.2, 0.2, by = -0.1)
    x <- simulate_spectrum("utest", cs = cs, x0 = c(0.7), A = 2.5)
    if (show) plot_spectrum(x, sub_show = FALSE)
    x$si <- round(x$si, 3)
    decon1 <- generate_lorentz_curves(x, sfr = c(0.9, 0.5), wshw = 0, smopts = c(0, 3), ask = FALSE)
    if (show) plot_spectrum(decon1, sub_show = FALSE)
    dput(decon1)
}

decon0_exp <- list(
    number_of_files = 1L,
    filename = "utest",
    x_values = c(0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001, 0),
    x_values_ppm = c(1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
    y_values = c(0.001553, 0.000263, 0.001942, 8.9e-05, 1e-08, 0.02249, 0.001455, 2e-05, 0.000335, 0.00015, 0.000701),
    spectrum_superposition = structure(c(6.5573545883911e-10, 9.75607020027042e-10, 1.59999712260121e-09, 3.07692334450644e-09, 8.00003725816359e-09, 4.00008603255252e-08, 3.99998317941017e-08, 7.99991383728281e-09, 3.07689291540482e-09, 1.59998560337295e-09, 9.75601513472325e-10), dim = c(1L, 11L)),
    mse_normed = 0.0257992661017059,
    index_peak_triplets_middle = 6,
    index_peak_triplets_left = 7,
    index_peak_triplets_right = 5,
    peak_triplets_middle = 0.7,
    peak_triplets_left = 0.6,
    peak_triplets_right = 0.8,
    integrals = structure(-1.18167676229858e-10, dim = c(1L, 1L)),
    signal_free_region = c(0.0087, 0.0043),
    range_water_signal_ppm = 0,
    A = 4.00003460607913e-11,
    lambda = -0.000499992682537264,
    x_0 = 0.00450000642817171
)

decon1 <- structure(list(
    number_of_files = 1L,
    filename = "utest",
    x_values = c(0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001, 0),
    x_values_ppm = c(1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
    y_values = c(0.001553, 0.000263, 0.001942, 8.9e-05, 1e-08, 0.02249, 0.001455, 2e-05, 0.000335, 0.00015, 0.000701),
    spectrum_superposition = structure(c(6.5573545883911e-10, 9.75607020027042e-10, 1.59999712260121e-09, 3.07692334450644e-09, 8.00003725816359e-09, 4.00008603255252e-08, 3.99998317941017e-08, 7.99991383728281e-09, 3.07689291540482e-09, 1.59998560337295e-09, 9.75601513472325e-10), dim = c(1L, 11L)),
    mse_normed = 0.0257992661017059,
    index_peak_triplets_middle = 6,
    index_peak_triplets_left = 7,
    index_peak_triplets_right = 5,
    peak_triplets_middle = 0.7,
    peak_triplets_left = 0.6,
    peak_triplets_right = 0.8,
    integrals = structure(-1.18167676229858e-10, dim = c(1L, 1L)),
    signal_free_region = c(0.0087, 0.0043),
    range_water_signal_ppm = 0,
    A = 4.00003460607913e-11,
    lambda = -0.000499992682537264,
    x_0 = 0.00450000642817171,
    y_values_raw = c(1553, -263, 1942, 89, 316, 22490, -1455, 20, 335, 150, 701),
    x_values_hz = c(600252086.646632, 600252146.671912, 600252206.697193, 600252266.722474, 600252326.747754, 600252386.773035, 600252446.798316, 600252506.823597, 600252566.848877, 600252626.874158, 600252686.899439),
    mse_normed_raw = 0.0407322097888289,
    signal_free_region_ppm = c(0.9, 0.5),
    x_0_hz = 600252416.78529,
    x_0_dp = 4.50000642817171,
    x_0_ppm = 0.650000642817171,
    A_hz = -2.40103200000386e-06,
    A_dp = 4.00003460607913e-08,
    A_ppm = 4.00003460607913e-09,
    lambda_hz = 30.0122011123419,
    lambda_dp = -0.499992682537264,
    lambda_ppm = -0.0499992682537264
), class = "decon1")

# CONTINUE HERE

decon2 <- structure(list(
    cs = c(1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
    si = c(1553, -263, 1942, 89, 316, 22490, -1455, 20, 335, 150, 701),
    meta = list(
        name = "utest",
        path = NULL,
        type = NULL,
        fq = c(600252086.646632, 600252146.671912, 600252206.697193, 600252266.722474, 600252326.747754, 600252386.773035, 600252446.798316, 600252506.823597, 600252566.848877, 600252626.874158, 600252686.899439),
        mfs = NULL,
        simpar = NULL
    ), args = list(
        nfit = NA,
        smopts = NA,
        delta = NA,
        sfr = c(0.0087, 0.0043),
        wshw = 0
    ), sit = list(
        wsrm = NA,
        nvrm = NA,
        sm = c(0.001553, 0.000263, 0.001942, 8.9e-05, 1e-08, 0.02249, 0.001455, 2e-05, 0.000335, 0.00015, 0.000701),
        sup = structure(c(6.5573545883911e-10, 9.75607020027042e-10, 1.59999712260121e-09, 3.07692334450644e-09, 8.00003725816359e-09, 4.00008603255252e-08, 3.99998317941017e-08, 7.99991383728281e-09, 3.07689291540482e-09, 1.59998560337295e-09, 9.75601513472325e-10), dim = c(1L, 11L)),
        al = NULL
    ), peak = list(
        center = 6,
        left = 7,
        right = 5
    ), lcpar = list(
        A = 4.00003460607913e-11,
        lambda = -0.000499992682537264,
        x0 = 0.00450000642817171
    ),
    mse = list(
        raw = NA,
        norm = 0.0407322097888289,
        sm = NA,
        smnorm = 0.0257992661017059
    )
), class = "decon2")

idecon <- structure(list(), class = "idecon") # TODO

result_from_decon1 <- testthat::test_that("as_decon0 works", {
    decon0 <- as_decon0(x = decon1)
    testthat::expect_identical(decon0, decon0_exp)
})

result_from_decon2 <- testthat::test_that("as_decon0 works", {
    decon0 <- as_decon0(x = decon2)
    testthat::expect_identical(decon0, decon0_exp)
})
