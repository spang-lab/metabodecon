testthat::skip_on_cran()
testthat::skip_on_ci()

sap <- test_that("MetaboDecon1D works for single spectrum", {
    decon0 <- MetaboDecon1D_silent(
        filepath = metabodecon_file("bruker/sap"),
        filename = "sap_01",
        number_iterations = 3,
        range_water_signal_ppm = 0,
        signal_free_region = c(3.2, -3.2),
        smoothing_param = c(1, 3),
        delta = 3
    )
    expect_identical(object = names(decon0), expected = decon0_members_mandatory)
    expect_identical(object = class(decon0), expected = "list")
    decon2 <- as_decon2(decon0, spectrum = sap[[1]], sfr = c(3.2, -3.2), wshw = 0)
    obj2 <- calc_prarp(x = decon2, truepar = sap[[1]]$meta$simpar)
    expect_true(obj2$prarpx >= 0.507)
})

sim_subset <- test_that("MetaboDecon1D works for multiple spectra", {
    decons0 <- MetaboDecon1D_silent(
        filepath = metabodecon_file("bruker/sim_subset"),
        number_iterations = 3,
        range_water_signal_ppm = 0,
        signal_free_region = c(3.55, 3.35)
    )
    expect_identical(names(decons0), c("sim_01", "sim_02"))
    expect_identical(class(decons0), "list")
    expect_identical(names(decons0[[1]]), decon0_members)
    expect_identical(class(decons0[[1]]), "list")
    expect_identical(names(decons0[[2]]), decon0_members)
    expect_identical(class(decons0[[2]]), "list")
    decons2 <- as_decons2(decons0, spectra = sim[1:2])
    if (identical(environment(), .GlobalEnv)) {
        plot_spectrum(
            decons2[[1]],
            truepar = sim[[1]]$meta$simpar,
            sub1 = list(tp_verts = TRUE)
        )
        plot_spectrum(
            decons2[[2]],
            foc_frac = c(1, 0),
            truepar = sim[[2]]$meta$simpar,
            sub1 = list(tp_verts = TRUE),
            sub3 = FALSE
        )
    }
    obj1 <- calc_prarp(decons2[[1]], truepar = sim[[1]]$meta$simpar)
    obj2 <- calc_prarp(decons2[[2]], truepar = sim[[2]]$meta$simpar)
    expect_true(obj1$prarpx >= 0.732)
    expect_true(obj2$prarpx >= 0.710)
})

