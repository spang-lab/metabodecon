sap1 <- get_sap1()

check_decon1 <- test_that("sap1$decon1 == generate_lorentz_curves(sap1$spectrum)", {
    decon1 <- generate_lorentz_curves(
        data_path = sap1$spectrum,
        sfr = c(1, -1),
        wshw = 0,
        smopts = c(0, 3),
        ask = FALSE,
        verbose = FALSE,
        force = TRUE
    )
    if (identical(environment(), globalenv())) {
        str(decon1)
        plot_spectrum(decon1, sub_show = FALSE)
    }
    expect_equal(decon1, sap1$decon1)
})
