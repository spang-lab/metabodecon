library(testthat)

test_that("align works", {

    # Prepare inputs
    sap_01 <- sap[[1]]
    sap_01_shifted <- simulate_spectrum(
        name = "sap_01_shifted",
        cs = sap_01$meta$simpar$cs,
        x0 = sap_01$meta$simpar$x0 + 0.3,
        A  = sap_01$meta$simpar$A,
        lambda = sap_01$meta$simpar$lambda,
        noise = sap_01$meta$simpar$noise
    )
    spectra <- as_spectra(list(sap_01, sap_01_shifted))
    plot_spectra(spectra)
    decons <- deconvolute(spectra, smopts=c(1,3), delta=3, sfr=c(3.2,-3.2))

    # Do the alignment
    aligns <- align(decons)

    # Check structure of returned object. Strategy: add all fields to the
    # decons object that we expect [align()] to add. At the end the objects
    # should be equal.
    decons_copy <- decons
    for (i in seq_along(aligns)) {
        decons_copy[[i]]$sit$al <- aligns[[i]]$sit$al
        decons_copy[[i]]$sit$supal <- aligns[[i]]$sit$supal
        decons_copy[[i]]$lcpar$x0_al <- aligns[[i]]$lcpar$x0_al
        class(decons_copy[[i]]) <- "align"
    }
    class(decons_copy) <- "aligns"
    expect_equal(object = aligns, expected = decons_copy)

    # Check that the alignment worked, our expectations are:
    # 1. x0_al     is shifted rougly 0.3 to the right compared to x0
    # 2. sit$al    is equal to the integral at peak-center-datapoints
    # 3. sit$al    is zero at non-peak-center-datapoints
    # 4. sit$supal is the superposition the aligned Lorentz Curves
    x0 <- aligns$sap_01_shifted$lcpar$x0
    x0_al <- aligns$sap_01_shifted$lcpar$x0_al
    shifts <- x0 - x0_al
    expect_true(all(shifts > 0.2 & shifts < 0.4))

    cs <- aligns$sap_01_shifted$cs
    dp <- seq_along(cs)
    pc <- match(x0_al, cs)
    non_pc <- setdiff(dp, pc)
    A <- aligns$sap_01_shifted$lcpar$A
    al <- aligns$sap_01_shifted$sit$al
    expect_equal(al[pc], A * pi)
    expect_equal(al[non_pc], rep(0, length(non_pc)))

    supal <- aligns$sap_01_shifted$sit$supal
    lambda <- aligns$sap_01_shifted$lcpar$lambda
    expect_equal(supal, lorentz_sup(cs, x0_al, A, lambda))
})
