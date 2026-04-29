library(testthat)

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
decons <- deconvolute(
    spectra,
    smopts = c(1, 3),
    delta = 3,
    sfr = c(3.2, -3.2),
    verbose = FALSE
)

test_that("align works", {

    skip_if_speaq_deps_missing()

    aligns <- align(decons, verbose = FALSE)

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

test_that("align gives same result for 1 vs multiple workers", {

    skip_if_speaq_deps_missing()

    sap_01_shifted_2 <- simulate_spectrum(
        name = "sap_01_shifted_2",
        cs = sap_01$meta$simpar$cs,
        x0 = sap_01$meta$simpar$x0 + 0.15,
        A = sap_01$meta$simpar$A,
        lambda = sap_01$meta$simpar$lambda,
        noise = sap_01$meta$simpar$noise
    )
    spectra3 <- as_spectra(list(sap_01, sap_01_shifted, sap_01_shifted_2))
    decons3 <- deconvolute(
        spectra3,
        smopts = c(1, 3),
        delta = 3,
        sfr = c(3.2, -3.2),
        verbose = FALSE
    )

    al1 <- align(decons3, verbose = FALSE, nworkers = 1)
    al2 <- align(decons3, verbose = FALSE, nworkers = 2)
    expect_equal(al2, al1)
})

test_that("built-in backend matches speaq backend", {

    skip_if_speaq_deps_missing()

    al_builtin <- align(decons, verbose = FALSE, method = 2, .internal = TRUE)
    al_speaq <- align(decons, verbose = FALSE, method = 1, .internal = TRUE)
    expect_equal(al_builtin, al_speaq)
})

test_that("fast peak backend aligns simple shifts", {

    al_fast <- align(decons, verbose = FALSE, method = 3)

    x0 <- al_fast$sap_01_shifted$lcpar$x0
    x0_al <- al_fast$sap_01_shifted$lcpar$x0_al
    shifts <- x0 - x0_al
    expect_true(all(shifts > 0.2 & shifts < 0.4))

    cs <- al_fast$sap_01_shifted$cs
    dp <- seq_along(cs)
    pc <- match(x0_al, cs)
    non_pc <- setdiff(dp, pc)
    A <- al_fast$sap_01_shifted$lcpar$A
    al <- al_fast$sap_01_shifted$sit$al
    expect_equal(al[pc], A * pi)
    expect_equal(al[non_pc], rep(0, length(non_pc)))

    supal <- al_fast$sap_01_shifted$sit$supal
    lambda <- al_fast$sap_01_shifted$lcpar$lambda
    expect_equal(supal, lorentz_sup(cs, x0_al, A, lambda))
})

test_that("fast peak backend gives same result for 1 vs multiple workers", {

    sap_01_shifted_2 <- simulate_spectrum(
        name = "sap_01_shifted_2",
        cs = sap_01$meta$simpar$cs,
        x0 = sap_01$meta$simpar$x0 + 0.15,
        A = sap_01$meta$simpar$A,
        lambda = sap_01$meta$simpar$lambda,
        noise = sap_01$meta$simpar$noise
    )
    spectra3 <- as_spectra(list(sap_01, sap_01_shifted, sap_01_shifted_2))
    decons3 <- deconvolute(
        spectra3,
        smopts = c(1, 3),
        delta = 3,
        sfr = c(3.2, -3.2),
        verbose = FALSE
    )

    al1 <- align(decons3, verbose = FALSE, method = 3, nworkers = 1)
    al2 <- align(decons3, verbose = FALSE, method = 3, nworkers = 2)
    expect_equal(al2, al1)
})

test_that("align raises error for invalid method", {
    expect_error(align(decons, verbose = FALSE, method = 99, .internal = TRUE),
                 "`method` must be 1, 2 or 3")
})

test_that("align with full=FALSE omits supal", {
    skip_if_speaq_deps_missing()
    al_nofull <- align(decons, verbose = FALSE, full = FALSE, .internal = TRUE)
    for (i in seq_along(al_nofull)) {
        expect_null(al_nofull[[i]]$sit$supal)
    }
    # full=TRUE (default) should include supal
    al_full <- align(decons, verbose = FALSE, full = TRUE, .internal = TRUE)
    for (i in seq_along(al_full)) {
        expect_false(is.null(al_full[[i]]$sit$supal))
    }
})

test_that("align with external ref returns only input spectra", {
    skip_if_speaq_deps_missing()
    # When ref is supplied the reference is prepended internally but must be
    # stripped from the result, so the caller gets exactly length(x) spectra.
    al_auto <- align(decons, verbose = FALSE)
    al_ref  <- align(decons, verbose = FALSE, ref = decons[[1]])
    expect_equal(length(al_ref), length(decons))
    expect_equal(names(al_ref), names(decons))
})

test_that("align raises error for spectra with mismatched data-point counts", {
    skip_if_speaq_deps_missing()
    short <- simulate_spectrum(ndp = 128, npk = 2)
    long  <- simulate_spectrum(ndp = 256, npk = 2)
    d_short <- deconvolute(short, sfr = c(Inf, -Inf), force = TRUE, verbose = FALSE)
    d_long  <- deconvolute(long,  sfr = c(Inf, -Inf), force = TRUE, verbose = FALSE)
    mixed <- structure(list(d_short, d_long), class = "decons2")
    expect_error(align(mixed, verbose = FALSE))
})


skip_if_slow_tests_disabled()
