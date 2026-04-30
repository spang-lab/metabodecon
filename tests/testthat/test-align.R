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
    smit = 1, smws = 3,
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
        decons_copy[[i]]$sit$supal    <- aligns[[i]]$sit$supal
        decons_copy[[i]]$lcpar$x0al  <- aligns[[i]]$lcpar$x0al
        decons_copy[[i]]$lcpar$pcial  <- aligns[[i]]$lcpar$pcial
        class(decons_copy[[i]]) <- c("align", "decon2", "spectrum")
    }
    class(decons_copy) <- c("aligns", "decons2", "spectra")
    expect_equal(object = aligns, expected = decons_copy)

    # Check that the alignment worked, our expectations are:
    # 1. x0al     is shifted roughly 0.3 to the right compared to x0
    # 2. pcial     indexes cs at the aligned peak centers (cs[pcial] == x0al)
    # 3. sit$supal is the superposition of the aligned Lorentz curves
    x0 <- aligns$sap_01_shifted$lcpar$x0
    x0al <- aligns$sap_01_shifted$lcpar$x0al
    shifts <- x0 - x0al
    expect_true(all(shifts > 0.2 & shifts < 0.4))

    cs <- aligns$sap_01_shifted$cs
    pcial <- aligns$sap_01_shifted$lcpar$pcial
    expect_equal(cs[pcial], x0al)

    A <- aligns$sap_01_shifted$lcpar$A
    supal <- aligns$sap_01_shifted$sit$supal
    lambda <- aligns$sap_01_shifted$lcpar$lambda
    expect_equal(supal, lorentz_sup(cs, x0al, A, lambda))
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
        smit = 1, smws = 3,
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

    al_builtin <- align_decons(decons, verbose = FALSE, use_speaq = FALSE)
    al_speaq   <- align_decons(decons, verbose = FALSE, use_speaq = TRUE)
    expect_equal(al_builtin, al_speaq)
})

test_that("align with full=FALSE omits supal", {
    skip_if_speaq_deps_missing()
    decons_h <- decons
    attr(decons_h, "hash") <- rlang::hash(decons_h)
    al_nofull <- align_decons(decons_h, verbose = FALSE, full = FALSE)
    for (i in seq_along(al_nofull)) {
        expect_null(al_nofull[[i]]$sit$supal)
    }
    # full=TRUE (default) should include supal
    al_full <- align_decons(decons_h, verbose = FALSE, full = TRUE)
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
    # Build two decon2 objects with different cs lengths directly to bypass
    # deconvolute() (which would itself fail on a degenerate sfr).
    short <- decons[[1]]
    long  <- decons[[1]]
    long$cs <- c(long$cs, long$cs[length(long$cs)] - 0.001)
    long$si <- c(long$si, 0)
    long$sit <- rbind(long$sit, long$sit[nrow(long$sit), ])
    mixed <- structure(list(short, long), class = c("decons2", "spectra"))
    expect_error(align(mixed, verbose = FALSE))
})


skip_if_slow_tests_disabled()
