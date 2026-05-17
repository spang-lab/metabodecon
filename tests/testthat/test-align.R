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
    # 2. pcial    indexes cssh at the aligned peak centers (cssh[pcial] == x0al)
    # 3. sit$supal is the superposition of the aligned Lorentz curves on cssh
    x0 <- aligns$sap_01_shifted$lcpar$x0
    x0al <- aligns$sap_01_shifted$lcpar$x0al
    shifts <- x0 - x0al
    expect_true(all(shifts > 0.2 & shifts < 0.4))

    cssh <- aligns$sap_01_shifted$cssh
    pcial <- aligns$sap_01_shifted$lcpar$pcial
    expect_equal(cssh[pcial], x0al)

    A <- aligns$sap_01_shifted$lcpar$A
    supal <- aligns$sap_01_shifted$sit$supal
    lambda <- aligns$sap_01_shifted$lcpar$lambda
    expect_equal(supal, lorentz_sup(cssh, x0al, A, lambda))
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

    al_builtin <- clupa(decons, verbose = FALSE, use_speaq = FALSE)
    al_speaq   <- clupa(decons, verbose = FALSE, use_speaq = TRUE)
    expect_equal(al_builtin, al_speaq)
})

test_that("built-in backend matches speaq backend on sim2 across maxShift", {

    # Verifies that the underlying CluPA / FFT shift search is bit-equivalent
    # to speaq for every positive maxShift used in the supervised parameter
    # grid. maxShift = 0 is handled at the clupa() wrapper level (no shift;
    # x0al = x0) and so does not exercise speaq, hence it is not in the loop.

    skip_if_speaq_deps_missing()
    skip_if_slow_tests_disabled()

    # Two sim2 spectra are sufficient: the alignment is per-spectrum and
    # sim2's per-peak jitter is what drives non-trivial shifts here.
    s <- sim2[1:2]
    d <- deconvolute(
        s, sfr = NULL, smit = 2, smws = 3, delta = 1.6,
        nfit = 10, npmax = 0, verbose = FALSE
    )
    for (ms in c(3, 5, 10, 50)) {
        al_builtin <- clupa(
            d, maxShift = ms, verbose = FALSE, use_speaq = FALSE
        )
        al_speaq <- clupa(
            d, maxShift = ms, verbose = FALSE, use_speaq = TRUE
        )
        for (i in seq_along(al_builtin)) {
            expect_equal(
                al_builtin[[i]]$lcpar$pcial,
                al_speaq[[i]]$lcpar$pcial,
                info = sprintf("maxShift = %d, spectrum %d", ms, i)
            )
            expect_equal(
                al_builtin[[i]]$lcpar$x0al,
                al_speaq[[i]]$lcpar$x0al,
                info = sprintf("maxShift = %d, spectrum %d", ms, i)
            )
        }
    }
})

test_that("maxShift = 0 short-circuits to a no-shift alignment", {

    # clupa must treat maxShift = 0 as a no-op: x0al = x0, pcial set to the
    # column nearest x0, returned as an `aligns` object. Lets mdm.R include
    # 0 as a valid grid-search point alongside positive shifts.

    skip_if_speaq_deps_missing()

    al <- clupa(decons, maxShift = 0L, verbose = FALSE)
    expect_s3_class(al, "aligns")
    for (i in seq_along(al)) {
        # No shift: aligned center equals the raw fitted center.
        expect_equal(al[[i]]$lcpar$x0al, al[[i]]$lcpar$x0)
        # pcial is the nearest cssh grid column for each (off-grid) x0.
        cssh <- al[[i]]$cssh
        nc <- length(cssh)
        expected_pcial <- pmin(nc, pmax(
            1L, round(metabodecon:::convert_pos(al[[i]]$lcpar$x0, cssh, seq_len(nc)))
        ))
        expect_equal(al[[i]]$lcpar$pcial, expected_pcial)
    }
})

test_that("align() chains CluPA and RefPA when maxCombine > 0", {
    skip_if_speaq_deps_missing()
    a0  <- align(decons, maxShift = 50, maxCombine = 0,  verbose = FALSE)
    a20 <- align(decons, maxShift = 50, maxCombine = 20, verbose = FALSE)
    # Every non-NA snapped peak must sit on a reference column.
    ref <- find_ref(a0)
    pp <- ref$lcpar$pcial
    for (i in seq_along(a20)) {
        pcisn <- a20[[i]]$lcpar$pcisn
        ok <- !is.na(pcisn)
        expect_true(all(pcisn[ok] %in% pp))
    }
    # RefPA never drops rows; it only annotates with pcisn / x0sn.
    expect_equal(
        vapply(a20, function(s) nrow(s$lcpar), integer(1)),
        vapply(a0,  function(s) nrow(s$lcpar), integer(1))
    )
})

test_that("align with full=FALSE omits supal", {
    skip_if_speaq_deps_missing()
    decons_h <- decons
    attr(decons_h, "hash") <- rlang::hash(decons_h)
    al_nofull <- clupa(decons_h, verbose = FALSE, full = FALSE)
    for (i in seq_along(al_nofull)) {
        expect_null(al_nofull[[i]]$sit$supal)
    }
    # full=TRUE (default) should include supal
    al_full <- clupa(decons_h, verbose = FALSE, full = TRUE)
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

test_that("align raises error for spectra with mismatched cssh grids", {
    skip_if_speaq_deps_missing()
    # Build two decon2 objects that already carry different cssh grids.
    # ensure_cssh() must refuse to align them since cssh indices would
    # not be comparable across spectra.
    a <- decons[[1]]
    b <- decons[[1]]
    b$cssh <- b$cssh + 0.5
    mixed <- structure(list(a, b), class = c("decons2", "spectra"))
    expect_error(align(mixed, verbose = FALSE))
})


skip_if_slow_tests_disabled()
