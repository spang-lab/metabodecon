library(testthat)

test_that("align works", {

    decons <- deconvolute(sim, sfr = c(3.55, 3.35))

    # Do the alignment
    aligns <- align(decons)

    # Then make all changes to the decons object that we expect [align()] to
    # make. At the end the objects should be equal. This way we can ensure that
    # [align()] adds exactly the fields we expect it to add (and no others).
    # However, this does NOT check, whether the added values contain sensible
    # values. For now we can use the downstream unit tests (e.g. [dohCluster()],
    # [combine_peaks()]), visual inspection and some sanity checks (see below)
    # for validation of the resulting values, but at some point it would be good
    # to setup a simulated dataset with known shifts and test against those.
    for (i in seq_along(aligns)) {
        decons[[i]]$sit$al <- aligns[[i]]$sit$al
        decons[[i]]$sit$supal <- aligns[[i]]$sit$supal
        decons[[i]]$lcpar$x0_al <- aligns[[i]]$lcpar$x0_al
        class(decons[[i]]) <- "align"
    }
    class(decons) <- "aligns"
    expect_equal(object = aligns, expected = decons)

    # Sanity checks
    align <- aligns[[1]]
    cs <- align$cs
    al <- align$sit$al
    supal <- align$sit$supal
    A <- align$lcpar$A
    x0_al <- align$lcpar$x0_al
    lamdba <- align$lcpar$lambda

    expect_equal(al[al > 0], A * pi)
    expect_equal(supal, lorentz_sup(cs, x0_al, A, lamdba))

})
