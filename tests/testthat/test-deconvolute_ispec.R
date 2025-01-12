
sap2 <- test_that("deconvolute_ispec works", {

    args0 <- args1 <- args2 <- list(
        ispec   = as_ispec(sap2),
        nfit    = 3,
        sfr     = c(3.2, -3.2),
        smopts  = c(1, 3),
        delta   = 3,
        bwc     = 0,
        verbose = FALSE
    )
    args1$bwc <- 1
    args2$bwc <- 2

    # Don't use a loop, so we can see which test failed in case of errors.
    idecon0 <- do.call(deconvolute_ispec, args0)
    idecon1 <- do.call(deconvolute_ispec, args1)
    idecon2 <- do.call(deconvolute_ispec, args2)

    expect_identical(object = names(idecon0), expected = idecon_members)
    expect_identical(object = names(idecon1), expected = idecon_members)
    expect_identical(object = names(idecon2), expected = idecon_members)

    expect_identical(object = class(idecon0), expected = "idecon")
    expect_identical(object = class(idecon1), expected = "idecon")
    expect_identical(object = class(idecon2), expected = "idecon")

    obj0 <- calc_prarp(x = idecon0, truepar = sap2$meta$simpar)
    obj1 <- calc_prarp(x = idecon1, truepar = sap2$meta$simpar)
    obj2 <- calc_prarp(x = idecon2, truepar = sap2$meta$simpar)

    expect_true(obj0$prarpx >= 0.507) # MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
    expect_true(obj1$prarpx >= 0.961)
    expect_true(obj2$prarpx >= 0.961)

    expect_true(obj2$prarpx >= obj1$prarpx)
    expect_true(obj1$prarpx >= obj0$prarpx)
})
