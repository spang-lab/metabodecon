library(testthat)

# IMPORTANT: we don't test for a good PRARP and/or different values of BWC here,
# as these checks are already done in `test-deconvolute_spectrum_r` (singular). Here,
# we only want to assert that the "glue code" in `test-deconvolute_spectra_r`
# (plural), orchestrating multiple invocations of `deconvolute_spectrum_r`, works as
# expected. We do this by checking that the returned objects have the correct
# types and the same PRARPs as returned `test-deconvolute_spectrum_r`.

# Call #####
idecons <- deconvolute_spectra_r(sim[1:2], nfit = 3, sfr = c(3.55, 3.35), verbose = FALSE, bwc = 0)
idecon1 <- deconvolute_spectrum_r(sim[[1]], nfit = 3, sfr = c(3.55, 3.35), verbose = FALSE, bwc = 0)
idecon2 <- deconvolute_spectrum_r(sim[[2]], nfit = 3, sfr = c(3.55, 3.35), verbose = FALSE, bwc = 0)

# Checks #####

types <- test_that("Types are ok", {
    expect_identical(class(idecons), "idecons")
    expect_identical(class(idecons[[1]]), "idecon")
    expect_identical(class(idecons[[2]]), "idecon")
    expect_identical(names(idecons), c("sim_01", "sim_02"))
    expect_identical(names(idecons[[1]]), idecon_members)
    expect_identical(names(idecons[[2]]), idecon_members)
})

prarps <- test_that("PRARPs are good", {
    truepars <- lapply(sim[1:2], function(x) x$meta$simpar)
    prarpxs <- c(
        calc_prarp(idecons[[1]], truepar = sim[[1]]$meta$simpar)$prarpx,
        calc_prarp(idecons[[2]], truepar = sim[[2]]$meta$simpar)$prarpx
    )
    prarpx1 <- calc_prarp(idecon1, truepar = sim[[1]]$meta$simpar)$prarpx
    prarpx2 <- calc_prarp(idecon2, truepar = sim[[2]]$meta$simpar)$prarpx
    expect_equal(prarpxs[1], prarpx1)
    expect_equal(prarpxs[2], prarpx2)
})
