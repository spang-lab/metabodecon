# Inputs for all Tests
spectra <- sim[1:2]
sfr <- c(3.55, 3.35)

# Expected output if defaults are accepted
expected <- list(sim_01 = c(3.55, 3.35), sim_02 = c(3.55, 3.35))

test_that("get_sfr works", {
    obj <- get_sfr(spectra, sfr)
    expect_equal(obj, expected)
})
