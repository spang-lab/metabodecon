# Inputs for all Tests
spectra <- sim[1:2]
sfr <- c(3.55, 3.35)

# Expected output if defaults are accepted
expected <- list(sim_01 = c(3.55, 3.35), sim_02 = c(3.55, 3.35))

test_that("get_sfr works", {
    # Inputs: ask=FALSE, adjno=0
    # Expectation: Defaults are returned
    obj0 <- get_sfr(spectra, sfr, ask = FALSE, adjno = 0)
    obj1 <- get_sfr(spectra, sfr, ask = FALSE, adjno = 1)
    expect_equal(obj0, expected)
    expect_equal(obj1, expected)
})
