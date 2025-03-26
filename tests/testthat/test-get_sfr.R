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

    # Inputs: ask=TRUE, adjno=1
    # Expectation: user is asked to confirm/update values for the first spectrum
    evalwith(
        answers = c(SFRok = "y"),
        message = "captured",
        obj0 <- get_sfr(spectra, sfr, ask = TRUE, adjno = 1)
    )
    evalwith(
        answers = c(SFRok = "n", SFRleft = 3.52, Invalid = "bb", SFRright = 3.38, SFRok = "y"),
        message = "captured",
        obj1 <- get_sfr(spectra, sfr, ask = TRUE, adjno = 1)
    )
    expect_equal(obj0, expected)
    expect_equal(obj1, list(sim_01 = c(3.52, 3.38), sim_02 = c(3.52, 3.38)))

    # Inputs: ask=TRUE, adjno=0
    # Expectation: user is asked to confirm/update values for all spectra
    evalwith(
        answers = c(SFRok = "y", SFRok = "y"),
        message = "captured",
        obj1 <- get_sfr(spectra, sfr, ask = TRUE, adjno = 0)
    )
    expect_equal(obj0, expected)
    expect_equal(obj1, expected)
})
