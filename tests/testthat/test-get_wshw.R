# Inputs for all Tests
x <- sim[1:2]
wshw <- 0.01

# Expected output if defaults are accepted
xpct <- list(sim_01 = wshw, sim_02 = wshw)

test_that("get_wshw works", {
    # Inputs: ask=FALSE, adjno=0
    # Expectation: Defaults are returned
    obj0 <- get_wshw(x, wshw, ask = FALSE, adjno = 0)
    obj1 <- get_wshw(x, wshw, ask = FALSE, adjno = 1)
    expect_equal(obj0, xpct)
    expect_equal(obj1, xpct)

    # Inputs: ask=TRUE, adjno=1
    # Expectation: user is asked to confirm/update values for the first spectrum
    evalwith(
        answers = c(WSHWok = "y"),
        plot = "captured",
        message = "captured",
        obj0 <- get_wshw(x, wshw, ask = TRUE, adjno = 1)
    )
    evalwith(
        answers = c(WSHWok = "n", WSHWval = 0.02, Invalid = "bb", WSHWok = "y"),
        plot = "captured",
        message = "captured",
        obj1 <- get_wshw(x, wshw, ask = TRUE, adjno = 1)
    )
    expect_equal(obj0, xpct)
    expect_equal(obj1, list(sim_01 = 0.02, sim_02 = 0.02))

    # Inputs: ask=TRUE, adjno=0
    # Expectation: user is asked to confirm/update values for all spectra
    evalwith(
        answers = c(WSHWok = "y", WSHWok = "y"),
        plot = "captured",
        message = "captured",
        obj <- get_wshw(x, wshw, ask = TRUE, adjno = 0)
    )
    expect_equal(obj, xpct)
})
