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
})
