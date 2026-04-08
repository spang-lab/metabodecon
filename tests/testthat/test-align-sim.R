library(testthat)

skip("Manual Test for r4 only")

aki <- read_aki_data()
decons <- deconvolute(aki$spectra, sfr = c(11, -2), verbose = FALSE, nworkers=30)
test_that("built-in backend matches speaq backend", {
    system.time(al_speaq <- align(decons, verbose=TRUE, method=1, full=FALSE))  # 177 seconds
    system.time(al_builtin <- align(decons, verbose=TRUE, method=2, full=FALSE)) # 95 seconds
    system.time(al_fast <- align(decons, verbose=TRUE, method=3, full=FALSE)) # 6.4 seconds
    expect_equal(abi, asp)
})
