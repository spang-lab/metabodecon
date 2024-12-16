check_prarp <- test_that("deconvolute_ispecs works", {
    decons2 <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    truepars <- lapply(sim[1:2], function(x) x$meta$simpar)
    prarps <- mapply(calc_prarp, decons2, truepars)
    expect_true(all(prarps > 0.8))
    expect_equal(class(decons2), "decons2")
    expect_equal(length(decons2), 2)
    if (environment() %==% .GlobalEnv) {
        plot_spectrum(decons2[[1]])
        plot_spectrum(decons2[[2]])
    }
})