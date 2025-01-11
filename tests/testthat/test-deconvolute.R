check_prarp <- test_that("deconvolute_ispecs works", {
    xx <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))

    obj1 <- calc_prarp(xx[[1]], truepar = sim[[1]]$meta$simpar)
    expect_true(obj$prarp > 0.84)
    expect_true(obj$prarpx >= 0.78)

    truepars <- lapply(sim[1:2], function(x) x$meta$simpar)
    prarp_objs <- mapply(calc_prarp, xx, truepars, SIMPLIFY = FALSE)
    prarpxs <- lapply(prarp_objs, function(obj) obj$prarpx)
    expect_true(all(prarpxs > 0.7))
    expect_equal(class(xx), "decons2")
    expect_equal(length(xx), 2)


    if (identical(environment(), .GlobalEnv)) {
        plot_spectrum(xx[[1]])
        plot_spectrum(xx[[2]])
    }
})