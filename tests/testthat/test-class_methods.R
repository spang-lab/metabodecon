test_that("format and summary work for spectrum and spectra", {
    s1 <- sim[[1]]
    ss <- sim[1:2]

    expect_true(is.character(format(s1)))
    expect_match(format(s1), "spectrum object")

    sum1 <- summary(s1)
    expect_true(is.list(sum1))
    expect_equal(sum1$n_dp, length(s1$cs))

    expect_true(is.character(format(ss)))
    expect_match(format(ss), "spectra object")

    sum2 <- summary(ss)
    expect_true(is.data.frame(sum2))
    expect_equal(nrow(sum2), 2)
})

test_that("format and summary work for decon2 and align", {
    d2 <- deconvolute(sim[[1]], sfr = c(3.55, 3.35), verbose = FALSE)
    a1 <- d2
    class(a1) <- c("align", "decon2", "spectrum")

    expect_match(format(d2), "decon2 object")
    expect_match(format(a1), "align object")

    s2 <- summary(d2)
    sa <- summary(a1)
    expect_true(is.list(s2))
    expect_true(is.list(sa))
    expect_equal(s2$n_peaks, length(d2$lcpar$A))
    expect_equal(sa$n_peaks, length(a1$lcpar$A))
})

test_that("c.spectra combines spectra subsets", {
    x <- c(head(sim, 2), tail(sim, 2))
    expect_true(is_spectra(x))
    expect_equal(length(x), 4)
    expect_equal(names(x), c("sim_01", "sim_02", "sim_15", "sim_16"))
})

test_that("c.spectra combines spectrum and spectra", {
    x <- c(sim[[1]], sim[2:3])
    expect_true(is_spectra(x))
    expect_equal(length(x), 3)
    expect_equal(names(x), c("sim_01", "sim_02", "sim_03"))
})

test_that("c.spectra rejects unsupported inputs", {
    expect_error(c(sim, 123), "spectrum or spectra")
})
