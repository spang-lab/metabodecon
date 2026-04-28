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

test_that("c methods work for decon1/decon2/align collections", {
    d1 <- generate_lorentz_curves_sim(sim[[1]])
    d2 <- as_decon2(d1)
    a1 <- d2
    class(a1) <- "align"

    cc1 <- c(d1, d1)
    expect_true(is_decons1(cc1))
    expect_equal(length(cc1), 2)

    cc2 <- c(d2, d2)
    expect_true(is_decons2(cc2))
    expect_equal(length(cc2), 2)

    cca <- c(a1, a1)
    expect_true(is_aligns(cca))
    expect_equal(length(cca), 2)
})

test_that("format and summary work for decon2 and align", {
    d1 <- generate_lorentz_curves_sim(sim[[1]])
    d2 <- as_decon2(d1)
    a1 <- d2
    class(a1) <- "align"

    expect_match(format(d2), "decon2 object")
    expect_match(format(a1), "align object")

    s2 <- summary(d2)
    sa <- summary(a1)
    expect_true(is.list(s2))
    expect_true(is.list(sa))
    expect_equal(s2$n_peaks, length(d2$lcpar$A))
    expect_equal(sa$n_peaks, length(a1$lcpar$A))
})

test_that("private c methods combine rdecon objects", {
    rd1 <- structure(list(meta = list(name = "rd1")), class = "rdecon")
    rd2 <- structure(list(meta = list(name = "rd2")), class = "rdecon")

    xrd <- c(rd1, rd2)

    expect_true(inherits(xrd, "rdecons"))
    expect_equal(length(xrd), 2)
})

test_that("private format methods return readable labels for rdecon", {
    rd1 <- structure(list(meta = list(name = "rd1")), class = "rdecon")

    expect_match(format(rd1), "rdecon object")
})

test_that("private summary methods return expected structures for rdecon", {
    rd1 <- structure(list(meta = list(name = "rd1")), class = "rdecon")

    sr <- summary(rd1)

    expect_true(is.list(sr))
    expect_equal(sr$name, "rd1")
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
