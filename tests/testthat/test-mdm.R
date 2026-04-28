# Minimal smoke tests for fit_mdm(), cv_mdm() and benchmark_mdm().
# Uses tiny simulated spectra (512 points, 3 peaks, 12 samples) to keep
# runtime short enough for CI.

set.seed(1)
n <- 32
npk <- 3
cs <- seq(from = 3.6, length.out = 512, by = -0.0006)
x0 <- sort(runif(npk, 3.42, 3.56))
A <- runif(npk, 8, 14) * 1e3
lam <- runif(npk, 0.9, 1.3) / 1e3
y <- factor(rep(c("A", "B"), each = n / 2))
sp <- vector("list", n)
for (i in seq_len(n)) {
    xi <- x0 + rnorm(npk, sd = 0.0003)
    Ai <- A * runif(npk, 0.8, 1.2)
    li <- lam * runif(npk, 0.9, 1.1)
    Ai[1] <- Ai[1] * (if (y[i] == "A") 1.3 else 0.7)
    sp[[i]] <- simulate_spectrum(
        name = sprintf("s_%02d", i), cs = cs,
        x0 = sort(xi), A = Ai, lambda = li,
        noise = rnorm(length(cs), sd = 500)
    )
}
class(sp) <- "spectra"

testthat::test_that("fit_mdm returns mdm object", {
    m <- fit_mdm(
        sp, y,
        nfit = 1, smit = 1, smws = 3, delta = 8, npmax = 0,
        maxShift = 50, maxCombine = 20, nfolds = 3,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(!is.null(m$model) || !is.null(m$meta))
})

testthat::test_that("cv_mdm returns mdm with pgrid", {
    pg <- data.frame(
        smit = 1, smws = 3, delta = c(6, 8),
        nfit = 1, npmax = 0,
        maxShift = 50, maxCombine = 20,
        stringsAsFactors = FALSE
    )
    pg$acc <- NA_real_
    pg$auc <- NA_real_
    m <- cv_mdm(
        sp, y,
        pgrid = pg, nfolds = 3,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_true(is.data.frame(m$pgrid))
    testthat::expect_equal(nrow(m$pgrid), 2)
})

testthat::test_that("benchmark_mdm returns predictions", {

    pg <- data.frame(
        smit = 1, smws = 3, delta = 8,
        nfit = 1, npmax = 0,
        maxShift = 50, maxCombine = 20,
        stringsAsFactors = FALSE
    )
    pg$acc <- NA_real_
    pg$auc <- NA_real_
    bm <- benchmark_mdm(
        sp, y,
        pgrid = pg, nfo = 4, nfl = 3,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_true(is.data.frame(bm$predictions))
    testthat::expect_true("true" %in% names(bm$predictions))
    testthat::expect_equal(nrow(bm$predictions), length(y))
})

testthat::test_that("fit_mdm returns mdm with NULL model when all spectra have zero peaks", {
    empty_lcpar <- data.frame(x0 = numeric(0), lambda = numeric(0), A = numeric(0))
    make_empty_decons <- function() {
        ds <- lapply(sprintf("s_%02d", seq_len(n)), function(nm) {
            structure(list(lcpar = empty_lcpar, meta = list(name = nm)), class = "decon2")
        })
        structure(setNames(ds, sprintf("s_%02d", seq_len(n))), class = "decons2")
    }
    with_mocked_bindings(
        deconvolute = function(...) make_empty_decons(),
        code = {
            m <- fit_mdm(
                sp, y,
                nfit = 1, smit = 1, smws = 3, delta = 8, npmax = 0,
                maxShift = 50, maxCombine = 20, nfolds = 3,
                use_rust = 0.5, nworkers = 1, verbosity = 0
            )
        }
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_null(m$model)
})

testthat::test_that("cv_mdm respects the ignore parameter", {
    pg <- data.frame(
        smit = 1, smws = 3, delta = 8,
        nfit = 1, npmax = 0,
        maxShift = 50, maxCombine = 20,
        stringsAsFactors = FALSE
    )
    pg$acc <- NA_real_
    pg$auc <- NA_real_
    m <- cv_mdm(
        sp, y,
        pgrid = pg, nfolds = 3,
        ignore = 1:4,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    testthat::expect_s3_class(m, "mdm")
    testthat::expect_equal(nrow(m$pgrid), 1)
    testthat::expect_false(is.na(m$pgrid$auc[1]))
})

testthat::test_that("fit_mdm caches deconvolution result in RAM when hash attribute is set", {
    options(metabodecon.fit_mdm.cache = NULL)
    sp_hashed <- sp
    attr(sp_hashed, "hash") <- rlang::hash(sp)
    fit_mdm(
        sp_hashed, y,
        nfit = 1, smit = 1, smws = 3, delta = 8, npmax = 0,
        maxShift = 50, maxCombine = 20, nfolds = 3,
        use_rust = 0.5, nworkers = 1, verbosity = 0
    )
    cache <- getOption("metabodecon.fit_mdm.cache")
    testthat::expect_false(is.null(cache))
    testthat::expect_true("decons" %in% names(cache))
    testthat::expect_true("aligns" %in% names(cache))
    options(metabodecon.fit_mdm.cache = NULL)
})
