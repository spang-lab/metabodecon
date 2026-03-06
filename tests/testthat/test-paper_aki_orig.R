test_that("get_test_ids creates stratified folds", {
    y <- factor(rep(c("Control", "AKI"), each = 30))
    ids <- get_test_ids(seq_len(length(y)), nfolds = 3, seed = 1, y = y)

    expect_equal(length(ids), 3)
    expect_equal(length(sort(unlist(ids))), length(y))
    expect_equal(sort(unlist(ids)), seq_len(length(y)))

    for (i in seq_along(ids)) {
        yi <- y[ids[[i]]]
        expect_true(any(yi == "Control"))
        expect_true(any(yi == "AKI"))
    }
})

test_that("binary label coercion supports arbitrary factor levels", {
    y <- factor(c("Control", "AKI", "Control", "AKI"),
      levels = c("Control", "AKI"))
    expect_equal(as_binary01(y), c(0L, 1L, 0L, 1L))
})

test_that("binary label coercion rejects non-binary factors", {
    y <- factor(c("A", "B", "C"))
    expect_error(as_binary01(y), "exactly 2")
})

test_that("bin_spectra bins toy integer example correctly", {
    sp <- list(list(
        cs = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
        si = c(2, 4, 1, 3, 5, 2, 6, 1, 4, 2)
    ))
    reg <- matrix(c(10, 0), ncol = 2)

    x <- bin_spectra(sp, regions = reg, binwidth = 2)
    expect_equal(dim(x), c(1L, 5L))
    expect_equal(as.numeric(x[1, ]), c(6, 4, 7, 7, 6))
    expect_equal(attr(x, "bin_centers"), c(9, 7, 5, 3, 1))
})

test_that("plot_prob_scatter runs on toy inputs", {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    expect_silent(plot_prob_scatter(c(0.1, 0.8, 0.4), c(0, 1, 0), "prob"))
    grDevices::dev.off()
})

test_that("single svm fit is below runtime threshold", {

    testthat::skip_if_not_installed("e1071")

    sys <- Sys.info()
    is_dev_machine <- (sys[["user"]] == "tobi" && sys[["sysname"]] == "Darwin")
    runtime_factor <- if (is_dev_machine) 1 else 100

    set.seed(1)
    n <- 120
    p <- 12
    x <- matrix(rnorm(n * p), nrow = n)
    y <- factor(sample(c(0, 1), n, replace = TRUE), levels = c(0, 1))

    t <- system.time({
        fit <- e1071::svm(
            x = x,
            y = y,
            type = "C-classification",
            kernel = "radial",
            cost = 0.5,
            gamma = 2^-10,
            probability = FALSE
        )
        expect_s3_class(fit, "svm")
    })[["elapsed"]]

    expect_lte(t, 1 * runtime_factor)
})

test_that("estimate_svm_performance returns pooled predictions frame", {

    testthat::skip_if_not_installed("e1071")

    set.seed(1)
    n <- 24
    p <- 10
    X <- matrix(rnorm(n * p), nrow = n)
    y <- factor(rep(c("Control", "AKI"), each = n / 2),
      levels = c("Control", "AKI"))
    te <- get_test_ids(seq_len(n), nfolds = 4, seed = 1, y = y)

    Y <- expect_silent(estimate_svm_performance(
        X = X,
        y = y,
        test_ids = te,
        costs = 2^(0:1),
        gammas = 2^(-3:-2),
        nfeats = 2:3,
        inner_nfolds = 3,
        verbose = FALSE,
        ncores = 1
    ))

    expect_s3_class(Y, "data.frame")
    expect_equal(nrow(Y), n)
    expect_setequal(colnames(Y), c("sample", "fold", "cost", "gamma",
        "nfeat", "true", "prob", "pred"))
})


