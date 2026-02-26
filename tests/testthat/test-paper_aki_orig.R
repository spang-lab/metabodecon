test_that("get_test_ids creates stratified folds", {
    y <- rep(c(0, 1), each = 30)
    ids <- get_test_ids(nfolds = 3, nsamples = length(y), seed = 1, y = y)

    expect_equal(length(ids), 3)
    expect_equal(length(sort(unlist(ids))), length(y))
    expect_equal(sort(unlist(ids)), seq_len(length(y)))

    for (i in seq_along(ids)) {
        yi <- y[ids[[i]]]
        expect_true(any(yi == 0))
        expect_true(any(yi == 1))
    }
})

test_that("paper_aki cache helpers store and retrieve values", {
    cache <- paper_aki_cache_new()
    key <- paper_aki_cache_key("grid", 1:3, c(0.2, 0.5, 0.8), 2^c(-11, -10, -9))

    expect_null(paper_aki_cache_get(cache, key))
    x <- data.frame(a = 1:3, b = letters[1:3])
    paper_aki_cache_set(cache, key, x)

    y <- paper_aki_cache_get(cache, key)
    expect_equal(y, x)
})

test_that("paper_aki cache size and clear work", {
    cache <- paper_aki_cache_new(tempfile("aki-cache-"))
    key <- paper_aki_cache_key("demo", 1:3)
    val <- data.frame(a = 1:3)
    paper_aki_cache_set(cache, key, val)

    sz <- paper_aki_cache_size(cache)
    expect_true(sz$files >= 1)
    expect_true(sz$bytes > 0)

    nrm <- paper_aki_cache_clear(cache)
    expect_true(nrm >= 1)

    sz2 <- paper_aki_cache_size(cache)
    expect_equal(sz2$files, 0L)
    expect_equal(sz2$bytes, 0)
})

test_that("paper_aki_get_cache_dir prefers persistent cache", {
    pdir <- file.path(datadir_persistent(), "cache")
    dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    old <- getOption("metabodecon.aki_cache")
    options(metabodecon.aki_cache = tempfile("aki-cache-opt-"))
    on.exit(options(metabodecon.aki_cache = old), add = TRUE)

    got <- paper_aki_get_cache_dir()
    expect_equal(got, pdir)
    expect_equal(getOption("metabodecon.aki_cache"), pdir)
})

test_that("train_model caches repeated fits", {
    testthat::skip_if_not_installed("e1071")

    cache <- paper_aki_get_cache_dir()
    X <- matrix(rnorm(80), nrow = 20)
    y <- rep(0:1, each = 10)
    payload <- list(X = X, y = y)
    txt <- serialize(payload, connection = NULL, ascii = TRUE, version = 2)
    sig <- paper_aki_cache_hash(rawToChar(txt))
    key <- paper_aki_cache_key("train_model", sig, 0.5, 2^-10, 3, TRUE)
    file <- paper_aki_cache_file(cache, key)
    if (file.exists(file)) {
        file.remove(file)
    }

    n0 <- paper_aki_cache_size(cache)$files
    m1 <- train_model(X, y, cost = 0.5, gamma = 2^-10, nfeat = 3)
    n1 <- paper_aki_cache_size(cache)$files
    m2 <- train_model(X, y, cost = 0.5, gamma = 2^-10, nfeat = 3)
    n2 <- paper_aki_cache_size(cache)$files

    expect_true(file.exists(file))
    expect_true(n1 >= n0)
    expect_equal(n2, n1)
    expect_equal(m2$idx, m1$idx)
})

test_that("plot helpers run on toy inputs", {
    ids <- list(c(1, 4), c(2, 5), c(3, 6))
    grid <- expand.grid(cost = c(0.2, 0.5, 0.8),
                        gamma = 2^c(-11, -10, -9),
                        nfeat = 2:3)
    grid$auc <- seq(0.6, 0.9, length.out = nrow(grid))

    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    expect_silent(plot_folds(ids, nsamples = 6, main = "outer"))
    expect_silent(plot_auc_grid(grid, nfeat = 2, main = "grid"))
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

test_that("grid apply gives same result on 1 vs 2 cores", {
    testthat::skip_if_not_installed("e1071")

    set.seed(1)
    n <- 120
    p <- 12
    x <- matrix(rnorm(n * p), nrow = n)
    y <- factor(sample(c(0, 1), n, replace = TRUE), levels = c(0, 1))
    grid <- expand.grid(
        cost = c(0.2, 0.5, 0.8),
        gamma = 2^c(-11, -10, -9),
        nfeat = 2:3
    )

    run_one <- function(g, i) {
        idx <- seq_len(min(g$nfeat[[1]], ncol(x)))
        fit <- e1071::svm(
            x = x[, idx, drop = FALSE],
            y = y,
            type = "C-classification",
            kernel = "radial",
            cost = g$cost[[1]],
            gamma = g$gamma[[1]],
            probability = FALSE
        )
        c(i = i, nsv = length(fit$index))
    }

    out1 <- paper_aki_grid_apply(grid, run_one, ncores = 1)
    out2 <- paper_aki_grid_apply(grid, run_one, ncores = 2)

    out1 <- do.call(rbind, out1)
    out2 <- do.call(rbind, out2)

    expect_equal(out1, out2)
    expect_equal(attr(out1, "dim"), attr(out2, "dim"))
})
