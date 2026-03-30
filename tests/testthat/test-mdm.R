# Tests for the planned fit_mdm(X, y, model, normalisation, ...) API.
# Currently skipped because fit_mdm has a different signature.

testthat::test_that("fit_mdm supports svm", {
  testthat::skip("fit_mdm(X, y, model=...) API not yet implemented")
  testthat::skip_if_not_installed("e1071")
  set.seed(1)
  X <- matrix(rnorm(30 * 20), nrow = 30)
  y <- factor(rep(c("Control", "AKI"), each = 15))

  m <- fit_mdm(X, y, model = "svm", normalisation = "quantile", nfeat = 5)

  testthat::expect_s3_class(m, "mdm")
  testthat::expect_true(is.list(m))
  testthat::expect_true(all(c("model", "normalisation", "ref") %in% names(m)))

  p <- predict(m, X, type = "prob")
  c <- predict(m, X, type = "class")
  testthat::expect_length(p, nrow(X))
  testthat::expect_length(c, nrow(X))
  testthat::expect_true(all(c %in% c(0L, 1L)))
})

testthat::test_that("fit_mdm supports lasso", {
  testthat::skip("fit_mdm(X, y, model=...) API not yet implemented")
  testthat::skip_if_not_installed("glmnet")
  set.seed(1)
  X <- matrix(rnorm(40 * 25), nrow = 40)
  y <- factor(rep(c("Control", "AKI"), each = 20))

  m <- fit_mdm(X, y, model = "lasso", normalisation = "median")

  testthat::expect_s3_class(m, "mdm")
  p <- predict(m, X, type = "prob")
  testthat::expect_length(p, nrow(X))

  cf <- coef(m)
  testthat::expect_true(nrow(as.matrix(cf)) >= 1)
})

testthat::test_that("mdm methods run", {
  testthat::skip("fit_mdm(X, y, model=...) API not yet implemented")
  testthat::skip_if_not_installed("e1071")
  X <- matrix(rnorm(20 * 10), nrow = 20)
  y <- factor(rep(c("Control", "AKI"), each = 10))
  m <- fit_mdm(X, y, model = "svm", normalisation = "sum", nfeat = 3)

  testthat::expect_invisible(print(m))
  s <- summary(m)
  testthat::expect_s3_class(s, "summary.mdm")

  testthat::expect_invisible(plot(m))
})

testthat::test_that("extension benchmark returns expected columns", {
  testthat::skip("benchmark_extension_table not yet implemented")
  testthat::skip_if_not_installed("e1071")
  set.seed(3)
  X <- matrix(rnorm(30 * 15), nrow = 30)
  y <- factor(rep(c("Control", "AKI"), each = 15))

  tab <- benchmark_extension_table(
    X,
    y,
    models = "svm",
    norms = c("none", "sum"),
    k = 3,
    costs = 0.5,
    gammas = 2^-10,
    nfeats = 2
  )

  testthat::expect_true(all(c("model", "normalisation") %in% names(tab)))
  testthat::expect_true(all(c("auc_cv_norm", "auc_global_norm") %in% names(tab)))
  testthat::expect_equal(nrow(tab), 2)
})

testthat::test_that("extension benchmark uses cache", {
  testthat::skip("benchmark_extension_table not yet implemented")
  testthat::skip_if_not_installed("e1071")
  X <- matrix(rnorm(30 * 15), nrow = 30)
  y <- factor(rep(c("Control", "AKI"), each = 15))
  cache <- paper_aki_get_cache_dir()

  key <- paper_aki_cache_key(
    "benchmark_extension_table",
    mdm_data_signature(X, y),
    "svm",
    c("none", "sum"),
    3,
    0.5,
    2^-10,
    2
  )
  file <- paper_aki_cache_file(cache, key)
  if (file.exists(file)) {
    file.remove(file)
  }

  n0 <- paper_aki_cache_size(cache)$files
  t1 <- benchmark_extension_table(
    X,
    y,
    models = "svm",
    norms = c("none", "sum"),
    k = 3,
    costs = 0.5,
    gammas = 2^-10,
    nfeats = 2
  )
  n1 <- paper_aki_cache_size(cache)$files
  t2 <- benchmark_extension_table(
    X,
    y,
    models = "svm",
    norms = c("none", "sum"),
    k = 3,
    costs = 0.5,
    gammas = 2^-10,
    nfeats = 2
  )
  n2 <- paper_aki_cache_size(cache)$files

  testthat::expect_true(file.exists(file))
  testthat::expect_true(n1 >= n0)
  testthat::expect_equal(n2, n1)
  testthat::expect_equal(t2, t1)
})

testthat::test_that("fit_mdm rejects non-binary labels", {
  testthat::skip("fit_mdm(X, y, model=...) API not yet implemented")
  X <- matrix(rnorm(30), nrow = 10)
  y <- factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A"))
  testthat::expect_error(
    fit_mdm(X, y, model = "svm", normalisation = "none", nfeat = 2),
    "exactly 2"
  )
})
