skip_if_slow_tests_disabled()

test_that("cache_example_datasets(persistent = FALSE) works with existing temporary datadir", {
    y <- evalwith(datadir_persistent = "empty", datadir_temp = "empty", message = "captured", expr = {
        zip_returned <- cache_example_datasets(persistent = FALSE)
        zip_temp <- zip_temp()
        zip_persistent <- zip_persistent()
    })
    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    msg <- paste(y$message, collapse = "\n")
    expect_match(msg, paste0("Downloading ", xds$url, " as "))
    z1 <- normalizePath(zip_returned, "/", mustWork = FALSE)
    z2 <- normalizePath(zip_temp, "/", mustWork = FALSE)
    expect_equal(z1, z2)
})

test_that("cache_example_datasets(persistent = FALSE) works with empty datadirs", {
    y <- evalwith(datadir_persistent = "empty", datadir_temp = "filled", message = "captured", expr = {
        zip_returned <- cache_example_datasets(persistent = FALSE)
        zip_temp <- zip_temp()
        zip_persistent <- zip_persistent()
    })
    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    z1 <- normalizePath(zip_returned, "/", mustWork = FALSE)
    z2 <- normalizePath(zip_temp, "/", mustWork = FALSE)
    expect_equal(z1, z2)
})
