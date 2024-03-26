test_that("cache_example_datasets(persistent = FALSE) works with existing temporary datadir", {
    y <- evalwith(datadir_persistent = "empty", datadir_temp = "empty", message = "captured", expr = {
        zip_returned <- cache_example_datasets(persistent = FALSE)
        zip_temp <- zip_temp()
        zip_persistent <- zip_persistent()
    })
    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    expect_equal(y$message, paste("Downloading", xds$url, "as", zip_temp))
    expect_equal(zip_returned, zip_temp)
})

skip_if_slow_tests_disabled()

test_that("cache_example_datasets(persistent = FALSE) works with empty datadirs", {
    y <- evalwith(datadir_persistent = "empty", datadir_temp = "filled", message = "captured", expr = {
        zip_returned <- cache_example_datasets(persistent = FALSE)
        zip_temp <- zip_temp()
        zip_persistent <- zip_persistent()
    })
    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    expect_equal(y$message, character())
    expect_equal(zip_returned, zip_temp)
})
