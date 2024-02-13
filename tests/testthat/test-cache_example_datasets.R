test_that("1. cache_example_datasets(persistent = FALSE) with empty datadirs", {
    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    x <- with(
        datadir_persistent = "empty",
        datadir_temp = "empty",
        message = "captured",
        expr = {
            zip_returned <- cache_example_datasets(persistent = FALSE)
            zip_temp <- zip_temp()
            zip_persistent <- zip_persistent()
        }
    )

    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    expect_equal(x$message$text, paste("Downloading", xds$url, "as", zip_temp))
    expect_equal(zip_returned, zip_temp)
})

test_that("2. cache_example_datasets(persistent = FALSE) with existing temporary datadir", {
    x <- with(
        datadir_persistent = "empty",
        datadir_temp = "filled",
        message = "captured",
        expr = {
            zip_returned <- cache_example_datasets(persistent = FALSE)
            zip_temp <- zip_temp()
            zip_persistent <- zip_persistent()
        }
    )
    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    expect_equal(x$message$text, character())
    expect_equal(zip_returned, zip_temp)
})
