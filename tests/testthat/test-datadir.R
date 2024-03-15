test_that("1. datadir(persistent = TRUE, warn = FALSE)", {
    x <- datadir(persistent = TRUE, warn = FALSE)
    y <- datadir_persistent()
    expect_equal(x, y)
})

test_that("2. datadir(persistent = FALSE, warn = FALSE)", {
    x <- datadir(persistent = FALSE, warn = FALSE)
    y <- datadir_temp()
    expect_equal(x, y)
})

test_that("3. datadir() with existing <datadir_persistent>/example_datasets.zip", {
    x <- evalwith(datadir_persistent = "filled", {
        dir_returned <- datadir()
        dir_expected <- datadir(persistent = TRUE)
    })
    expect_equal(dir_returned, dir_expected)
})

test_that("4. datadir(file = 'non_existent_file')", {
    x <- evalwith(message = "captured", {
        datadir(file = "non_existent_file")
    })
    expected_warning <- paste("Warning:", x$rv, "does not exist. Please call `download_example_datasets()` first.")
    expect_equal(x$message$text, expected_warning)
})

test_that("5. datadir(warn = FALSE) with missing <datadir_persistent> and <datadir_temp>", {
    x <- evalwith(datadir_persistent = "missing", datadir_temp = "missing", message = "captured", {
        datadir(file = "non_existent_file", warn = FALSE)
    })
    expect_equal(x$message$text, character())
})
