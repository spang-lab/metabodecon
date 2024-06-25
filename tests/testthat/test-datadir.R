test_that("datadir works if persistent=T, warn=F", {
    x <- datadir(persistent = TRUE, warn = FALSE)
    y <- datadir_persistent()
    expect_equal(x, y)
})

test_that("datadir works if persistent=FALSE, warn=FALSE", {
    x <- datadir(persistent = FALSE, warn = FALSE)
    y <- datadir_temp()
    expect_equal(x, y)
})

test_that("datadir works if datadir_persistent=filled", {
    evalwith(datadir_persistent = "filled", {
        x <- datadir()
        y <- datadir(persistent = TRUE)
    })
    expect_equal(x, y)
})

test_that("datadir works if file=non_existent_file", {
    x <- evalwith(message = "captured", {
        datadir(file = "non_existent_file")
    })
    expected_warning <- paste("Warning:", x$rv, "does not exist. Please call `download_example_datasets()` first.")
    expect_equal(x$message, expected_warning)
})

test_that("datadir works if warn=F, datadirs=missing", {
    x <- evalwith(datadir_persistent = "missing", datadir_temp = "missing", message = "captured", {
        datadir(file = "non_existent_file", warn = FALSE)
    })
    expect_equal(x$message, character())
})
