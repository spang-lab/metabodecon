test_that("download_example_datasets works if xdszip=missing", {
    x <- evalwith(datadir_persistent = "missing", datadir_temp = "missing", message = "captured", {
        download_example_datasets()
        expected_path <- file.path(datadir(), "example_datasets.zip")
    })
    expected_message <- paste("Downloading", xds$url, "as", expected_path)
    expect_equal(file.exists(expected_path), TRUE)
    expect_equal(file.size(expected_path), xds$zip_size)
    expect_equal(x$message, expected_message)
})

test_that("download_example_datasets works if xdszip=missing and persistent=T", {
    x <- evalwith(datadir_persistent = "missing", datadir_temp = "missing", message = "captured", {
        download_example_datasets(persistent = TRUE)
        expected_path <- file.path(datadir(), "example_datasets.zip")
    })
    expected_message <- paste("Downloading", xds$url, "as", expected_path)
    expect_equal(file.exists(expected_path), TRUE)
    expect_equal(file.size(expected_path), xds$zip_size)
    expect_equal(x$message, expected_message)
})

test_that("download_example_datasets works if xdszip=cached", {
    x <- evalwith(datadir_persistent = "filled", datadir_temp = "missing", message = "captured", {
        download_example_datasets(persistent = TRUE)
        expected_path <- file.path(datadir(), "example_datasets.zip")
    })
    expect_true(file.exists(expected_path))
    expect_equal(file.size(expected_path), xds$zip_size)
    expect_equal(x$message, character())
})
