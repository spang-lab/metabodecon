xds$url <- "https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip"

test_that("1. download_example_datasets() with example_datasets.zip missing on disk", {

    skip_if(Sys.getenv("SKIP_SLOW_TESTS") == "TRUE", "Skipped because SKIP_SLOW_TESTS=TRUE")

    x <- with(datadir_persistent = "missing", datadir_temp = "missing", message = "captured", {
        download_example_datasets()
        expected_path <- file.path(datadir(), "example_datasets.zip")
    })
    expected_message <- paste("Downloading", xds$url, "as", expected_path)

    expect_equal(file.exists(expected_path), TRUE)
    expect_equal(file.size(expected_path), xds$zip_size)
    expect_equal(x$message$text, expected_message)
})

test_that("2. download_example_datasets(persistent = TRUE) with example_datasets.zip missing on disk", {

    skip_if(Sys.getenv("SKIP_SLOW_TESTS") == "TRUE", "Skipped because SKIP_SLOW_TESTS=TRUE")

    x <- with(datadir_persistent = "missing", datadir_temp = "missing", message = "captured", {
        download_example_datasets(persistent = TRUE)
        expected_path <- file.path(datadir(), "example_datasets.zip")
    })
    expected_message <- paste("Downloading", xds$url, "as", expected_path)

    expect_equal(file.exists(expected_path), TRUE)
    expect_equal(file.size(expected_path), xds$zip_size)
    expect_equal(x$message$text, expected_message)
})

test_that("3. download_example_datasets() with example_datasets.zip already cached on disk", {

    skip_if(Sys.getenv("SKIP_SLOW_TESTS") == "TRUE", "Skipped because SKIP_SLOW_TESTS=TRUE")

    x <- with(datadir_persistent = "filled", datadir_temp = "missing", message = "captured", {
        download_example_datasets(persistent = TRUE)
        expected_path <- file.path(datadir(), "example_datasets.zip")
    })

    expect_true(file.exists(expected_path))
    expect_equal(file.size(expected_path), xds$zip_size)
    expect_equal(x$message$text, character())
})