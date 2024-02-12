# # Inputs:
# # persistent: TRUE/FALSE/NULL
# # ask: TRUE/FALSE
# # persdir_exists: TRUE/FALSE
# # perszip_exists: TRUE/FALSE
# # perszip_correct_size: TRUE/FALSE
# # tempzip_exists: TRUE/FALSE
# # tempzip_correct_size: TRUE/FALSE

# XDG_DATA_HOME_BAK <- mock_persistent_datadir()
# on.exit(restore_persistent_datadir(XDG_DATA_HOME_BAK), add = TRUE)

# datadir_persistent <- datadir_persistent()
# datadir_temp <- datadir_temp()
# zip_persistent <- file.path(datadir_persistent, "example_datasets.zip")
# zip_temp <- file.path(datadir_temp, "example_datasets.zip")
# xds$url <- "https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip"


test_that("1. cache_example_datasets(persistent = FALSE) with empty datadirs", {

    skip_if_not(Sys.getenv("RUN_SLOW_TESTS") == "TRUE", "Skipped because RUN_SLOW_TESTS != TRUE")

    x <- with(datadir_persistent = "empty", datadir_temp = "empty", message = "captured", {
        zip_returned <- cache_example_datasets(persistent = FALSE)
        zip_temp <- zip_temp()
        zip_persistent <- zip_persistent()
    })

    expected_message <- paste("Downloading", xds$url, "as", zip_temp)

    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    expect_equal(x$message$text, expected_message)
    expect_equal(zip_returned, zip_temp)
})

test_that("2. cache_example_datasets(persistent = FALSE) with existing temporary datadir", {

    x <- with(datadir_persistent = "empty", datadir_temp = "filled", message = "captured", {
        zip_returned <- cache_example_datasets(persistent = FALSE)
        zip_temp <- zip_temp()
        zip_persistent <- zip_persistent()
    })

    expected_message <- character()

    expect_true(file.exists(zip_temp))
    expect_false(file.exists(zip_persistent))
    expect_equal(x$message$text, expected_message)
    expect_equal(zip_returned, zip_temp)
})

# test_that("`cache_example_datasets(persistent = NULL, ask = TRUE)` works", {

#     mock_readline(answers = "n")
#     mock_stderr()
#     try(path <- cache_example_datasets())
#     restore_stderr()
#     restore_readline()

#     stderr <- .Options$metabodecon.mocks.stderr$vec
#     stderr_expected <- paste("Cache example datasets (38.3 MB) permanently on disk at", zip_persistent, "to speed up future calls to this function? (y/n) n")
#     expect_equal(stderr, stderr_expected)
#     expect_equal(path, zip_temp)
#     expect_true(file.exists(zip_temp))
#     expect_false(file.exists(zip_persistent))

#     mock_readline(answers = "y")
#     mock_stderr()
#     try(path <- cache_example_datasets())
#     restore_stderr()
#     restore_readline()

#     stderr <- .Options$metabodecon.mocks.stderr$vec
#     stderr_expected <- c()
#     stderr_expected[1] <- paste("Cache example datasets (38.3 MB) permanently on disk at", zip_persistent, "to speed up future calls to this function? (y/n) y")
#     stderr_expected[2] <- paste("Downloading", xds$url, "as", zip_persistent)
#     expect_equal(stderr, stderr_expected)
#     expect_equal(path, zip_persistent)
#     expect_true(file.exists(zip_temp))
#     expect_true(file.exists(zip_persistent))
# })

# test_that("`cache_example_datasets(persistent = NULL, ask = FALSE)` works", {

#     clear(datadir_persistent)
#     stdout_expected <- c(paste("Downloading", xds$url, "as", zip_temp))
#     stdout <- capture.output(type = "message", {path <- cache_example_datasets(ask = FALSE)})
#     expect_equal(stdout, stdout_expected)
#     expect_equal(path, zip_temp)
#     expect_true(file.exists(zip_temp))
#     expect_false(file.exists(zip_persistent))

#     # Test `cache_example_datasets(persistent = TRUE, ask = FALSE)`
#     stdout_expected <- c(paste("Downloading", xds$url, "as", zip_persistent))
#     stdout <- capture.output(type = "message", {path <- cache_example_datasets(ask = FALSE)})
#     expect_equal(stdout, stdout_expected)
#     expect_equal(path, zip_persistent)
#     expect_true(file.exists(zip_temp))
#     expect_false(file.exists(zip_persistent))
# })
