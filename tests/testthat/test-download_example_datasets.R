test_that("download_example_datasets downloads and extracts files correctly", {
    with_mocked_persistent_datadir({
        datadir_temp <- clear(datadir_temp())
        datadir_persistent <- clear(datadir_persistent())
        expect_equal(dir(datadir_temp), character())
        expect_equal(dir(datadir_persistent), character())
        download_example_datasets()

    })
})

test_that("download_example_datasets handles cache_persistent and ask parameters correctly", {
    # Mock the functions to avoid actual download and extraction
    mock_cache_example_datasets <- function(cache = NULL, ask = FALSE) {
        return(list(cache = cache, ask = ask))
    }
    withr::local_mock(list(cache_example_datasets = mock_cache_example_datasets))

    # Test with cache_persistent = TRUE and ask = FALSE
    result <- download_example_datasets(cache_persistent = TRUE, ask = FALSE)
    expect_equal(result$cache, TRUE)
    expect_equal(result$ask, FALSE)

    # Test with cache_persistent = FALSE and ask = TRUE
    result <- download_example_datasets(cache_persistent = FALSE, ask = TRUE)
    expect_equal(result$cache, FALSE)
    expect_equal(result$ask, TRUE)
})
