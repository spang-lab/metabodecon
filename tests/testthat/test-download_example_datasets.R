library(testthat)

skip_if_slow_tests_disabled()

test_that("download_example_datasets works if xdszip=cached", {
    system.time({
        x <- evalwith(
            datadir_persistent = "filled",
            datadir_temp = "missing",
            message = "captured",
            expr = {
                download_example_datasets(persistent = TRUE)
                expected_path <- file.path(datadir(), "example_datasets.zip")
            }
        )
        expect_true(file.exists(expected_path))
        expect_equal(file.size(expected_path), xds$zip_size)
        expect_equal(object = x$message, expected = character())
    })
})

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

test_that("example datasets v1.1.0 and v1.6.3 differ only by aki dataset", {
    testthat::skip_on_cran()

    xds_new <- xds
    xds_old <- list(
        url = "https://github.com/spang-lab/metabodecon/releases/download/v1.1.0/example_datasets.zip",
        zip_size = 38425397,
        dir_size = 56378684,
        n_files = 1018
    )

    dst_dir <- file.path(tempdir(), "example-dataset-comparison")
    dst_dir_old <- file.path(dst_dir, "old")
    xds_zip_old <- file.path(dst_dir, "old/example_datasets.zip")
    xds_dir_old <- file.path(dst_dir, "old/example_datasets")
    dst_dir_new <- file.path(dst_dir, "new")
    xds_zip_new <- file.path(dst_dir, "new/example_datasets.zip")
    xds_dir_new <- file.path(dst_dir, "new/example_datasets")

    download_example_datasets_zip(path = xds_zip_old, url = xds_old$url, zip_size = xds_old$zip_size, silent = TRUE)
    download_example_datasets_zip(path = xds_zip_new, url = xds_new$url, zip_size = xds_new$zip_size, silent = TRUE)
    utils::unzip(xds_zip_old, exdir = dst_dir_old)
    utils::unzip(xds_zip_new, exdir = dst_dir_new)

    old_files <- list.files(xds_dir_old, full.names = FALSE, recursive = TRUE)
    new_files <- list.files(xds_dir_new, full.names = FALSE, recursive = TRUE)
    old_files <- old_files[file.info(file.path(xds_dir_old, old_files))$isdir == FALSE]
    new_files <- new_files[file.info(file.path(xds_dir_new, new_files))$isdir == FALSE]
    old_files <- old_files[grepl("^(bruker|jcampdx)/", old_files)]
    new_files <- new_files[grepl("^(bruker|jcampdx)/", new_files)]

    removed <- setdiff(old_files, new_files)
    added <- setdiff(new_files, old_files)
    common <- sort(intersect(old_files, new_files))

    expect_length(removed, 0)
    expect_true(length(added) > 0)
    expect_true(all(grepl("^bruker/aki/", added)))

    old_common_files <- file.path(xds_dir_old, common)
    new_common_files <- file.path(xds_dir_new, common)
    old_hash <- unname(tools::md5sum(old_common_files))
    new_hash <- unname(tools::md5sum(new_common_files))
    changed <- common[old_hash != new_hash]
    expect_true(length(changed) == 0, info = paste("Changed files:", paste(changed, collapse = ", ")))
})
