test_that("pkg_file returns correct path to DESCRIPTION", {
    expected <- system.file("DESCRIPTION", package = "metabodecon")
    actual <- pkg_file("DESCRIPTION")
    expect_equal(actual, expected)
})

test_that("pkg_file returns correct path to file in inst folder", {
    expected <- system.file("WORDLIST", package = "metabodecon")
    actual <- pkg_file("WORDLIST")
    expect_equal(actual, expected)
})

test_that("pkg_file returns correct path to package root directory", {
    expected <- system.file(package = "metabodecon")
    actual <- pkg_file()
    expect_equal(actual, expected)
})

test_that("pkg_file raises error if file does not exist and mustWork is TRUE", {
    expect_error(pkg_file("nonexistent_file.txt", mustWork = TRUE))
})
