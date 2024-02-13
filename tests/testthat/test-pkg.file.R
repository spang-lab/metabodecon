test_that("pkg.file returns correct path to DESCRIPTION", {
    expected <- system.file("DESCRIPTION", package = "metabodecon")
    actual <- pkg.file("DESCRIPTION")
    expect_equal(actual, expected)
})

test_that("pkg.file returns correct path to file in inst folder", {
    expected <- system.file("WORDLIST", package = "metabodecon")
    actual <- pkg.file("WORDLIST")
    expect_equal(actual, expected)
})

test_that("pkg.file returns correct path to package root directory", {
    expected <- system.file(package = "metabodecon")
    actual <- pkg.file()
    expect_equal(actual, expected)
})

test_that("pkg.file raises error if file does not exist and mustWork is TRUE", {
    expect_error(pkg.file("nonexistent_file.txt", mustWork = TRUE))
})
