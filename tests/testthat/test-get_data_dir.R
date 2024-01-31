existing_temp_dir <- file.path(tempdir(), timestamp)
dir.create(existing_temp_dir, recursive = TRUE)

test_that("datadir returns correct persistent path", {
  x <- datadir(persistent = TRUE, warn = FALSE)
  y <- datadir_persistent()
  expect_equal(x, y)
})

test_that("datadir returns correct temp path", {
  x <- datadir(persistent = FALSE, warn = FALSE)
  y <- datadir_temp()
  expect_equal(x, y)
})

test_that("datadir returns correct path if '{datadir_persistent}' exists", {
  with_mock_datadirs({
    x <- datadir()
    y <- datadir_persistent()
    expect_equal(x, y)
  })
})

test_that("datadir returns correct path if '{datadir_persistent}' does not exist", {
  with_mock_datadirs(persistent = FALSE, {
    x <- datadir()
    y <- datadir_temp()
    expect_equal(x, y)
  })
})

test_that("datadir issues a warning when the file does not exist", {
  expect_warning(datadir(file = "non_existent_file"))
})

test_that("datadir does not issue a warning when warn is FALSE", {
  expect_silent(datadir(file = "non_existent_file", warn = FALSE))
})
