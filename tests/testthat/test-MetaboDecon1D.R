test_that("Deconvolution of 1 jcampdx works", {
  skip_if(Sys.getenv("SKIP_SLOW_TESTS") == "TRUE", "Skipped because SKIP_SLOW_TESTS=TRUE")
  output_dir <- prepare_test_dir(
    fn = "MetaboDecon1D",
    tc = "1",
    inputs = c(urine.dx = "example_datasets/jcampdx/urine/urine.dx")
  )
  runtime_txt <- file.path(output_dir, "runtime.txt")
  readline_mock <- get_readline_mock(texts = c(
    "y", # Signal free region borders correct selected? ... (y/n)
    "y" # Water artefact fully inside red vertical lines? (y/n)
  ))
  owd <- setwd(output_dir)
  on.exit(setwd(owd), add = TRUE)
  with_redirects(stdout = "out.txt", stderr = "err.txt", plots = "plots.pdf", expr = {
    testthat::with_mocked_bindings(readline = readline_mock, code = {
      start <- Sys.time()
      set.seed(123)
      MetaboDecon1D(
        filepath = output_dir,
        filename = "urine.dx",
        file_format = "jcampdx",
        number_iterations = 1
      )
      end <- Sys.time()
      runtime <-  round(as.numeric(end - start, units = "secs"), 2)
      write(paste(runtime, "seconds"), file = runtime_txt)
    })
  })
  output_checksum <- checksum(output_dir, ignore = "runtime.txt")
  expect_checksum <- c()
  expect_identical(output_checksum, expect_checksum)
  expect_equal(1, 1)
})
