VERBOSE <- FALSE

test_that("with_redirects works", {
  output_dir <- prepare_test_dir("with_redirects", "1", verbose = VERBOSE)
  owd <- setwd(output_dir)
  on.exit(setwd(owd), add = TRUE)
  with_redirects(stdout = "out.txt", stderr = "err.txt", plots = "plots.pdf", expr = {
    cat2("Starting plotting")
    plot(1:10)
    cat2("Done plotting")
    message("This goes to stderr")
  })
  cks_out <- checksum(output_dir)
  cks_exp <- c(err.txt=21, out.txt=34, plots.pdf=4836)
  expect_identical(cks_out, cks_exp)
})
