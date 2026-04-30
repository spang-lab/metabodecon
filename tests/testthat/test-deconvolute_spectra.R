library(testthat)

# IMPORTANT: in this file we don't test for good PRARPs, as these checks are
# already done in `test-deconvolute_spectrum` (singular). Here, we only want to
# assert that the "glue code" in `test-deconvolute_spectra` (plural),
# orchestrating multiple invocations of `deconvolute_spectrum`, works as
# expected. We do this by checking that the returned objects have the correct
# types.

# Inputs #####
defaults <- list(
    x = sap,
    nfit = 3, smit = 1, smws = 3, delta = 3, sfr = c(3.2, -3.2),
    force = FALSE, verbose = FALSE,
    use_rust = FALSE, nworkers = 1, igrs = list(), rtyp = "decon2"
)
args <- list(
    decons0_R = set(defaults, rtyp = "decon0"),
    decons1_R = set(defaults, rtyp = "decon1"),
    decons2_R = set(defaults, rtyp = "decon2"),
    rdecons_rust = set(defaults, rtyp = "rdecon", use_rust = TRUE)
)
mdrb_available <- check_mdrb()

# Helpers #####
try_deconvolute_spectra <- function(args) {
    try(do.call(deconvolute_spectra, args), silent = TRUE)
}

# Calls #####
obj <- sapply(args, try_deconvolute_spectra, simplify = FALSE)

# Checks #####
r_return_types <- test_that("R return types are ok", {

    expect_identical(class(obj$decons0_R), "list")
    expect_identical(class(obj$decons1_R), "decons1")
    expect_identical(class(obj$decons2_R), "decons2")

    expect_identical(class(obj$decons0_R[[1]]), "list")
    expect_identical(class(obj$decons1_R[[1]]), "decon1")
    expect_identical(class(obj$decons2_R[[1]]), "decon2")

    expect_identical(names(obj$decons1_R[[1]]), decon1_members)
    expect_identical(names(obj$decons2_R[[1]]), decon2_members)
})

# Rust Checks #####

skip_on_cran()
skip_if(getRversion() < numeric_version("4.2"))
skip_if_not(mdrb_available) # (1)
# (1) If we reach this point, we're not on CRAN and our R version is greater
# equal 4.2. I.e., mdrb should be available. If it is not, the "MDRB is
# available" check from `test-deconvolute.R` will fail and that's enough for us
# to see that something is wrong. I.e, in such as scenario, there is no need to
# execute the following tests and spam the log file.

rust_return_types <- test_that("Rust return types are ok", {
    expect_identical(class(obj$rdecons_rust), "rdecons")
    expect_identical(class(obj$rdecons_rust[[1]]), "rdecon")
    expect_identical(names(obj$rdecons_rust[[1]]), rdecon_members)
})
