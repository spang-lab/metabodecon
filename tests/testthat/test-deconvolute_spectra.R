library(testthat)

# IMPORTANT: in this file we don't test for good PRARPs, as these checks are
# already done in `test-deconvolute_spectrum` (singular). Here, we only want to
# assert that the "glue code" in `test-deconvolute_spectra` (plural),
# orchestrating multiple invocations of `deconvolute_spectrum`, works as
# expected. We do this by checking that the returned objects have the correct
# types.

# Inputs #####
defaults <- list(
    x=sap,
    nfit=3, smopts=c(1,3), delta=3, sfr=c(3.2,-3.2), wshw=0,
    ask=FALSE, force=FALSE, verbose=FALSE, bwc=0,
    use_rust=FALSE, nw=1, igr=list(), rtyp="idecon"
)

args_idecons_bwc0_R <- set(defaults)
args_idecons_bwc1_R <- set(defaults, bwc=1)
args_idecons_bwc2_R <- set(defaults, bwc=2)
args_decons0_bwc2_R <- set(defaults, bwc=2, rtyp="decon0")
args_decons1_bwc2_R <- set(defaults, bwc=2, rtyp="decon1")
args_decons2_bwc2_R <- set(defaults, bwc=2, rtyp="decon2")
args_rdecons_bwc2_R <- set(defaults, bwc=2, rtyp="rdecon")
args_rdecons_bwc2_rust <- set(defaults, bwc=2, rtyp="rdecon", use_rust=TRUE)

# Calls #####
idecons_bwc0_R <- try(do.call(deconvolute_spectra, args_idecons_bwc0_R), silent = TRUE)
idecons_bwc1_R <- try(do.call(deconvolute_spectra, args_idecons_bwc1_R), silent = TRUE)
idecons_bwc2_R <- try(do.call(deconvolute_spectra, args_idecons_bwc2_R), silent = TRUE)
decons0_bwc2_R <- try(do.call(deconvolute_spectra, args_decons0_bwc2_R), silent = TRUE)
decons1_bwc2_R <- try(do.call(deconvolute_spectra, args_decons1_bwc2_R), silent = TRUE)
decons2_bwc2_R <- try(do.call(deconvolute_spectra, args_decons2_bwc2_R), silent = TRUE)
rdecons_bwc2_R <- try(do.call(deconvolute_spectra, args_rdecons_bwc2_R), silent = TRUE)
rdecons_bwc2_rust <- try(do.call(deconvolute_spectra, args_rdecons_bwc2_rust), silent = TRUE)

# Checks #####
types <- test_that("Types are ok", {
    expect_identical(class(idecons_bwc0_R), "idecons")
    expect_identical(class(idecons_bwc1_R), "idecons")
    expect_identical(class(idecons_bwc2_R), "idecons")
    expect_identical(class(decons0_bwc2_R), "list")
    expect_identical(class(decons1_bwc2_R), "decons1")
    expect_identical(class(decons2_bwc2_R), "decons2")
    expect_identical(class(rdecons_bwc2_R), "try-error")
    expect_identical(class(rdecons_bwc2_rust), "rdecons")

    expect_identical(class(idecons_bwc0_R[[1]]), "idecon")
    expect_identical(class(idecons_bwc1_R[[1]]), "idecon")
    expect_identical(class(idecons_bwc2_R[[1]]), "idecon")
    expect_identical(class(decons0_bwc2_R[[1]]), "list")
    expect_identical(class(decons1_bwc2_R[[1]]), "decon1")
    expect_identical(class(decons2_bwc2_R[[1]]), "decon2")
    expect_identical(class(rdecons_bwc2_rust[[1]]), "rdecon")

    expect_identical(names(idecons_bwc0_R[[1]]), idecon_members)
    expect_identical(names(idecons_bwc1_R[[1]]), idecon_members)
    expect_identical(names(idecons_bwc2_R[[1]]), idecon_members)
    expect_identical(names(decons1_bwc2_R[[1]]), decon1_members)
    expect_identical(names(decons2_bwc2_R[[1]]), decon2_members)
    expect_identical(names(rdecons_bwc2_rust[[1]]), rdecon_members)
})
