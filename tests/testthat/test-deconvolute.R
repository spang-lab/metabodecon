library(testthat)

# Inputs #####

sap1 <- sap[[1]]

# Helpers #####

deconvolute_sap1 <- function(use_rust=FALSE) {
    deconvolute(sap1, smopts=c(1,3), delta=3, sfr=c(3.2,-3.2), use_rust=use_rust)
}

expect_sap1_deconvolution_worked <- function(decon2) {
    expect_identical(object=names(decon2), expected=decon2_members)
    expect_identical(object=class(decon2), expected="decon2")
    expect_true(calc_prarpx(decon2) >= 0.961) # (1)
    # (1) MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
}

expect_sap1_deconvolution_failed <- function(decon2, expected_error) {
    expect_identical(object=class(decon2), expected="try-error")
    observed_error <- attributes(decon2)$condition$message
    expect_match(observed_error, expected_error)
}

get_zero_version <- function() {
    package_version("0.0.0")
}

# Tests #####
test_sap1 <- test_that("deconvolute works for single spectrum", {
    decon2 <- deconvolute_sap1()
    expect_sap1_deconvolution_worked(decon2)
})

test_sim_subset <- test_that("deconvolute works for multiple spectra", {
    decons2 <- deconvolute(x = sim[1:2], sfr = c(3.55, 3.35))
    expect_identical(class(decons2), "decons2")
    expect_identical(class(decons2[[1]]), "decon2")
    expect_identical(class(decons2[[2]]), "decon2")
    expect_identical(names(decons2), c("sim_01", "sim_02"))
    expect_identical(names(decons2[[1]]), decon2_members)
    expect_identical(names(decons2[[2]]), decon2_members)
    obj1 <- calc_prarp(decons2[[1]], truepar = sim[[1]]$meta$simpar)
    obj2 <- calc_prarp(decons2[[2]], truepar = sim[[2]]$meta$simpar)
    expect_true(obj1$prarpx >= 0.777) # MetaboDecon1D has a PRARPX of 0.732. See test-MetaboDecon1d.R.
    expect_true(obj2$prarpx >= 0.750) # MetaboDecon1D has a PRARPX of 0.710. See test-MetaboDecon1d.R.
})

test_wrong_sfr <- test_that("deconvolute works when no peaks are filtered out", {
    x <- simulate_spectrum(ndp=256, npk=3)
    expect_error(deconvolute(x, sfr=c(Inf,-Inf), smopts=c(0,3)))
    deconForc <- deconvolute(x, sfr=c(Inf,-Inf), smopts=c(0,3), force=TRUE)
})

# Local Tests #####

skip_on_cran()

rust_backend <- test_that("deconvolute works with rust backend", {
    # Mdrb missing
    with_mocked_bindings(get_mdrb_version=get_zero_version, code = {
        rt_mm <- try(deconvolute_sap1(use_rust=TRUE), silent=TRUE)
        rn_mm <- try(deconvolute_sap1(use_rust=NULL), silent=TRUE)
        rf_mm <- try(deconvolute_sap1(use_rust=FALSE), silent=TRUE)
    })
    expect_sap1_deconvolution_failed(rt_mm, "Using.*Rust.*requires mdrb.*")
    expect_sap1_deconvolution_worked(rn_mm)
    expect_sap1_deconvolution_worked(rf_mm)
    # Mdrb available
    if (!check_mdrb()) {
        message("\nPackage 'mdrb' is missing. Attempting to install it.")
        install_mdrb(ask = FALSE, verbose = FALSE)
    }
    if (mdrb_available <- check_mdrb()) {
        rt_ma <- try(deconvolute_sap1(use_rust=TRUE), silent=TRUE)
        rn_ma <- try(deconvolute_sap1(use_rust=NULL), silent=TRUE)
        rf_ma <- try(deconvolute_sap1(use_rust=FALSE), silent=TRUE)
        expect_sap1_deconvolution_worked(rn_mm)
        expect_sap1_deconvolution_worked(rn_mm)
        expect_sap1_deconvolution_worked(rn_mm)
    } else {
        warning(paste(
            "Cannot test Rust backend because 'mdrb' is missing.",
            "Call `install_mdrb()` first and then and try again."
        ))
        expect_true(mdrb_available)
    }
})
