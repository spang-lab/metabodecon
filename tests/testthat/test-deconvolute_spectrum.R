# Inputs #####
defaults <- list(
    x = sap[[1]],
    nfit = 3, smit = 1, smws = 3, delta = 3, sfr = c(3.2, -3.2),
    force = FALSE, verbose = FALSE,
    use_rust = FALSE, nworkers = 1, igrs = list(), rtyp = "decon2"
)
args <- list(
    decon2_r    = set(defaults),
    rdecon_r    = set(defaults, rtyp = "rdecon"),
    rdecon_rust = set(defaults, rtyp = "rdecon", use_rust = TRUE)
)
mdrb_available <- check_mdrb()

# Helpers #####
try_deconvolute_spectrum <- function(args) {
    try(do.call(deconvolute_spectrum, args), silent = TRUE)
}

try_calc_prarpx <- function(obj) {
    try(calc_prarpx(obj), silent = TRUE)
}

# Calls #####
obj <- sapply(args, try_deconvolute_spectrum, simplify = FALSE)
prarp <- sapply(obj, try_calc_prarpx, simplify = FALSE)

# Checks #####
r_structures <- test_that(
    "Returned decon objects have correct structure with R backend", {
    expect_identical(object = names(obj$decon2_r), expected = decon2_members)
    expect_identical(object = class(obj$decon2_r), expected = "decon2")
    expect_identical(object = class(obj$rdecon_r), expected = "try-error")
})

r_prarps <- test_that(
    "PRARPs are good with R backend", {
    expect_true(prarp$decon2_r >= 0.961) # (1)
    expect_true(inherits(prarp$rdecon_r, "try-error"))
    # (1) MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
})

igrs_test <- test_that(
    "Peaks in ignore regions are excluded", {
    # sap[[1]] has 4 peaks (delta=3); 1 peak is near ppm ~0.1
    igrs <- list(c(-0.5, 0.5))
    d_r <- deconvolute_spectrum(
        sap[[1]], nfit = 3, smit = 1, smws = 3, delta = 3,
        sfr = c(3.2, -3.2), force = FALSE,
        verbose = FALSE, use_rust = FALSE, igrs = igrs, rtyp = "decon2"
    )
    n_no_igrs  <- nrow(obj$decon2_r$lcpar)
    n_with_igrs <- nrow(d_r$lcpar)
    expect_equal(n_no_igrs, 4)
    expect_equal(n_with_igrs, 3)
    if (mdrb_available) {
        d_rust <- deconvolute_spectrum(
            sap[[1]], nfit = 3, smit = 1, smws = 3, delta = 3,
            sfr = c(3.2, -3.2), force = FALSE,
            verbose = FALSE, use_rust = TRUE, igrs = igrs, rtyp = "decon2"
        )
        expect_equal(nrow(d_rust$lcpar), 3)
    }
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

rust_structures <- test_that(
    "Returned decon objects have correct structure with Rust backend", {
    expect_identical(object = names(obj$rdecon_rust), expected = rdecon_members)
    expect_identical(object = class(obj$rdecon_rust), expected = "rdecon")
})

rust_prarps <- test_that(
    "PRARPs are good with Rust backend", {
    expect_true(prarp$rdecon_rust >= 0.961)
    expect_true(prarp$decon2_r <= prarp$rdecon_rust) # Rust >= R
})

r_rust_comparison <- test_that(
    "R and Rust produce similar quality results", {
    # Make sure the objects are plottable
    expect_no_error(evalwith(plot = "captured", {
        plot_spectrum(obj$decon2_r)
        plot_spectrum(obj$rdecon_rust)
    }))
})

