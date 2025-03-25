# Inputs #####
defaults <- list(
    x=sap[[1]],
    nfit=3, smopts=c(1,3), delta=3, sfr=c(3.2,-3.2), wshw=0,
    ask=FALSE, force=FALSE, verbose=FALSE, bwc=0,
    use_rust=FALSE, nw=1, igr=list(), rtyp="idecon"
)
args <- list(
    idecon_bwc0 = set(defaults),
    idecon_bwc1 = set(defaults, bwc=1),
    idecon_bwc2 = set(defaults, bwc=2),
    idecon_rust = set(defaults, bwc=2, use_rust=TRUE),
    decon2_bwc2 = set(defaults, bwc=2, rtyp="decon2"),
    decon2_rust = set(defaults, bwc=2, rtyp="decon2", use_rust=TRUE)
)
mdrb_available <- check_mdrb()

# Helpers #####
try_deconvolute_spectrum <- function(args) {
    try(do.call(deconvolute_spectrum, args), silent=TRUE)
}

try_calc_prarpx <- function(obj) {
    try(calc_prarpx(obj), silent=TRUE)
}

# Calls #####
obj <- sapply(args, try_deconvolute_spectrum, simplify=FALSE)
prarp <- sapply(obj, try_calc_prarpx, simplify=FALSE)

# Checks #####
r_structures <- test_that(
    "Returned decon objects have correct structure with R backend", {
    # We only test 'idecon' and 'decon2' objects, as that's enough to
    # cover the R and Rust backends. Conversions to other types are tested
    # in test-as_decon.R.
    expect_identical(object = names(obj$idecon_bwc0), expected = idecon_members)
    expect_identical(object = names(obj$idecon_bwc1), expected = idecon_members)
    expect_identical(object = names(obj$idecon_bwc2), expected = idecon_members)
    expect_identical(object = names(obj$decon2_bwc2), expected = decon2_members)
    expect_identical(object = class(obj$idecon_bwc0), expected = "idecon")
    expect_identical(object = class(obj$idecon_bwc1), expected = "idecon")
    expect_identical(object = class(obj$idecon_bwc2), expected = "idecon")
    expect_identical(object = class(obj$decon2_bwc2), expected = "decon2")
})

r_prarps <- test_that(
    "PRARPs are good with R backend", {
    expect_true(prarp$idecon_bwc0 >= 0.507) # (1)
    expect_true(prarp$idecon_bwc1 >= 0.961)
    expect_true(prarp$idecon_bwc2 >= 0.961)
    expect_true(prarp$decon2_bwc2 >= 0.961)
    expect_true(prarp$idecon_bwc0 <= prarp$idecon_bwc1) # (2)
    expect_true(prarp$idecon_bwc1 <= prarp$idecon_bwc2) # (2)
    expect_true(prarp$idecon_bwc2 == prarp$decon2_bwc2) # (3)
    # (1) MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
    # (2) Higher bwc versions should lead to higher PRARPs
    # (3) Different return types should not affect PRARPs
})

bwc <- test_that(
    "Deconvolute with bwc=0 returns the same result as MetaboDecon1D", {
    decon0_deconvolute <- as_decon0(obj$idecon_bwc0, optional = FALSE)
    decon0_MetaboDecon1D <- MetaboDecon1D_silent(
        filepath = metabodecon_file("bruker/sap"),
        filename = "sap_01",
        number_iterations = 3,
        range_water_signal_ppm = 0,
        signal_free_region = c(3.2, -3.2),
        smoothing_param = c(1, 3),
        delta = 3
    )
    spec <- read_spectrum(metabodecon_file("bruker/sap/sap_01"))
    expect_equal(decon0_deconvolute, decon0_MetaboDecon1D)
})

# Rust Checks #####

skip_on_cran()
skip_if(getRversion() < numeric_version("4.2"))

mdrb <- test_that(
    "MDRB is available", {
    expect_true(mdrb_available)
})

skip_if_not(mdrb_available) # (1)
# (1) If we reach this point, we're not on CRAN and our R version is greater
# equal 4.2. I.e., mdrb should be available. If it is not, the "MDRB is
# available" check from above will fail and that's enough for us to see that
# something is wrong. I.e, in such as scenario, there is no need to execute the
# following tests and spam the log file.

rust_structures <- test_that(
    "Returned decon objects have correct structure with Rust backend", {
    expect_identical(object = names(obj$idecon_rust), expected = idecon_members)
    expect_identical(object = names(obj$decon2_rust), expected = decon2_members)
    expect_identical(object = class(obj$idecon_rust), expected = "idecon")
    expect_identical(object = class(obj$decon2_rust), expected = "decon2")
})

rust_prarps <- test_that(
    "PRARPs are good with Rust backend", {
    expect_true(prarp$idecon_rust >= 0.961)
    expect_true(prarp$decon2_rust >= 0.961)
    expect_true(prarp$idecon_bwc2 <= prarp$idecon_rust) # (1)
    expect_true(prarp$idecon_bwc2 <= prarp$decon2_rust) # (1)
    # (1) Rust should be better or equal to R
})

rust_r_equality <- test_that(
    "Rust and R backend produce equal results", {

    # Update all fields from the idecon_bwc2 object for which differences
    # are to be expected. After patching these fields, the objects should
    # be equal.
    idecon_bwc2_mod <- obj$idecon_bwc2
    idecon_rust_mod <- obj$idecon_rust
    idecon_bwc2_mod$args$use_rust <- TRUE
    idecon_bwc2_mod[c("peak", "Z", "lci", "lca", "lcr")] <- NULL
    idecon_rust_mod[c("peak", "Z", "lci", "lca", "lcr")] <- NULL
    expect_equal(idecon_rust_mod, idecon_bwc2_mod)

    # Altough the estimated Lorentzian Parameters are expected to be slightly
    # different, the PRARPs should be roughly the same. In fact, the rust prarp
    # should be higher, as it already uses the raw SIs for the approximation
    # steps, which is not yet implemented in the R version.
    expect_true(prarp$idecon_rust >= prarp$idecon_bwc2)

    # We also expect the rust version to produce an object that is plottable
    expect_no_error(evalwith(plot = "captured",
        plot_spectrum(obj$idecon_rust)
    ))
})