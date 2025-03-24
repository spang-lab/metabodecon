# Inputs #####
defaults <- list(
    x=sap[[1]],
    nfit=3, smopts=c(1,3), delta=3, sfr=c(3.2,-3.2), wshw=0,
    ask=FALSE, force=FALSE, verbose=FALSE, bwc=0,
    use_rust=FALSE, nw=1, igr=list(), rtyp="idecon"
)
args_idecon_bwc0 <- set(defaults)
args_idecon_bwc1 <- set(defaults, bwc=1)
args_idecon_bwc2 <- set(defaults, bwc=2)
args_idecon_rust <- set(defaults, bwc=2, use_rust=TRUE)
args_decon2_bwc2 <- set(defaults, bwc=2, rtyp="decon2")
args_decon2_rust <- set(defaults, bwc=2, rtyp="decon2", use_rust=TRUE)

# Calls #####
idecon_bwc0 <- do.call(deconvolute_spectrum, args_idecon_bwc0)
idecon_bwc1 <- do.call(deconvolute_spectrum, args_idecon_bwc1)
idecon_bwc2 <- do.call(deconvolute_spectrum, args_idecon_bwc2)
idecon_rust <- do.call(deconvolute_spectrum, args_idecon_rust)
decon2_bwc2 <- do.call(deconvolute_spectrum, args_decon2_bwc2)
decon2_rust <- do.call(deconvolute_spectrum, args_decon2_rust)

# Checks #####
structure <- test_that("Returned idecon objects have correct structure", {
    # Objects of other types are tested in other files:
    # decon[0-2],idecon: test-as_decon.R
    # decon1: test-generate_lorentz_curves.R
    # decon2: test-deconvolute.R
    expect_identical(object = names(idecon_bwc0), expected = idecon_members)
    expect_identical(object = names(idecon_bwc1), expected = idecon_members)
    expect_identical(object = names(idecon_bwc2), expected = idecon_members)
    expect_identical(object = names(idecon_rust), expected = idecon_members)
    expect_identical(object = names(decon2_bwc2), expected = decon2_members)
    expect_identical(object = names(decon2_rust), expected = decon2_members)
    expect_identical(object = class(idecon_bwc0), expected = "idecon")
    expect_identical(object = class(idecon_bwc1), expected = "idecon")
    expect_identical(object = class(idecon_bwc2), expected = "idecon")
    expect_identical(object = class(idecon_rust), expected = "idecon")
    expect_identical(object = class(decon2_bwc2), expected = "decon2")
    expect_identical(object = class(decon2_rust), expected = "decon2")
})

prarps <- test_that("PRARPs are good", {
    prarp_idecon_bwc0 <- calc_prarp(x = idecon_bwc0)$prarpx
    prarp_idecon_bwc1 <- calc_prarp(x = idecon_bwc1)$prarpx
    prarp_idecon_bwc2 <- calc_prarp(x = idecon_bwc2)$prarpx
    prarp_idecon_rust <- calc_prarp(x = idecon_rust)$prarpx
    prarp_decon2_bwc2 <- calc_prarp(x = decon2_bwc2)$prarpx
    prarp_decon2_rust <- calc_prarp(x = decon2_rust)$prarpx

    expect_true(prarp_idecon_bwc0 >= 0.507) # MetaboDecon1D has a PRARPX of 0.507. See test-MetaboDecon1d.R.
    expect_true(prarp_idecon_bwc1 >= 0.961)
    expect_true(prarp_idecon_bwc2 >= 0.961)
    expect_true(prarp_idecon_rust >= 0.961)
    expect_true(prarp_decon2_bwc2 >= 0.961)
    expect_true(prarp_decon2_rust >= 0.961)


    expect_true(prarp_idecon_bwc0 <= prarp_idecon_bwc1)
    expect_true(prarp_idecon_bwc1 <= prarp_idecon_bwc2)
    expect_true(prarp_idecon_bwc2 <= prarp_idecon_rust)

    expect_true(prarp_idecon_bwc2 == prarp_decon2_bwc2)
    expect_true(prarp_idecon_bwc2 <= prarp_decon2_rust)
})

bwc <- test_that("Deconvolute with bwc=0 returns the same as MetaboDecon1D", {
    decon0_deconvolute <- as_decon0(idecon_bwc0, optional = FALSE)
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

backends <- test_that("Rust and R backend produce equal results", {
    # Update all fields from the idecon_bwc2 object for which differences
    # are to be expected. After patching these fields, the objects should
    # be equal.
    idecon_bwc2_mod <- idecon_bwc2
    idecon_rust_mod <- idecon_rust
    idecon_bwc2_mod$args_use_rust <- TRUE
    idecon_bwc2_mod[c("peak", "Z", "lci", "lca", "lcr")] <- NULL
    idecon_rust_mod[c("peak", "Z", "lci", "lca", "lcr")] <- NULL
    expect_equal(idecon_rust_mod, idecon_bwc2_mod)
    # Altough the estimated Lorentzian Parameters are expected to be slightly
    # different, the PRARPs should be roughly the same. In fact, the rust prarp
    # should be higher, as it already uses the raw SIs for the approximation
    # steps, which is not yet implemented in the R version.
    prarp_rust <- calc_prarp(idecon_rust)$prarpx
    prarp_bwc2 <- calc_prarp(idecon_bwc2)$prarpx
    expect_true(prarp_rust >= prarp_bwc2)
    # We also expect the rust version to produce an object that is plottable
    expect_no_error(evalwith(plot = "captured",
        plot_spectrum(idecon_rust)
    ))
})
