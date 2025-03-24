library(testthat)

# Test Strategy #####

# 1. Get the sap spectrum.
# 2. Deconvolute the sap spectrum using [deconvolute_spectrum()] into an idecon
#    object.
# 3. Convert the idecon object to decon0, decon1 and decon2 objects.
# 4. Verify the values of the idecon object and the decon[0-2] objects by
#    plotting and printing their structure. After that, we know that the idecon
#    objects as well as decon[0-2] objects obtained by direct conversion are
#    correct.
# 5. Convert to decon from each typ to each other type and verify that the
#    objects obtained from these conversions are equal to the original
#    (verified) objects.

# Inputs #####

idecon <- deconvolute_spectrum(
    x       = sap$sap_01,
    sfr     = c(3.2, -3.2),
    smopts  = c(2, 3)
)
rdecon <- deconvolute_spectrum(
    x        = sap$sap_01,
    sfr      = c(3.2, -3.2),
    smopts   = c(2, 3),
    rtyp     = "rdecon",
    use_rust = TRUE
)
decon2i <- as_decon2(idecon) # (1)
decon1i <- as_decon1(idecon) # (1)
decon0i <- as_decon0(idecon) # (1)
decon2r <- as_decon2(rdecon) # (1)
decon1r <- try(as_decon1(rdecon), silent = TRUE) # (1)
decon0r <- try(as_decon0(rdecon), silent = TRUE) # (1)
# (1) Name objects as decon[0-2][ir] (with trailing i or r) to indicate the
# origin of these objects (idecon or rdecon).

# Tests #####

test_structure <- test_that("([ir]decon -> decon[0-2]) works", {
    expect_equal(names(idecon), idecon_members)
    expect_equal(names(decon2i), decon2_members)
    expect_equal(names(decon1i), decon1_members)
    expect_equal(names(decon0i), decon0_members)
    expect_equal(names(rdecon),  rdecon_members)
    expect_equal(names(decon2r), decon2_members)
    expect_equal(names(decon1r), NULL) # (1)
    expect_equal(names(decon0r), NULL) # (1)
    expect_equal(class(idecon), "idecon")
    expect_equal(class(decon2i), "decon2")
    expect_equal(class(decon1i), "decon1")
    expect_equal(class(decon0i), "list")
    expect_equal(class(rdecon), "rdecon")
    expect_equal(class(decon2r), "decon2")
    expect_equal(class(decon1r), "try-error") # (1)
    expect_equal(class(decon0r), "try-error") # (1)
    # (1) "Converting rdecon to decon[01] is not supported"
})

test_interactively <- if (identical(environment(), globalenv())) {
    str(idecon,  2, digits.d = 10)
    str(decon2i, 2, digits.d = 10)
    str(decon1i, 2, digits.d = 10)
    str(decon0i, 2, digits.d = 10)
    str(rdecon,  2, digits.d = 10)
    str(decon2r, 2, digits.d = 10)
    plot_spectrum(idecon,  sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    plot_spectrum(decon2i, sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    plot_spectrum(decon1i, sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    plot_spectrum(rdecon,  sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    plot_spectrum(decon2r, sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    # decon0 objects are not supported by plot_spectrum(), however, since decon0
    # is a strict subset of decon1, it's enough to plot decon1. Later on, the
    # decon0 and decon1 objects will be compared to verify that the decon0
    # object is valid as well.
}

test_decon22 <- test_that("(decon2 -> decon2) == (idecon -> decon2)", {
    decon22i <- as_decon2(decon2i)
    expect_equal(decon22i, decon2i)
})

test_decon21 <- test_that("(decon1 -> decon2) == (idecon -> decon2)", {
    decon21 <- as_decon2(decon1i)
    # Diffs in `meta$simpar`, `args`, `sit$wsrm`, `sit$nvrm` are expected, so we
    # patch them first to make the comparsion possible.
    decon21$sit$wsrm <- decon2i$sit$wsrm
    decon21$sit$nvrm <- decon2i$sit$nvrm
    decon21$meta$simpar <- decon2i$meta$simpar
    decon21$args <- decon2i$args
    expect_equal(decon21, decon2i)
})

test_decon20 <- test_that("(decon0 -> decon2) == (idecon -> decon2)", {
    decon20 <- as_decon2(decon0i, spectrum = sap[[1]])
    # Diffs in `meta$simpar`, `args`, `sit$wsrm`, `sit$nvrm` are expected, so we
    # patch them first to make the comparsion possible.
    decon20$sit$wsrm <- decon2i$sit$wsrm
    decon20$sit$nvrm <- decon2i$sit$nvrm
    decon20$meta$simpar <- decon2i$meta$simpar
    decon20$args <- decon2i$args
    expect_equal(decon20, decon2i)
})

test_decon12 <- test_that("(decon2 -> decon1) == (idecon -> decon1)", {
    decon12 <- as_decon1(decon2i);
    expect_equal(decon12, decon1i)
})

test_decon11 <- test_that("(decon1 -> decon1) == (idecon -> decon1)", {
    decon11 <- as_decon1(decon1i)
    expect_equal(decon11, decon1i)
})

test_decon10 <- test_that("(decon0 -> decon1) == (idecon -> decon1)", {
    decon10 <- as_decon1(decon0i, spectrum = sap[[1]])
    expect_equal(decon10, decon1i)
})

test_decon02 <- test_that("(decon2 -> decon0) == (idecon -> decon0)", {
    decon02 <- as_decon0(decon2i);
    expect_equal(decon02, decon0i)
})

test_decon01 <- test_that("(decon1 -> decon0) == (idecon -> decon0)", {
    decon01 <- as_decon0(decon1i)
    expect_equal(decon01, decon0i)
})

test_decon00 <- test_that("(decon0 -> decon0) == (idecon -> decon0)", {
    decon00 <- as_decon0(decon0i)
    expect_equal(decon00, decon0i)
})

test_reversibility <- test_that("conversion are reversible", {

    decon222 <- as_decon2(as_decon2(decon2i))
    decon212 <- as_decon2(as_decon1(decon2i))
    decon202 <- as_decon2(as_decon1(decon2i))
    decon121 <- as_decon1(as_decon2(decon1i))
    decon111 <- as_decon1(as_decon1(decon1i))
    decon101 <- as_decon1(as_decon0(decon1i), spectrum = sap[[1]])
    decon020 <- as_decon0(as_decon2(decon0i, spectrum = sap[[1]]))
    decon010 <- as_decon0(as_decon1(decon0i, spectrum = sap[[1]]))
    decon000 <- as_decon0(as_decon0(decon0i))

    # We know that conversion from decon2 to decon[01] is lossy, so when going
    # back from decon[01] to decon2, this information cannot be recovered
    # without additional user input. I.e., in order to allow the comparison of
    # decon212 and decon202 with decon2, we need to restore the missing
    # information manually.
    decon212$args        <- decon2i$args        # Optional element. Ok.
    decon212$meta$simpar <- decon2i$meta$simpar # Optional element. Ok.
    decon212$sit$wsrm    <- decon2i$sit$wsrm    # Required element. TODO.
    decon212$sit$nvrm    <- decon2i$sit$nvrm    # Required element. TODO.

    decon202$args        <- decon2i$args
    decon202$meta$simpar <- decon2i$meta$simpar
    decon202$sit$wsrm    <- decon2i$sit$wsrm
    decon202$sit$nvrm    <- decon2i$sit$nvrm

    expect_equal(decon020, decon0i)
    expect_equal(decon010, decon0i)
    expect_equal(decon000, decon0i)
    expect_equal(decon121, decon1i)
    expect_equal(decon111, decon1i)
    expect_equal(decon101, decon1i)
    expect_equal(decon222, decon2i)
    expect_equal(decon212, decon2i)
    expect_equal(decon202, decon2i)
})

