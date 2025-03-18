library(testthat)

# Prepare Test Environment #####

# 1. Get the sap spectrum
# 2. Deconvolute the sap spectrum using deconvolute_spectrum
# 3. Convert the spectrum to decon[0-2]
# 4. Verify the values within decon[0-2] by plotting and printing their structure
# ==> we know that the decon[0-2] objects obtained by direct conversion are correct

idecon <- deconvolute_spectrum(sap[[1]], sfr = c(3.2, -3.2), smopts = c(2, 3))
decon2 <- as_decon2(idecon)
decon1 <- as_decon1(idecon)
decon0 <- as_decon0(idecon)

test_that("(idecon -> decon[0-2]) works", {
    expect_equal(names(idecon), idecon_members)
    expect_equal(names(decon2), decon2_members)
    expect_equal(names(decon1), decon1_members)
    expect_equal(names(decon0), decon0_members)
    expect_equal(class(idecon), "idecon")
    expect_equal(class(decon2), "decon2")
    expect_equal(class(decon1), "decon1")
    expect_equal(class(decon0), "list")
})

# In addition to above minimal checks, the decon objects should be verified
# manually by checking their values for plausibility using below plot and
# printouts.
if (identical(environment(), globalenv())) {
    str(idecon, 2, digits.d = 10)
    str(decon2, 2, digits.d = 10)
    str(decon1, 2, digits.d = 10)
    str(decon0, 2, digits.d = 10)
    plot_spectrum(idecon, sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
}

# Test conversions to decon2 ####

test_that("(decon2 -> decon2) == (idecon -> decon2)", {
    decon2_from2 <- as_decon2(decon2)
    expect_equal(decon2_from2, decon2)
})

test_that("(decon1 -> decon2) == (idecon -> decon2)", {

    decon21 <- as_decon2(decon1)

    # Diffs in `meta$simpar`, `args`, `sit$wsrm`, `sit$nvrm` are expected, so we
    # patch them first to make the comparsion possible.
    decon21$sit$wsrm <- decon2$sit$wsrm
    decon21$sit$nvrm <- decon2$sit$nvrm
    decon21$meta$simpar <- decon2$meta$simpar
    decon21$args <- decon2$args

    expect_equal(decon21, decon2)
})

test_that("(decon0 -> decon2) == (idecon -> decon2)", {

    decon20 <- as_decon2(decon0, spectrum = sap[[1]])

    # Diffs in `meta$simpar`, `args`, `sit$wsrm`, `sit$nvrm` are expected, so we
    # patch them first to make the comparsion possible.
    decon20$sit$wsrm <- decon2$sit$wsrm
    decon20$sit$nvrm <- decon2$sit$nvrm
    decon20$meta$simpar <- decon2$meta$simpar
    decon20$args <- decon2$args

    expect_equal(decon20, decon2)
})

# Test conversions to decon1 ####

test_that("(decon2 -> decon1) == (idecon -> decon1)", {
    decon12 <- as_decon1(decon2);
    expect_equal(decon12, decon1)
})

test_that("(decon1 -> decon1) == (idecon -> decon1)", {
    decon11 <- as_decon1(decon1)
    expect_equal(decon11, decon1)
})

test_that("(decon0 -> decon1) == (idecon -> decon1)", {
    decon10 <- as_decon1(decon0, spectrum = sap[[1]])
    expect_equal(decon10, decon1)
})

# Test conversions to decon0 ####

test_that("(decon2 -> decon0) == (idecon -> decon0)", {
    decon02 <- as_decon0(decon2);
    expect_equal(decon02, decon0)
})

test_that("(decon1 -> decon0) == (idecon -> decon0)", {
    decon01 <- as_decon0(decon1)
    expect_equal(decon01, decon0)
})

test_that("(decon0 -> decon0) == (idecon -> decon0)", {
    decon00 <- as_decon0(decon0)
    expect_equal(decon00, decon0)
})

# Test conversions are reversible ####

test_that("conversions are reversible", {

    decon020 <- as_decon0(as_decon2(decon0, spectrum = sap[[1]]))
    decon010 <- as_decon0(as_decon1(decon0, spectrum = sap[[1]]))
    decon000 <- as_decon0(as_decon0(decon0))
    decon121 <- as_decon1(as_decon2(decon1))
    decon111 <- as_decon1(as_decon1(decon1))
    decon101 <- as_decon1(as_decon0(decon1), spectrum = sap[[1]])
    decon222 <- as_decon2(as_decon2(decon2))
    decon212 <- as_decon2(as_decon1(decon2))
    decon202 <- as_decon2(as_decon1(decon2))

    # We know that conversion from decon2 to decon[01] is lossy, so when going
    # back from decon[01] to decon2, this information cannot be recovered by
    # default. I.e., in order to allow the comparison of decon212 and decon202
    # with decon2, we need to restore the missing information manually.
    decon212$sit$wsrm <- decon2$sit$wsrm
    decon212$sit$nvrm <- decon2$sit$nvrm
    decon212$meta$simpar <- decon2$meta$simpar
    decon212$args <- decon2$args
    decon202$sit$wsrm <- decon2$sit$wsrm
    decon202$sit$nvrm <- decon2$sit$nvrm
    decon202$meta$simpar <- decon2$meta$simpar
    decon202$args <- decon2$args

    expect_equal(decon020, decon0)
    expect_equal(decon010, decon0)
    expect_equal(decon000, decon0)
    expect_equal(decon121, decon1)
    expect_equal(decon111, decon1)
    expect_equal(decon101, decon1)
    expect_equal(decon222, decon2)
    expect_equal(decon212, decon2)
    expect_equal(decon202, decon2)
})
