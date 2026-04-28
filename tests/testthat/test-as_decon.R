library(testthat)

# Test Strategy #####

# 1. Get the sap spectrum.
# 2. Deconvolute the sap spectrum using [deconvolute_spectrum()] into a decon2
#    object.
# 3. Convert the decon2 object to decon0, decon1 objects.
# 4. Verify the structure of all decon[0-2] objects.
# 5. Convert each type to each other type and verify equality.

# Inputs #####
trySilent <- function(expr) try(expr, silent = TRUE)
decon2 <- trySilent(deconvolute_spectrum(
    x = sap$sap_01,
    sfr = c(3.2, -3.2),
    smopts = c(2, 3),
    verbose = FALSE
))
rdecon <- trySilent(deconvolute_spectrum(
    x = sap$sap_01,
    sfr = c(3.2, -3.2),
    smopts = c(2, 3),
    rtyp = "rdecon",
    use_rust = TRUE,
    verbose = FALSE
))
decon12 <- try(as_decon1(decon2), silent = TRUE) # decon1 from decon2
decon02 <- try(as_decon0(decon2), silent = TRUE) # decon0 from decon2
decon2r <- try(as_decon2(rdecon), silent = TRUE) # decon2 from rdecon
decon1r <- try(as_decon1(rdecon), silent = TRUE) # expected try-error
decon0r <- try(as_decon0(rdecon), silent = TRUE) # expected try-error
mdrb_available <- check_mdrb()

# Checks #####

test_derivatives_decon2 <- test_that("decon2 -> decon[0-2] works", {
    withr::local_output_sink(nullfile())
    expect_equal(names(decon2),  decon2_members)
    expect_equal(names(decon12), decon1_members)
    expect_equal(names(decon02), decon0_members)
    expect_equal(class(decon2),  "decon2")
    expect_equal(class(decon12), "decon1")
    expect_equal(class(decon02), "list")
})

test_interactively_decon2 <- if (identical(environment(), globalenv())) {
    str(decon2,  2, digits.d = 10)
    str(decon12, 2, digits.d = 10)
    str(decon02, 2, digits.d = 10)
    plot_spectrum(decon2,  sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    plot_spectrum(decon12, sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
}

test_decon22 <- test_that("(decon2 <- decon2) is a roundtrip", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon2(decon2), decon2)
})

test_decon21 <- test_that("(decon2 <- decon1) == (decon2 <- decon2)", {
    withr::local_output_sink(nullfile())
    x <- as_decon2(decon12)
    # Diffs in `meta$simpar`, `args`, `sit$wsrm`, `sit$nvrm` are expected, so we
    # patch them first to make the comparison possible. Row names of `peak` are
    # also reset during roundtrip, so patch those too.
    x$sit$wsrm    <- decon2$sit$wsrm
    x$sit$nvrm    <- decon2$sit$nvrm
    x$meta$simpar <- decon2$meta$simpar
    x$args        <- decon2$args
    row.names(x$peak) <- as.integer(row.names(decon2$peak))
    expect_equal(x, decon2)
})

test_decon20 <- test_that("(decon2 <- decon0) == (decon2 <- decon2)", {
    withr::local_output_sink(nullfile())
    x <- as_decon2(decon02, spectrum = sap[[1]])
    # Diffs in `meta$simpar`, `args`, `sit$wsrm`, `sit$nvrm` are expected, so we
    # patch them first to make the comparison possible. Row names of `peak` are
    # also reset during roundtrip, so patch those too.
    x$sit$wsrm    <- decon2$sit$wsrm
    x$sit$nvrm    <- decon2$sit$nvrm
    x$meta$simpar <- decon2$meta$simpar
    x$args        <- decon2$args
    row.names(x$peak) <- as.integer(row.names(decon2$peak))
    expect_equal(x, decon2)
})

test_decon12 <- test_that("(decon1 <- decon2) == (decon1 <- decon2)", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon1(decon2), decon12)
})

test_decon11 <- test_that("(decon1 <- decon1) is a roundtrip", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon1(decon12), decon12)
})

test_decon10 <- test_that("(decon1 <- decon0) == (decon1 <- decon2)", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon1(decon02, spectrum = sap[[1]]), decon12)
})

test_decon02 <- test_that("(decon0 <- decon2) == (decon0 <- decon2)", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon0(decon2), decon02)
})

test_decon01 <- test_that("(decon0 <- decon1) == (decon0 <- decon2)", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon0(decon12), decon02)
})

test_decon00 <- test_that("(decon0 <- decon0) is a roundtrip", {
    withr::local_output_sink(nullfile())
    expect_equal(as_decon0(decon02), decon02)
})

test_decon2_reversibility <- test_that("conversions are reversible", {
    withr::local_output_sink(nullfile())

    decon222 <- as_decon2(as_decon2(decon2))
    decon212 <- as_decon2(as_decon1(decon2))
    decon202 <- as_decon2(as_decon1(decon2))
    decon121 <- as_decon1(as_decon2(decon12))
    decon111 <- as_decon1(as_decon1(decon12))
    decon101 <- as_decon1(as_decon0(decon12), spectrum = sap[[1]])
    decon020 <- as_decon0(as_decon2(decon02, spectrum = sap[[1]]))
    decon010 <- as_decon0(as_decon1(decon02, spectrum = sap[[1]]))
    decon000 <- as_decon0(as_decon0(decon02))

    # We know that conversion from decon2 to decon[01] is lossy, so when going
    # back from decon[01] to decon2, this information cannot be recovered
    # without additional user input. I.e., in order to allow the comparison of
    # decon212 and decon202 with decon2, we need to restore the missing
    # information manually.
    decon212$args        <- decon2$args        # Optional element. Ok.
    decon212$meta$simpar <- decon2$meta$simpar # Optional element. Ok.
    decon212$sit$wsrm    <- decon2$sit$wsrm    # Required element. TODO.
    decon212$sit$nvrm    <- decon2$sit$nvrm    # Required element. TODO.
    row.names(decon212$peak) <- as.integer(row.names(decon2$peak)) # Reset by roundtrip.

    decon202$args        <- decon2$args
    decon202$meta$simpar <- decon2$meta$simpar
    decon202$sit$wsrm    <- decon2$sit$wsrm
    decon202$sit$nvrm    <- decon2$sit$nvrm
    row.names(decon202$peak) <- as.integer(row.names(decon2$peak)) # Reset by roundtrip.

    expect_equal(decon020, decon02)
    expect_equal(decon010, decon02)
    expect_equal(decon000, decon02)
    expect_equal(decon121, decon12)
    expect_equal(decon111, decon12)
    expect_equal(decon101, decon12)
    expect_equal(decon222, decon2)
    expect_equal(decon212, decon2)
    expect_equal(decon202, decon2)
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

test_derivatives_rdecon <- test_that("rdecon -> decon[0-2] works", {
    withr::local_output_sink(nullfile())
    expect_equal(names(rdecon),  rdecon_members)
    expect_equal(names(decon2r), decon2_members)
    expect_equal(names(decon1r), NULL) # (1)
    expect_equal(names(decon0r), NULL) # (1)
    expect_equal(class(rdecon),  "rdecon")
    expect_equal(class(decon2r), "decon2")
    expect_equal(class(decon1r), "try-error") # (1)
    expect_equal(class(decon0r), "try-error") # (1)
    # (1) "Converting rdecon to decon[01] is not supported"
})

test_interactively_rdecon <- if (identical(environment(), globalenv())) {
    str(rdecon,  2, digits.d = 10)
    str(decon2r, 2, digits.d = 10)
    plot_spectrum(rdecon,  sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
    plot_spectrum(decon2r, sub1 = list(lt_axis = list(sf = 100)), sub2 = TRUE)
}
