library(testthat)

test_that("install_mdrb works", {

    # Remove all site libs from to make sure mdrb is not pre-installed.
    # Load packages required for the test before doing this.
    requireNamespace("waldo", quietly = TRUE)   # required by `expect_false(<bool>)`
    requireNamespace("diffobj", quietly = TRUE) # required by `expect_false(<try-error>)`
    tmp_lib <- norm_path(tmpfile("tmp-library/"))
    dir.create(tmp_lib, recursive = TRUE)
    old_libs <- .libPaths()
    n_libs <- length(old_libs)
    defer(.libPaths(old_libs))
    .libPaths(c(tmp_lib, old_libs[n_libs]))
    mdrb_available <- check_mdrb()
    expect_false(mdrb_available)

    # Try installation from binary
    obj <- evalwith(
        output = "captured",
        message = "captured",
        expr = x <- try(install_mdrb(ask = FALSE))
    )
    mdrb_available <- check_mdrb()
    if (getRversion() <  numeric_version("4.2")) {
        expect_true(inherits(x, "try-error"))
        expect_false(mdrb_available)
    } else {
        expect_true(is.null(x))
        expect_true(mdrb_available)
    }

    skip_if_slow_tests_disabled()

    # Test source installation
    if (mdrb_available) remove.packages("mdrb", lib = tmp_lib)
    obj <- evalwith(
        answers = "y",
        output = "captured",
        message = "captured",
        expr = x <- try(install_mdrb(type = "source", keep_outputs = TRUE))
    )
    mdrb_available <- check_mdrb()
    if (getRversion() <  numeric_version("4.2")) {
        expect_true(inherits(x, "try-error"))
        expect_false(mdrb_available)
    } else {
        expect_true(is.null(x))
        expect_true(mdrb_available)
    }
})
