library(testthat)

test_that("install_mdrb works", {

    # Remove all site libs to make sure mdrb is not pre-installed.
    # Load packages required for the test before doing this.
    requireNamespace("waldo", quietly = TRUE)   # required by `expect_false(<bool>)`
    requireNamespace("diffobj", quietly = TRUE) # required by `expect_false(<try-error>)`
    requireNamespace("glue", quietly = TRUE)    # required by testthat for test failure messages

    tmp_lib <- norm_path(tmpfile("tmp-library/"))
    dir.create(tmp_lib, recursive = TRUE)
    old_libs <- .libPaths()
    n_libs <- length(old_libs)
    defer(.libPaths(old_libs))
    .libPaths(tmp_lib, include.site = FALSE)
    unloadNamespace("mdrb")

    mdrb_available <- check_mdrb()
    if ((mdrb_available)) {
        skip("Skipping test as mdrb is installed in Sys Library and cannot be uninstalled.")
    }

    expect_false(mdrb_available)

    # Skip on Linux, to prevent source install, which takes up to 2 minutes.
    if (Sys.info()[["sysname"]] == "Linux") {
        skip_if_slow_tests_disabled()
    }

    obj <- evalwith(
        output = "captured",
        message = "captured",
        expr = {
            x <- try(install_mdrb(ask = FALSE))
            unlink("mdrb.out") # If executed by testthat
            unlink("tests/testthat/mdrb.out") # If executed interactively
        }
    )
    mdrb_available <- check_mdrb()
    if (getRversion() <  numeric_version("4.2")) {
        expect_true(inherits(x, "try-error"))
        expect_match(paste(obj$message, collapse = "\n"),
                     "requires R version 4.2 or greater")
        expect_false(mdrb_available)
    } else {
        # install_mdrb() no longer errors when no pre-built binary is available
        # (older R / unsupported platform): it prints guidance and returns
        # FALSE. When a binary IS available it installs mdrb and returns TRUE.
        # Either way the return value reflects mdrb's availability afterwards.
        expect_false(inherits(x, "try-error"))
        expect_type(x, "logical")
        expect_equal(x, mdrb_available)
    }
})
