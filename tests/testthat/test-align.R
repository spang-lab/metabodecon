library(testthat)

# Prepare inputs
sap_01 <- sap[[1]]
sap_01_shifted <- simulate_spectrum(
    name = "sap_01_shifted",
    cs = sap_01$meta$simpar$cs,
    x0 = sap_01$meta$simpar$x0 + 0.3,
    A  = sap_01$meta$simpar$A,
    lambda = sap_01$meta$simpar$lambda,
    noise = sap_01$meta$simpar$noise
)
spectra <- as_spectra(list(sap_01, sap_01_shifted))
decons <- deconvolute(spectra, smopts=c(1,3), delta=3, sfr=c(3.2,-3.2), verbose = FALSE)

test_that("align works", {

    aligns <- align(decons, verbose = FALSE)

    # Check structure of returned object. Strategy: add all fields to the
    # decons object that we expect [align()] to add. At the end the objects
    # should be equal.
    decons_copy <- decons
    for (i in seq_along(aligns)) {
        decons_copy[[i]]$sit$al <- aligns[[i]]$sit$al
        decons_copy[[i]]$sit$supal <- aligns[[i]]$sit$supal
        decons_copy[[i]]$lcpar$x0_al <- aligns[[i]]$lcpar$x0_al
        class(decons_copy[[i]]) <- "align"
    }
    class(decons_copy) <- "aligns"
    expect_equal(object = aligns, expected = decons_copy)

    # Check that the alignment worked, our expectations are:
    # 1. x0_al     is shifted rougly 0.3 to the right compared to x0
    # 2. sit$al    is equal to the integral at peak-center-datapoints
    # 3. sit$al    is zero at non-peak-center-datapoints
    # 4. sit$supal is the superposition the aligned Lorentz Curves
    x0 <- aligns$sap_01_shifted$lcpar$x0
    x0_al <- aligns$sap_01_shifted$lcpar$x0_al
    shifts <- x0 - x0_al
    expect_true(all(shifts > 0.2 & shifts < 0.4))

    cs <- aligns$sap_01_shifted$cs
    dp <- seq_along(cs)
    pc <- match(x0_al, cs)
    non_pc <- setdiff(dp, pc)
    A <- aligns$sap_01_shifted$lcpar$A
    al <- aligns$sap_01_shifted$sit$al
    expect_equal(al[pc], A * pi)
    expect_equal(al[non_pc], rep(0, length(non_pc)))

    supal <- aligns$sap_01_shifted$sit$supal
    lambda <- aligns$sap_01_shifted$lcpar$lambda
    expect_equal(supal, lorentz_sup(cs, x0_al, A, lambda))
})

skip_if_slow_tests_disabled()

test_that("align can install its dependencies", {

    # TEST PURPOSE: If `metabodecon` or `speaq` is installed via
    # `install.packages()`, the dependencies `impute` and `MassSpecWavelet`
    # may not be installed. In this case, the missing dependencies should be
    # automatically detected and installed by [align()] when `install_deps` is TRUE.
    #
    # TEST STRATEGY: Remove all site libraries and unload the `speaq` dependencies.
    # Then, add a new (empty) site library and attempt to call `align()`. The
    # function should detect the missing dependencies in the new site library and
    # install them.
    #
    # IMPLEMENTATION DETAILS: Since `speaq` depends on `impute` and
    # `MassSpecWavelet`, we cannot unload these dependencies while `speaq` is
    # still loaded. Therefore, we need to unload `speaq` as well. However, after
    # removing the site libraries, R will no longer be able to find `speaq`. To
    # address this, we need to copy the `speaq` installation to the new site
    # library before calling `align()`. This ensures the expected scenario where
    # `speaq` is available but its dependencies are missing.

    # Load packages required by testthat before removing site libs
    requireNamespace("waldo", quietly = TRUE)   # required by `expect_false(<bool>)`
    requireNamespace("diffobj", quietly = TRUE) # required by `expect_false(<try-error>)`
    requireNamespace("glue", quietly = TRUE)    # required to print failure messages

    # Remember installation path of speaq
    speaq_path <- system.file(package = "speaq")

    # Remove existing site libs and instead add new (empty) site lib
    tmp_lib <- norm_path(tmpfile("tmp-library/"))
    dir.create(tmp_lib, recursive = TRUE)
    old_libs <- .libPaths()
    n_libs <- length(old_libs)
    defer(.libPaths(old_libs))
    .libPaths(c(tmp_lib, old_libs[n_libs]), include.site = FALSE)

    # Unload speaq, impute and MassSpecWavelet
    deps <- c("impute", "MassSpecWavelet")
    unloadNamespace("speaq")
    unloadNamespace(deps[1])
    unloadNamespace(deps[2])

    # Copy speaq to new site lib
    file.copy(speaq_path, tmp_lib, recursive = TRUE)

    deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
    expect_false(any(deps_installed))

    # Call align with auto-installation of dependencies enabled
    obj <- evalwith(
        output = "captured",
        message = "captured",
        expr = aligns <- align(decons, install_deps = TRUE, verbose = FALSE)
    )

    deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
    expect_true(all(deps_installed))
})
