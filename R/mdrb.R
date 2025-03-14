#' @export
#'
#' @title Check Rust Backend Requirements
#'
#' @description
#' `check_mdrb()` returns a boolean indicating whether a suitable version of
#' the metabodecon Rust backend [mdrb](https://github.com/spang-lab/mdrb) is
#' currently installed.
#'
#' `check_mdrb_deps()` returns a list with detailed information about the
#' installation status of mdrb and its dependencies.
#'
#' @return
#' `check_mdrb()` returns TRUE if a suitable version of mdrb is installed,
#' else FALSE.
#'
#' `check_mdrb_deps()` returns a data.frame with the following columns:
#'
#' - `check`: description of the checked dependency as a string
#' - `passed`: boolean indicating whether the check passed
#' - `comment`: free text string describing the check result
#'
#' The rownames of the dataframe are: `r`, `rtools`, `cargo` and `rustc` and
#' correspond to the checked dependencies.
#'
#' Example:
#'
#' ```R
#' check_mdrb_deps()
#'
#' ##                check passed                                      comment
#' ## r           R >= 4.2   TRUE                             Current: R 4.4.2
#' ## rtools  Rtools exist   TRUE     Tested with: pkgbuild::has_build_tools()
#' ## cargo  cargo >= 1.78   TRUE Current: cargo 1.84.1 (66221abde 2024-11-19)
#' ## rustc  rustc >= 1.78   TRUE Current: rustc 1.84.1 (e71f9a9a9 2025-01-27)
#' ```
#'
#' @examples
#' check_mdrb()
#' check_mdrb_deps()
check_mdrb <- function() {
    tryCatch(
        packageVersion("mdrbx") >= package_version("0.0.1"),
        error = function(e) FALSE
    )
}

#' @export
#' @rdname check_mdrb
check_mdrb_deps <- function() {
    check <- c("R >= 4.2", "Rtools exist", "cargo >= 1.78", "rustc >= 1.78")
    df <- data.frame(check, passed = NA, comment = NA)
    rownames(df) <- c("r", "rtools", "cargo", "rustc")

    # Check R version
    r_cur <- package_version(getRversion())
    r_req <- package_version("4.2.0")
    df["r", "passed"] <- r_cur >= r_req
    df["r", "comment"] <- paste("Current: R", r_cur)

    # Check buildtools exist
    df["rtools", "passed"] <- pkgbuild::has_build_tools()
    df["rtools", "comment"] <- "Testcall: pkgbuild::has_build_tools()"

    # Add $HOME/.cargo/bin to PATH before checking cargo and rustc
    win <- .Platform$OS.type == "windows"
    home_cargo_bin <- paste(home(), ".cargo", "bin", sep = .Platform$file.sep)
    PATH <- Sys.getenv("PATH")
    PATH <- paste(PATH, home_cargo_bin, sep = if (win) ";" else ":")
    Sys.setenv(PATH = PATH)

    # Check cargo version
    cargo_out <- "not found in PATH"
    cargo_str <- "0.0.0"
    try(silent = TRUE, {
        local_options(warn = -1)
        cargo_out <- system("cargo --version", intern = TRUE, ignore.stderr = TRUE)
        cargo_str <- strsplit(cargo_out, " ")[[1]][2] # cargo_out == 'cargo x.y.z (commitsha yyyy-mm-dd)'
    })
    cargo_cur <- package_version(cargo_str)
    cargo_req <- package_version("1.78.0")
    df["cargo", "passed"] <- cargo_cur >= cargo_req
    df["cargo", "comment"] <- paste("Current:", cargo_out)

    # Check rustc version
    rustc_out <- "not found in PATH"
    rustc_str <- "0.0.0"
    try(silent = TRUE, {
        local_options(warn = -1)
        rustc_out <- system("rustc --version", intern = TRUE, ignore.stderr = TRUE)
        rustc_str <- strsplit(rustc_out, " ")[[1]][2] # rustc_out == 'rustc x.y.z (commitsha yyyy-mm-dd)'
    })
    rustc_cur <- package_version(rustc_str)
    rustc_req <- package_version("1.78.0")
    df["rustc", "passed"] <- rustc_cur >= rustc_req
    df["rustc", "comment"] <- paste("Current:", rustc_out)

    df
}


#' @export
#'
#' @title Install Rust Backend
#'
#' @description
#' Installs metabodecon's Rust backend [mdrb](https://github.com/spang-lab/mdrb)
#' via [remotes::install_github()]. If remotes is not available, it is installed
#' first via [utils::install.packages()]
#'
#' @param ask
#' Whether to ask for confirmation before attempting installion. Default is TRUE.
#'
#' @param args_remotes
#' Additional arguments to pass to [utils::install.packages()] when attemping
#' installation of [remotes](https://remotes.r-lib.org/index.html).
#'
#' @param args_mdrb
#' Additional arguments to pass to [remotes::install_github()] when attempting
#' installation of mdrb.
#'
#' @return NULL. Called for side effect of installing the Rust backend.
install_mdrb <- function(ask = TRUE,
                         args_remotes = list(),
                         args_mdrb = list()) {
    if (check_mdrb()) return()
    checks <- check_mdrb_deps()
    if (!all(checks$passed)) {
        tbl <- capture.output2(print(checks, row.names = FALSE, right = FALSE))
        msg <- paste(
            "cannot install rust backend due to missing system dependencies:\n",
            tbl, "",
            "Please install the missing dependencies and try again.",
            "For installation instructions, see the following links:\n",
            "R & Rtools:    https://cran.r-project.org/",
            "cargo & rustc: https://www.rust-lang.org/tools/install",
            sep = "\n"
        )
        stop(msg, call. = FALSE)
    }
    remotes_available <- requireNamespace("remotes", quietly = TRUE)
    pkg_vec <- if (remotes_available) "mdrb" else c("remotes", "mdrb")
    pkg_str <- paste(pkg_vec, collapse = " and ")
    pkg_word <- if (length(pkg_vec) == 1) "package" else "packages"
    msg <- sprintf("Proceeding will install the following %s: %s. Continue?", pkg_word, pkg_str)
    cont <- get_yn_input(msg)
    if (!cont) return()
    args_remotes$pkgs <- "remotes"
    if (!remotes_available) do.call(utils::install.packages, args_remotes)
    args_mdrb$repo <- "spang-lab/mdrb"
    do.call(remotes::install_github, args_mdrb)
    invisible()
}

#' @title Deconvolute Spectrum
#' @description Deconvolutes a given spectrum using the specified backend.
#' @param spectrum The spectrum to be deconvoluted.
#' @param sfr The spectral frequency range.
#' @param nfit Number of iterations for the analytical fitter. Default is 3.
#' @param smopts Smoothing options as a vector of two elements: iterations and window size. Default is c(2, 5).
#' @param delta Noise score threshold. Default is 6.4.
#' @param ignore_regions Regions to ignore during deconvolution. Default is NULL.
#' @param parallel Whether to perform deconvolution in parallel. Default is FALSE.
#' @param optimize_settings Whether to optimize settings. Default is FALSE.
#' @param backend The backend to use for deconvolution. Default is "R".
#' @return A list containing the deconvoluted spectrum and additional information.
deconvolute <- function(spectrum,
                        sfr,
                        nfit = 3,
                        smopts = c(2, 5),
                        delta = 6.4,
                        ignore_regions = NULL,
                        parallel = FALSE,
                        optimize_settings = FALSE,
                        backend = "R") {
    if (backend == "rust") {
        if (!check_rust_backend_requirements()) {
            install_mdrb(deps = NULL)
        }
        if (check_rust_backend_requirements()) {
            return(mdrb::deconvolute(spectrum, sfr, nfit, smopts, delta, ignore_regions, parallel, optimize_settings))
        } else {
            stop("Rust backend requirements are not met.")
        }
    }

    # Default R backend implementation
    rust_spectrum <- as_rust_spectrum(spectrum, sfr)
    stopifnot(is_spectrum(spectrum), is_num(sfr, 2))
    deconvoluter <- set_up_deconvoluter(nfit, smopts, delta, ignore_regions)
    if (parallel) {
        if (optimize_settings) {
            deconvoluter$optimize_settings(rust_spectrum)
        }
        deconvolution <- deconvoluter$par_deconvolute_spectrum(rust_spectrum)
    } else {
        deconvolution <- deconvoluter$deconvolute_spectrum(rust_spectrum)
    }
    decon2 <- convert_to_decon2(spectrum, deconvoluter, deconvolution, sfr)
    decon2$rust_spectrum <- rust_spectrum
    decon2$rust_deconvolution <- deconvolution
    decon2
}

#' @title Deconvolute Spectrum using Rust Backend
#' @description Deconvolutes a given spectrum using the Rust backend.
#' @param spectrum The spectrum to be deconvoluted.
#' @param sfr The spectral frequency range.
#' @param nfit Number of iterations for the analytical fitter. Default is 3.
#' @param smopts Smoothing options as a vector of two elements: iterations and window size. Default is c(2, 5).
#' @param delta Noise score threshold. Default is 6.4.
#' @param ignore_regions Regions to ignore during deconvolution. Default is NULL.
#' @param parallel Whether to perform deconvolution in parallel. Default is FALSE.
#' @param optimize_settings Whether to optimize settings. Default is FALSE.
#' @return A list containing the deconvoluted spectrum and additional information.
deconvolute_spectrum_rust <- function(spectrum,
                                      sfr,
                                      nfit = 3,
                                      smopts = c(2, 5),
                                      delta = 6.4,
                                      ignore_regions = NULL,
                                      parallel = FALSE,
                                      optimize_settings = FALSE) {
    rust_spectrum <- as_rust_spectrum(spectrum, sfr)
    stopifnot(is_spectrum(spectrum), is_num(sfr, 2))
    deconvoluter <- set_up_deconvoluter(nfit, smopts, delta, ignore_regions)
    if (parallel) {
        if (optimize_settings) {
            deconvoluter$optimize_settings(rust_spectrum)
        }
        deconvolution <- deconvoluter$par_deconvolute_spectrum(rust_spectrum)
    } else {
        deconvolution <- deconvoluter$deconvolute_spectrum(rust_spectrum)
    }
    decon2 <- convert_to_decon2(spectrum, deconvoluter, deconvolution, sfr)
    decon2$rust_spectrum <- rust_spectrum
    decon2$rust_deconvolution <- deconvolution
    decon2
}

#' @title Multi Deconvolute Spectra using Rust Backend
#' @description Deconvolutes a list of spectra using the Rust backend.
#' @param spectra A list of spectra to be deconvoluted.
#' @param sfr The spectral frequency range.
#' @param nfit Number of iterations for the analytical fitter. Default is 3.
#' @param smopts Smoothing options as a vector of two elements: iterations and window size. Default is c(2, 5).
#' @param delta Noise score threshold. Default is 6.4.
#' @param ignore_regions Regions to ignore during deconvolution. Default is NULL.
#' @param parallel Whether to perform deconvolution in parallel. Default is FALSE.
#' @param optimize_settings Whether to optimize settings. Default is FALSE.
#' @return A list of deconvoluted spectra and additional information.
deconvolute_spectra_rust <- function(spectra,
                                     sfr,
                                     nfit = 3,
                                     smopts = c(2, 5),
                                     delta = 6.4,
                                     ignore_regions = NULL,
                                     parallel = FALSE,
                                     optimize_settings = FALSE) {
    rust_spectra <- lapply(spectra, function(spectrum) as_rust_spectrum(spectrum, sfr))
    stopifnot(is_spectrum(spectrum), is_num(sfr, 2))
    deconvoluter <- set_up_deconvoluter(nfit, smopts, delta, ignore_regions)
    if (parallel) {
        if (optimize_settings) {
            deconvoluter$optimize_settings(rust_spectra[[1]])
        }
        deconvolutions <- deconvoluter$par_deconvolute_spectra(rust_spectra)
    } else {
        deconvolutions <- deconvoluter$deconvolute_spectra(rust_spectra)
    }
    decon2s <- lapply(seq_along(spectra), function(i) convert_to_decon2(spectra[[i]], deconvoluter, deconvolutions[[i]], sfr))
    for (i in seq_along(decon2s)) {
        decon2s[[i]]$rust_spectrum <- rust_spectra[[i]]
        decon2s[[i]]$rust_deconvolution <- deconvolutions[[i]]
    }
    decon2s
}

#' @title Convert to Rust Spectrum
#' @description Converts a given spectrum to a Rust spectrum.
#' @param spectrum The spectrum to be converted.
#' @param sfr The spectral frequency range.
#' @return A Rust spectrum object.
as_rust_spectrum <- function(spectrum, sfr) {
    stopifnot(is_spectrum(spectrum), is_num(sfr, 2))
    spectrum <- Spectrum$new(spectrum$cs, spectrum$si, sfr)
    spectrum
}

#' @title Convert to Decon2
#' @description Converts a given spectrum and deconvolution results to Decon2 format.
#' @param spectrum The original spectrum.
#' @param deconvoluter The deconvoluter object.
#' @param deconvolution The deconvolution results.
#' @param sfr The spectral frequency range.
#' @return A Decon2 object.
convert_to_decon2 <- function(spectrum, deconvoluter, deconvolution, sfr) {
    lorentzians <- deconvolution$lorentzians()
    spectrum$args <- list(
        nfit = deconvoluter$fitting_settings()$iterations,
        smopts = c(deconvoluter$smoothing_settings()$iterations, deconvoluter$smoothing_settings()$window_size),
        delta = deconvoluter$selection_settings()$threshold,
        sfr = sfr,
        ignore_regions = deconvoluter$ignore_regions()
    )
    spectrum$sit <- list(
        wsrm = NA,
        nvrm = NA,
        sm = NA,
        sup = deconvolution$superposition_vec(spectrum$cs)
    )
    spectrum$peak <- NA
    spectrum$lcpar <- data.frame(
        A = lorentzians$A,
        lambda = lorentzians$lambda,
        x0 = lorentzians$x0
    )
    spectrum$mse <- list(
        raw = deconvolution$mse(),
        norm = deconvolution$mse() / sum(spectrum$sit$sup),
        sm = NA,
        smnorm = NA
    )
    spectrum
}

#' @title Set Up Deconvoluter
#' @description Sets up the deconvoluter with specified parameters.
#' @param nfit Number of iterations for the analytical fitter.
#' @param smopts Smoothing options as a vector of two elements: iterations and window size.
#' @param delta Noise score threshold.
#' @param ignore_regions Regions to ignore during deconvolution.
#' @return A Deconvoluter object.
set_up_deconvoluter <- function(nfit, smopts, delta, ignore_regions) {
    if (length(ignore_regions) %% 2 != 0) {
        stop("ignore_regions must have an even number of elements")
    }
    deconvoluter <- Deconvoluter$new()
    deconvoluter$set_moving_average_smoother(smopts[1], smopts[2])
    deconvoluter$set_noise_score_selector(delta)
    deconvoluter$set_analytical_fitter(nfit)
    if (length(ignore_regions) > 0) {
        for (i in seq(1, length(ignore_regions), by = 2)) {
            deconvoluter$add_ignore_region(ignore_regions[i], ignore_regions[i + 1])
        }
    }
    deconvoluter
}
