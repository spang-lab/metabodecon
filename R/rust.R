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
    stopifnot(is_spectrum(sap[[1]]), is_num(sfr, 2))
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
