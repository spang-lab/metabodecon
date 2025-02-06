deconvolute_rust <- function(spectrum,
                             sfr,
                             nfit = 3,
                             smopts = c(2, 5),
                             delta = 6.4,
                             ignore_regions = NULL,
                             parallel = FALSE,
                             optimize_settings = FALSE) {
    rust_spectrum <- convert_to_rust_spectrum(spectrum, sfr)
    deconvoluter <- set_up_deconvoluter(nfit, smopts, delta, ignore_regions)
    if (parallel) {
        # Optimize requires parallel feature of the Rust crate to be enabled, so
        # it should probably require it here as well.
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

multi_deconvolute_rust <- function(spectra,
                                   sfr,
                                   nfit = 3,
                                   smopts = c(2, 5),
                                   delta = 6.4,
                                   ignore_regions = NULL,
                                   parallel = FALSE,
                                   optimize_settings = FALSE) {
    rust_spectra <- lapply(spectra, function(spectrum) convert_to_rust_spectrum(spectrum, sfr))
    deconvoluter <- set_up_deconvoluter(nfit, smopts, delta, ignore_regions)
    if (parallel) {
        # Optimize requires parallel feature of the Rust crate to be enabled, so
        # it should probably require it here as well.
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

convert_to_rust_spectrum <- function(spectrum, sfr) {
    spectrum <- Spectrum$new(spectrum$cs, spectrum$si, sfr)

    spectrum
}

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
