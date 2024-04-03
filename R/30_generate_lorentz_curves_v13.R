#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spetra and generates a Lorentz curve for each detected signal within a spectra.
#' @param data_path Either the path to an existing directory containing measured NMR spectra or a dataframe with columns `ppm` (parts per million) and `si` (signal intensity) or a list of such dataframes.
#' @param file_format Format of the spectra files. Either `"bruker"` or `"jcampdx"`. Only relevant if `data_path` is a directory.
#' @param make_rds Store results as rds file on disk? Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param expno The experiment number for the spectra files. E.g. `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param procno The processing number for the spectra. E.g. `"10"`. Only relevant if `data_path` is a directory and `file_format` is `"bruker"`.
#' @param nfit Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param wshw Half width of the water artefact in ppm.
#' @param sfr Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smopts Vector with two entries consisting of the number of smoothing iterations and the number of data points to use for smoothing (must be uneven). TODO: add details.
#' @param delta Threshold value to distinguish between signal and noise. TODO: add details.
#' @param sf Vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask  Whether to ask for user input during the deconvolution process. If set to FALSE, the provided default values will be used.
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum. TODO: add details.
#' @examples \dontrun{
#' xds_path <- download_example_datasets()
#' data_path <- file.path(xds_path, "bruker/urine")
#' file_format <- "bruker"
#' ask <- FALSE
#' nfit <- 2
#' generate_lorentz_curves_v12(data_path, file_format, ask = ask, nfit = nfit)
#' }
#' @noRd
generate_lorentz_curves_v13 <- function(data_path = file.path(download_example_datasets(), "bruker/urine"),
                                        file_format = "bruker",
                                        make_rds = FALSE,
                                        expno = 10,
                                        procno = 10,
                                        nfit = 10,
                                        wshw = 0.1527692,
                                        sfr = c(11.44494, -1.8828),
                                        smopts = c(2, 5),
                                        delta = 6.4,
                                        sf = c(1000, 1000000),
                                        ask = TRUE,
                                        debug = FALSE,
                                        ncores = 2) {

    # Read spectra and ask user for parameters
    spectra <- read_spectra(data_path, file_format, expno, procno, ask, sf)
    adjno <- get_adjno(spectra, sfr, wshw, ask)
    spectra <- add_sfrs(spectra, sfr, ask, adjno)
    spectra <- add_wsrs(spectra, wshw, ask, adjno)

    # Deconvolute spectra
    n <- length(spectra)
    nams <- names(spectra)
    spectra <- lapply(seq_len(n), function(i) {
        nam <- nams[i]
        msg("Starting deconvolution of", nam)
        spec <- spectra[[i]]
        spec <- rm_water_signal_v12(spec)
        spec <- rm_negative_signals_v12(spec)
        spec <- smooth_signals_v12(spec, reps = smopts[1], k = smopts[2])
        spec <- find_peaks_v12(spec)
        spec <- rm_peaks_with_low_scores_v12(spec, delta)
        spec <- fit_lorentz_curves(spec, nfit)
        spec <- add_return_list_v13(spec, n, nam, debug)
    })
    names(spectra) <- nams

    # Create return lists
    if (debug) {
        return(spectra)
    } else {
        rets <- lapply(spectra, function(spec) spec$ret)
        names(rets) <- nams
        return(rets)
    }
}

#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param n Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @param debug Add debug info to the return list
add_return_list_v13 <- function(spec = glc_urine1_yy_ni3_dbg()$rv$urine_1,
                                n = 1,
                                nam = "urine_1",
                                debug = TRUE) {
    spec$ret <- list(
        number_of_files = n,
        filename = nam,
        x_values = spec$sdp,
        x_values_ppm = spec$ppm,
        y_values = spec$y_smooth,
        spectrum_superposition = spec$lcr$spectrum_approx[1, ],
        mse_normed = spec$lcr$mse_normed,
        index_peak_triplets_middle = spec$peak$center[spec$peak$high],
        index_peak_triplets_left = spec$peak$right[spec$peak$high],
        index_peak_triplets_right = spec$peak$left[spec$peak$high],
        peak_triplets_middle = spec$ppm[spec$peak$center[spec$peak$high]],
        peak_triplets_left = spec$ppm[spec$peak$right[spec$peak$high]],
        peak_triplets_right = spec$ppm[spec$peak$left[spec$peak$high]],
        integrals = spec$lcr$integrals,
        signal_free_region = c(spec$sfr$left_sdp, spec$sfr$right_sdp),
        range_water_signal_ppm = spec$wsr$hwidth_ppm,
        A = spec$lcr$A,
        lambda = spec$lcr$lambda,
        x_0 = spec$lcr$w
    )
    spec
}