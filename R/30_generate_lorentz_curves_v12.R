# Private API Functions #####

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
generate_lorentz_curves_v12 <- function(data_path = file.path(download_example_datasets(), "bruker/urine"),
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
        lc2 <- init_lorentz_curves_v13(x = spec$sdp, y = spec$y_smooth, pc = spec$peak$center[spec$peak$high], pl = spec$peak$right[spec$peak$high], pr = spec$peak$left[spec$peak$high])
        spec <- init_lorentz_curves_v12(spec)
        lc1 <- spec$lc
        spec <- refine_lorentz_curves_v12(spec, nfit)
        lcr1 <- spec$lcr
        spec <- add_return_list_v12(spec, n, nam, debug)
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

# Helpers for generate_lorentz_curves_v12 #####

#' @title Get the spectrum number for adjusting parameters
#' @description Asks the user which spectrum should be used for adjusting signal free region (SFR) and water signal half width (WSHW) and returns the corresponding spectrum number. If `ask` is `FALSE` or the user chooses to use different parameters for each spectrum, 0 is returned.
#' @param spectra A list of spectra.
#' @param sfr A numeric value or a list specifying the SFR for each spectrum.
#' @param wshw A numeric value or a list specifying the WSHW for each spectrum.
#' @param ask A logical value indicating whether to ask the user to confirm the spectrum number. Default is TRUE.
#' @return The spectrum number chosen by the user for adjusting parameters, or 0 if `ask` is `FALSE` or the user chooses to use different parameters for each spectrum.
#' @noRd
get_adjno <- function(spectra, sfr, wshw, ask) {
    if (!ask || length(spectra) == 1) {
        return(0)
    }
    same_param <- get_yn_input("Use same parameters for all spectra?")
    if (!same_param) {
        return(0)
    }
    namestr <- paste(seq_along(spectra), names(spectra), sep = ": ", collapse = ", ")
    prompt <- sprintf("Number of spectrum for adjusting parameters? (%s)", namestr)
    get_num_input(prompt, min = 1, max = length(spectra), int = TRUE)
}

#' @title Add signal free region to each spectrum
#' @description Adds details about each spectrum's signal free region (SFR) to the respective spectrum in the list.
#' @param spectra A list of spectra.
#' @param sfr A numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm, or a list of such vectors (one for each spectrum).
#' @param ask A logical value indicating whether to ask the user to confirm the SFR(s).
#' @param adjno The spectrum number for adjusting the SFR. If 0, the user will be asked to confirm/adjust the SFR for each spectrum. If > 0, the SFR of the spectrum with this number will be used for all spectra.
#' @return A list of spectra with the SFR added to each spectrum.
#' @noRd
add_sfrs <- function(spectra, sfr, ask, adjno) {
    n <- length(spectra)

    # Check initial, user provided SFR(s). Convert to list if necessary.
    sfrs <- if (is_num(sfr, 2)) {
        rep(list(sfr), n)
    } else if (is_list_of_nums(sfr, n, 2)) {
        sfr
    } else {
        stop(sprintf("Argument `sfr` should be either\n- a numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm or\n- a list of length %d of such vectors (one for each spectrum)", n))
    }

    # Ask user to confirm/adjust initial SFR(s).
    if (ask) {
        if (adjno == 0) {
            sfrs <- lapply(seq_len(n), function(i) confirm_sfr(spectra[[i]], sfrs[[i]], ask))
        } else if (adjno > 0) {
            sfr_adjno <- confirm_sfr(spectra[[adjno]], sfrs[[adjno]], ask)
            sfrs <- rep(list(sfr_adjno), n)
        }
    }

    # Add SFR in different units to each spectrum.
    for (i in seq_len(n)) {
        spectra[[i]] <- add_sfr(spectra[[i]], sfrs[[i]])
    }
    spectra
}

#' @title Add water signal region to each spectrum
#' @description Adds details about the water signal region (WSR) to each spectrum in the list.
#' @param spectra A list of spectra.
#' @param wshw A numeric value giving the half width of the water artefact in ppm, or a vector of such values (one for each spectrum).
#' @param ask A logical value indicating whether to ask the user to confirm the WSHW(s).
#' @param adjno The spectrum number for adjusting the WSHW. If 0, the user will be asked to confirm/adjust the WSHW for each spectrum. If > 0, the WSHW of the spectrum with this number will be used for all spectra.
#' @return A list of spectra with the WSR added to each spectrum.
#' @noRd
add_wsrs <- function(spectra, wshw, ask, adjno) {
    n <- length(spectra)

    # Check initial, user provided WSHW(s). Convert to list if necessary.
    if (is_num(wshw, 1)) {
        wshws <- rep(list(wshw), n)
    } else if (is_num(wshw, n) || is_list_of_nums(wshw, n, 1)) {
        wshws <- as.list(wshw)
    } else {
        stop(sprintf("Argument `wshw` should be either\n- a single value, giving the half width of the water artefact in ppm or\n- a vector of length %d of such values (one for each spectrum)", n))
    }

    # Ask user to confirm/adjust initial WSHW(s).
    if (ask) {
        if (adjno == 0) {
            wshws <- lapply(seq_len(n), function(i) confirm_wshw(spectra[[i]], wshws[[i]], ask))
        } else if (adjno > 0) {
            wshw_adjno <- confirm_wshw(spectra[[adjno]], wshws[[adjno]], ask)
            wshws <- rep(list(wshw_adjno), n)
        }
    }

    # Calculate WSR from WSHW for each spectrum.
    for (i in seq_len(n)) {
        spectra[[i]] <- add_wsr(spectra[[i]], wshws[[i]])
    }
    spectra
}

rm_water_signal_v12 <- function(spec) {
    msg("Removing water signal")
    y <- spec$y_scaled
    left <- spec$wsr$left_dp
    right <- spec$wsr$right_dp
    y[right:left] <- 0.01 / spec$sfy # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    spec$y_nows <- y
    spec
}

rm_negative_signals_v12 <- function(spec) {
    msg("Removing negative signals")
    if (is.null(spec$y_nows)) stop("Water signal not removed yet. Please call `rm_water_signal_v12()` first.")
    spec$y_pos <- abs(spec$y_nows)
    spec
}

#' @title Smooth signal intensities using a moving average
#' @description This function smooths signal intensities by applying a [moving average](https://en.wikipedia.org/wiki/Moving_average) filter with a window size of k.
#' @param spec A list representing the spectrum, which should include the scaled signal intensities, after removal of the water artefact and negative values (`spec$y_pos`).
#' @param reps The number of times to apply the moving average.
#' @param k The number of points within the moving average window. Must be odd, so the smoothed point is in the middle of the window.
#' @return A numeric vector of the smoothed values.
#' @details Old and slow version producing the same results as the implementation within `deconvolution` from `MetaboDecon1D_deconvolution.R`.
#' @noRd
smooth_signals_v12 <- function(spec, reps = 2, k = 5) {
    msg("Smoothing signals")
    if (k %% 2 == 0) stop("k must be odd")
    Z <- vector("list", length = reps)
    y <- spec$y_pos
    n <- length(y)
    for (i in 1:reps) {
        z <- y
        for (j in 1:(n)) {
            left_border <- j - floor(k / 2)
            right_border <- j + floor(k / 2)
            if (left_border <= 0) {
                left_border <- 1
                z[j] <- (1 / right_border) * sum(y[left_border:right_border])
            } else if (right_border >= n) {
                right_border <- n
                z[j] <- (1 / (right_border - left_border + 1)) * sum(y[left_border:right_border])
            } else {
                z[j] <- (1 / k) * sum(y[left_border:right_border])
            }
        }
        y <- Z[[i]] <- as.numeric(z)
    }
    spec$Z <- Z
    spec$y_smooth <- Z[[reps]]
    spec
}

#' @inherit find_peaks_v12
#' @param details Successor of `select_peaks_v0`, `find_left_positions_v0` and `find_right_positions_v0`. The new function `find_peaks_v12()` is a combination of these three functions with a corrected naming convention: what was incorrectly referred to as "left" is now correctly called "right" and vice versa.
#' @noRd
find_peaks_v12 <- function(spec) {
    msg("Starting peak selection")
    d <- spec$d <- calc_second_derivative_v12(y = spec$y_smooth)
    a <- abs(d)
    m <- length(d)
    dl <- c(NA, d[-m]) # dl[i] == d[i-1]
    dr <- c(d[-1], NA) # dr[i] == d[i+1]
    center <- which(d < 0 & d <= dl & d < dr)
    spec$peak <- data.frame(left = NA, center = center, right = NA, score = NA)
    for (i in seq_along(center)) {
        j <- center[i]
        l <- spec$peak$left[i] <- get_left_border_v12(j, d)
        r <- spec$peak$right[i] <- get_right_border_v12(j, d, m)
        spec$peak$score[i] <- get_peak_score_v12(j, l, r, a)
    }
    msg("Detected", length(center), "peaks")
    return(spec)
}

rm_peaks_with_low_scores_v13 <- function(ppm, # x values in ppm
                                         pc, # peak center indices
                                         ps, # peak scores
                                         sfrl, # signal free region left in ppm
                                         sfrr, # signal free region right in ppm
                                         delta = 6.4) { # threshold parameter to distinguish between "real" peaks from noise
    if (any(is.na(ps))) stop("Peak scores must never be NA")
    msg("Removing peaks with low scores")
    in_sfr <- which(ppm[pc] >= ppm || ppm[pc] <= ppm)
    mu <- mean(ps[in_sfr]) # mean (greek letter mu)
    sigma <- sd(ps[in_sfr]) # standard deviation (greek letter sigma)
    gt_tau <- (ps >= mu + delta * sigma) # greater than threshold value (greek letter tau)
    msg("Removed", sum(!gt_tau), "peaks")
    list(in_sfr, gt_tau) # gttho
}

rm_peaks_with_low_scores_v12 <- function(spec, delta = 6.4) {
    msg("Removing peaks with low scores")
    score <- spec$peak$score
    l <- which(spec$sdp[spec$peak$center] >= spec$sfr$left_sdp)
    r <- which(spec$sdp[spec$peak$center] <= spec$sfr$right_sdp)
    mu <- mean(score[c(l, r)])
    sigma <- sd(c(score[l], score[r]))
    spec$peak$high <- score >= mu + delta * sigma
    spec$peak$region <- "norm"
    spec$peak$region[l] <- "sfrl"
    spec$peak$region[r] <- "sfrr"
    msg("Removed", sum(!spec$peak$high), "peaks")
    spec
}

refine_lorentz_curves_v12 <- function(spec, nfit) {
    msg("Refining Lorentz curves")

    spectrum_x <- spec$sdp
    spectrum_y <- spec$y_smooth
    filtered_peaks <- as.integer(spec$peak$center[spec$peak$high] - 1)
    filtered_left_position <- spec$peak$right[spec$peak$high] - 1
    filtered_right_position <- spec$peak$left[spec$peak$high] - 1
    A <- spec$lc$A
    lambda <- spec$lc$lambda
    w <- spec$lc$w
    nfit <- nfit

    # Calculate all initial lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spectrum_x))
    for (i in seq_along(filtered_peaks)) {
        # If A = 0, then the lorentz curve is a zero line
        if (A[i] == 0) {
            lorentz_curves_initial[i, ] <- 0
        } else {
            lorentz_curves_initial[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (spectrum_x - w[i])^2)))
        }
    }

    # Approximation of lorentz curves
    for (b in 1:nfit) {
        # Calculate new heights of peak triplets
        w_1_new <- c()
        w_2_new <- c()
        w_3_new <- c()
        y_1_new <- c()
        y_2_new <- c()
        y_3_new <- c()
        w_1_2_new <- c()
        w_1_3_new <- c()
        w_2_3_new <- c()
        y_1_2_new <- c()
        y_1_3_new <- c()
        y_2_3_new <- c()
        w_delta_new <- c()
        w_new <- c()
        lambda_new <- c()
        A_new <- c()
        sum_left <- c()
        sum_peaks <- c()
        sum_right <- c()
        proportion_left <- c()
        proportion_peaks <- c()
        proportion_right <- c()

        for (i in seq_along(filtered_peaks)) {
            # Calculate the position of the peak triplets
            w_1_new <- c(w_1_new, spectrum_x[filtered_left_position[i] + 1])
            w_2_new <- c(w_2_new, spectrum_x[filtered_peaks[i] + 1])
            w_3_new <- c(w_3_new, spectrum_x[filtered_right_position[i] + 1])

            # Calculate the sum of all lorentz curves for each data point
            sum_left[i] <- sum(lorentz_curves_initial[seq_along(filtered_left_position), filtered_left_position[i] + 1])
            sum_peaks[i] <- sum(lorentz_curves_initial[seq_along(filtered_peaks), filtered_peaks[i] + 1])
            sum_right[i] <- sum(lorentz_curves_initial[seq_along(filtered_right_position), filtered_right_position[i] + 1])

            # Calculate the proprotion between original spectrum an the sum of the lorentz curves for each peak triplets position
            proportion_left[i] <- spectrum_y[filtered_left_position[i] + 1] / sum_left[i]
            proportion_peaks[i] <- spectrum_y[filtered_peaks[i] + 1] / sum_peaks[i]
            proportion_right[i] <- spectrum_y[filtered_right_position[i] + 1] / sum_right[i]

            # Calculate the new heights of the peak triplets
            y_1_new[i] <- lorentz_curves_initial[i, filtered_left_position[i] + 1] * proportion_left[i]
            y_2_new[i] <- lorentz_curves_initial[i, filtered_peaks[i] + 1] * proportion_peaks[i]
            y_3_new[i] <- lorentz_curves_initial[i, filtered_right_position[i] + 1] * proportion_right[i]

            # Calculate mirrored points if necesccary
            # For ascending shoulders
            if ((y_1_new[i] < y_2_new[i]) && (y_2_new[i] < y_3_new[i])) {
                w_3_new[i] <- 2 * w_2_new[i] - w_1_new[i]
                y_3_new[i] <- y_1_new[i]
            }
            # For descending shoulders
            if ((y_1_new[i] > y_2_new[i]) && (y_2_new[i] > y_3_new[i])) {
                w_1_new[i] <- 2 * w_2_new[i] - w_3_new[i]
                y_1_new[i] <- y_3_new[i]
            }

            # Move triplet to zero position
            w_delta_new[i] <- w_1_new[i]
            w_1_new[i] <- w_1_new[i] - w_delta_new[i]
            w_2_new[i] <- w_2_new[i] - w_delta_new[i]
            w_3_new[i] <- w_3_new[i] - w_delta_new[i]

            # Calculate difference of peak triplet positions
            w_1_2_new <- c(w_1_2_new, w_1_new[i] - w_2_new[i])
            w_1_3_new <- c(w_1_3_new, w_1_new[i] - w_3_new[i])
            w_2_3_new <- c(w_2_3_new, w_2_new[i] - w_3_new[i])

            # Calculate difference of new intensity values of peak triplets
            y_1_2_new <- c(y_1_2_new, y_1_new[i] - y_2_new[i])
            y_1_3_new <- c(y_1_3_new, y_1_new[i] - y_3_new[i])
            y_2_3_new <- c(y_2_3_new, y_2_new[i] - y_3_new[i])

            # Calculate w for each peak triplet
            w_result <- (w_1_new[i]^2 * y_1_new[i] * y_2_3_new[i] + w_3_new[i]^2 * y_3_new[i] * y_1_2_new[i] + w_2_new[i]^2 * y_2_new[i] * (-y_1_3_new[i])) / (2 * w_1_2_new[i] * y_1_new[i] * y_2_new[i] - 2 * (w_1_3_new[i] * y_1_new[i] + (-w_2_3_new[i]) * y_2_new[i]) * y_3_new[i])
            w_result <- w_result + w_delta_new[i]
            w_new <- c(w_new, w_result)

            # If y values are getting 0 after height adjustment, then w_new[i]=NaN
            if (is.nan(w_new[i])) {
                w_new[i] <- 0
            }

            # Calculate lambda for each peak triplet
            lambda_result <- -((sqrt(abs(((-w_2_new[i]^4 * y_2_new[i]^2 * y_1_3_new[i]^2 - w_1_new[i]^4 * y_1_new[i]^2 * y_2_3_new[i]^2 - w_3_new[i]^4 * y_1_2_new[i]^2 * y_3_new[i]^2 + 4 * w_2_new[i] * w_3_new[i]^3 * y_2_new[i] * ((-y_1_new[i]) + y_2_new[i]) * y_3_new[i]^2 + 4 * w_2_new[i]^3 * w_3_new[i] * y_2_new[i]^2 * y_3_new[i] * ((-y_1_new[i]) + y_3_new[i]) + 4 * w_1_new[i]^3 * y_1_new[i]^2 * y_2_3_new[i] * (w_2_new[i] * y_2_new[i] - w_3_new[i] * y_3_new[i]) + 4 * w_1_new[i] * y_1_new[i] * (w_2_new[i]^3 * y_2_new[i]^2 * y_1_3_new[i] - w_2_new[i] * w_3_new[i]^2 * y_2_new[i] * (y_1_new[i] + y_2_new[i] - 2 * y_3_new[i]) * y_3_new[i] + w_3_new[i]^3 * y_1_2_new[i] * y_3_new[i]^2 - w_2_new[i]^2 * w_3_new[i] * y_2_new[i] * y_3_new[i] * (y_1_new[i] - 2 * y_2_new[i] + y_3_new[i])) + 2 * w_2_new[i]^2 * w_3_new[i]^2 * y_2_new[i] * y_3_new[i] * (y_1_new[i]^2 - 3 * y_2_new[i] * y_3_new[i] + y_1_new[i] * (y_2_new[i] + y_3_new[i])) + 2 * w_1_new[i]^2 * y_1_new[i] * (-2 * w_2_new[i] * w_3_new[i] * y_2_new[i] * y_3_new[i] * (-2 * y_1_new[i] + y_2_new[i] + y_3_new[i]) + w_3_new[i]^2 * y_3_new[i] * (y_1_new[i] * (y_2_new[i] - 3 * y_3_new[i]) + y_2_new[i] * (y_2_new[i] + y_3_new[i])) + w_2_new[i]^2 * y_2_new[i] * (y_1_new[i] * (-3 * y_2_new[i] + y_3_new[i]) + y_3_new[i] * (y_2_new[i] + y_3_new[i]))))))))) / (2 * sqrt((w_1_new[i] * y_1_new[i] * y_2_3_new[i] + w_3_new[i] * y_1_2_new[i] * y_3_new[i] + w_2_new[i] * y_2_new[i] * ((-y_1_new[i]) + y_3_new[i]))^2))

            # If y and w are 0, then 0/0=NaN
            if (is.nan(lambda_result)) {
                lambda_result <- 0
            }
            lambda_new <- c(lambda_new, lambda_result)

            # Calculate scaling factor A for each peak triplet
            A_result <- (-4 * w_1_2_new[i] * w_1_3_new[i] * w_2_3_new[i] * y_1_new[i] * y_2_new[i] * y_3_new[i] * (w_1_new[i] * y_1_new[i] * y_2_3_new[i] + w_3_new[i] * y_3_new[i] * y_1_2_new[i] + w_2_new[i] * y_2_new[i] * (-y_1_3_new[i])) * lambda_new[i]) / (w_1_2_new[i]^4 * y_1_new[i]^2 * y_2_new[i]^2 - 2 * w_1_2_new[i]^2 * y_1_new[i] * y_2_new[i] * (w_1_3_new[i]^2 * y_1_new[i] + w_2_3_new[i]^2 * y_2_new[i]) * y_3_new[i] + (w_1_3_new[i]^2 * y_1_new[i] - w_2_3_new[i]^2 * y_2_new[i])^2 * y_3_new[i]^2)

            # If y and w are 0, then 0/0=NaN
            if (is.nan(A_result)) {
                A_result <- 0
            }
            A_new <- c(A_new, A_result)

            # Calculate new lorentz curves
            # If y values are zero, then lorentz curves should also be zero
            if ((w_new[i] == 0) || (lambda_new[i] == 0) || (A_new[i] == 0)) {
                lorentz_curves_initial[i, ] <- 0
            } else {
                lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum_x - w_new[i])^2)))
            }
        }

        # Calculate sum of lorentz curves
        spectrum_approx <- matrix(nrow = 1, ncol = length(spectrum_x))
        for (i in seq_along(spectrum_x)) {
            spectrum_approx[1, i] <- sum(lorentz_curves_initial[seq_along(filtered_peaks), i])
        }

        # Standardization of spectra so that total area equals 1
        spectrum_y_normed <- c()
        spectrum_approx_normed <- c()

        # Standardize the spectra
        spectrum_y_normed <- spectrum_y / sum(spectrum_y)
        spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

        # Calculate the difference between normed original spectrum and normed approximated spectrum
        difference_normed <- c()
        for (i in seq_along(spectrum_x)) {
            difference_normed[i] <- (spectrum_y_normed[i] - spectrum_approx_normed[i])^2
        }
        mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
        msgf("Normed MSE after iteration %d: %.22f", b, mse_normed)
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(lambda_new))
    for (i in seq_along(lambda_new)) {
        integrals[1, i] <- A_new[i] * (atan((-w_new[i] + (spec$n / spec$sfx)) / lambda_new[i]) - atan((-w_new[i]) / lambda_new[i]))
    }

    spec$lcr <- list(
        A_new = A_new,
        lambda_new = lambda_new,
        w_new = w_new,
        spectrum_approx = spectrum_approx,
        spectrum_y_normed = spectrum_y_normed,
        spectrum_approx_normed = spectrum_approx_normed,
        difference_normed = difference_normed,
        mse_normed = mse_normed,
        integrals = integrals
    )
    spec
}

#' @title create backwards compatible return list
#' @param spec Deconvoluted spectrum as returned by [refine_lorentz_curves_v12()].
#' @param n Number of deconvoluted spectrum.
#' @param nam Name of current spectrum.
#' @param debug Add debug info to the return list
add_return_list_v12 <- function(spec = glc_urine1_yy_ni3_dbg()$rv$urine_1,
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
        integrals = spec$lcr$integrals[1, ],
        signal_free_region = c(spec$sfr$left_sdp, spec$sfr$right_sdp),
        range_water_signal_ppm = spec$wsr$hwidth_ppm,
        A = spec$lcr$A_new,
        lambda = spec$lcr$lambda_new,
        x_0 = spec$lcr$w_new
    )
    spec
}

# Helpers for find_peaks_v12 #####

calc_second_derivative_v12 <- function(y) {
    n <- length(y)
    x <- c(NA, y[-n]) # x[i] == y[i-1]
    z <- c(y[-1], NA) # z[i] == y[i+1]
    d <- x + z - 2 * y
    d
}

get_right_border_v12 <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_right_border)) {
            return(r)
        }
        r <- r + 1
    }
    return(NA)
}

get_left_border_v12 <- function(j, d) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_left_border)) {
            return(l)
        }
        l <- l - 1
    }
    return(NA)
}

get_peak_score_v12 <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        0
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}

# Helpers for add_sfrs / add_wsrs #####

#' @title Confirm signal free region (SFR)
#' @description Repeatedly asks the user to refine and/or confirm the borders of the SFR.
#' @param spec A list containing the spectrum details including 'n', 'ppm_max', 'ppm_nstep', and 'sfx'.
#' @param sfr A vector of length 2 containing the initial left and right borders of the SFR in ppm.
#' @return A list containing the final left and right borders of the SFR in both ppm and dp units.
#' @noRd
confirm_sfr <- function(spec, sfr = c(11.44494, -1.8828), ask = TRUE) {
    plot_sfr(spec, sfr[1], sfr[2])
    sfr_ok <- get_yn_input("Signal free region correctly selected?")
    while (!sfr_ok && ask) {
        sfr[1] <- get_num_input("Choose another left border: [e.g. 12]", min = spec$ppm_min, max = spec$ppm_max)
        sfr[2] <- get_num_input("Choose another right border: [e.g. -2]", min = spec$ppm_min, max = spec$ppm_max)
        plot_sfr(spec, sfr[1], sfr[2])
        sfr_ok <- get_yn_input("Signal free region correctly selected?")
    }
    sfr
}

#' @title Confirm water signal half width (WSHW)
#' @description Repeatly asks the user to refine and/or confirm the half width of the water artefact.
#' @param spec A list representing the spectrum.
#' @param wshw The initial half width in ppm.
#' @return The confirmed half width of the water artefact in ppm units.
confirm_wshw <- function(spec, wshw, ask = TRUE) {
    plot_ws(spec, wshw)
    ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    while (!ws_ok) {
        wshw <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
        plot_ws(spec, wshw)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    }
    wshw
}

#' @title Add signal free region (SFR) to the spectrum
#' @description Adds info about the signal free region (SFR) to the spectrum.
#' @param spec A list representing the spectrum.
#' @param sfr A vector of length 2 containing the left and right borders of the SFR in ppm.
#' @return The spectrum with added SFR info.
#' @noRd
add_sfr <- function(spec, sfr) {
    spec$sfr <- within(list(), {
        left_ppm <- sfr[1]
        right_ppm <- sfr[2]
        left_dp <- (spec$n + 1) - (spec$ppm_max - left_ppm) / spec$ppm_nstep
        left_sdp <- left_dp / spec$sfx # nolint: object_usage_linter
        right_dp <- (spec$n + 1) - (spec$ppm_max - right_ppm) / spec$ppm_nstep
        right_sdp <- right_dp / spec$sfx # nolint: object_usage_linter
    })
    spec
}

#' @title Add water signal region
#' @description Adds info about water signal region (WSR) to the spectrum.
#' @param spec A list representing the spectrum.
#' @param wshw The half width of the water artefact in ppm.
#' @return The spectrum with added WSR info.
#' @noRd
add_wsr <- function(spec, wshw) {
    spec$wsr <- within(list(), {
        hwidth_ppm <- wshw
        hwidth_dp <- hwidth_ppm / spec$ppm_nstep # half width in dp
        center_dp <- spec$n / 2 # center line in dp
        right_dp <- center_dp + hwidth_dp # right border in dp
        left_dp <- center_dp - hwidth_dp # left border in dp
        center_ppm <- spec$ppm[center_dp] # center in ppm # nolint: object_usage_linter.
        right_ppm <- spec$ppm[right_dp] # right border in ppm # nolint: object_usage_linter.
        left_ppm <- spec$ppm[left_dp] # left border in ppm # nolint: object_usage_linter.
    })
    spec
}


# Deprecated #####

init_lorentz_curves_v12 <- function(spec) {
    msg("Initializing Lorentz curves")
    w_1 <- c()
    w_2 <- c()
    w_3 <- c()
    y_1 <- c()
    y_2 <- c()
    y_3 <- c()
    w_1_2 <- c()
    w_1_3 <- c()
    w_2_3 <- c()
    y_1_2 <- c()
    y_1_3 <- c()
    y_2_3 <- c()
    w_delta <- c()
    w <- c()
    lambda <- c()
    A <- c()

    spectrum_x <- spec$sdp
    spectrum_y <- spec$y_smooth
    filtered_peaks <- as.integer(spec$peak$center[spec$peak$high] - 1)
    filtered_left_position <- spec$peak$right[spec$peak$high] - 1
    filtered_right_position <- spec$peak$left[spec$peak$high] - 1

    # Calculate parameters w, lambda and A for the initial lorentz curves
    for (i in seq_along(filtered_peaks)) {
        # Calculate position of peak triplets
        w_1 <- c(w_1, spectrum_x[filtered_left_position[i] + 1])
        w_2 <- c(w_2, spectrum_x[filtered_peaks[i] + 1])
        w_3 <- c(w_3, spectrum_x[filtered_right_position[i] + 1])

        # Calculate intensity of peak triplets
        y_1 <- c(y_1, spectrum_y[filtered_left_position[i] + 1])
        y_2 <- c(y_2, spectrum_y[filtered_peaks[i] + 1])
        y_3 <- c(y_3, spectrum_y[filtered_right_position[i] + 1])

        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if (is.na(((y_1[i] < y_2[i]) & (y_2[i] < y_3[i])))) {
            1
        }
        if ((y_1[i] < y_2[i]) && (y_2[i] < y_3[i])) {
            w_3[i] <- 2 * w_2[i] - w_1[i]
            y_3[i] <- y_1[i]
        }
        # For descending shoulders
        if ((y_1[i] > y_2[i]) && (y_2[i] > y_3[i])) {
            w_1[i] <- 2 * w_2[i] - w_3[i]
            y_1[i] <- y_3[i]
        }

        # Move triplet to zero position
        w_delta[i] <- w_1[i]
        w_1[i] <- w_1[i] - w_delta[i]
        w_2[i] <- w_2[i] - w_delta[i]
        w_3[i] <- w_3[i] - w_delta[i]

        # Calculate difference of position of peak triplets
        w_1_2 <- c(w_1_2, w_1[i] - w_2[i])
        w_1_3 <- c(w_1_3, w_1[i] - w_3[i])
        w_2_3 <- c(w_2_3, w_2[i] - w_3[i])

        # Calculate difference of intensity values of peak triplets
        y_1_2 <- c(y_1_2, y_1[i] - y_2[i])
        y_1_3 <- c(y_1_3, y_1[i] - y_3[i])
        y_2_3 <- c(y_2_3, y_2[i] - y_3[i])

        # Calculate w for each peak triplet
        w_result <- (w_1[i]^2 * y_1[i] * y_2_3[i] + w_3[i]^2 * y_3[i] * y_1_2[i] + w_2[i]^2 * y_2[i] * (-y_1_3[i])) / (2 * w_1_2[i] * y_1[i] * y_2[i] - 2 * (w_1_3[i] * y_1[i] + (-w_2_3[i]) * y_2[i]) * y_3[i])
        w_result <- w_result + w_delta[i]
        w <- c(w, w_result)
        # Wenn y Werte nach der H?henanpassung 0 werden, so ist w_new[i] NaN
        if (is.nan(w[i])) {
            w[i] <- 0
        }

        # Calculate lambda for each peak triplet
        lambda_result <- -((sqrt(abs((-w_2[i]^4 * y_2[i]^2 * y_1_3[i]^2 - w_1[i]^4 * y_1[i]^2 * y_2_3[i]^2 - w_3[i]^4 * y_1_2[i]^2 * y_3[i]^2 + 4 * w_2[i] * w_3[i]^3 * y_2[i] * ((-y_1[i]) + y_2[i]) * y_3[i]^2 + 4 * w_2[i]^3 * w_3[i] * y_2[i]^2 * y_3[i] * ((-y_1[i]) + y_3[i]) + 4 * w_1[i]^3 * y_1[i]^2 * y_2_3[i] * (w_2[i] * y_2[i] - w_3[i] * y_3[i]) + 4 * w_1[i] * y_1[i] * (w_2[i]^3 * y_2[i]^2 * y_1_3[i] - w_2[i] * w_3[i]^2 * y_2[i] * (y_1[i] + y_2[i] - 2 * y_3[i]) * y_3[i] + w_3[i]^3 * y_1_2[i] * y_3[i]^2 - w_2[i]^2 * w_3[i] * y_2[i] * y_3[i] * (y_1[i] - 2 * y_2[i] + y_3[i])) + 2 * w_2[i]^2 * w_3[i]^2 * y_2[i] * y_3[i] * (y_1[i]^2 - 3 * y_2[i] * y_3[i] + y_1[i] * (y_2[i] + y_3[i])) + 2 * w_1[i]^2 * y_1[i] * (-2 * w_2[i] * w_3[i] * y_2[i] * y_3[i] * (-2 * y_1[i] + y_2[i] + y_3[i]) + w_3[i]^2 * y_3[i] * (y_1[i] * (y_2[i] - 3 * y_3[i]) + y_2[i] * (y_2[i] + y_3[i])) + w_2[i]^2 * y_2[i] * (y_1[i] * (-3 * y_2[i] + y_3[i]) + y_3[i] * (y_2[i] + y_3[i])))))))) / (2 * sqrt((w_1[i] * y_1[i] * y_2_3[i] + w_3[i] * y_1_2[i] * y_3[i] + w_2[i] * y_2[i] * ((-y_1[i]) + y_3[i]))^2))
        # If y and w are 0, then 0/0=NaN
        if (is.nan(lambda_result)) {
            lambda_result <- 0
        }
        lambda <- c(lambda, lambda_result)

        # Calculate scaling factor A for each peak triplet
        A_result <- (-4 * w_1_2[i] * w_1_3[i] * w_2_3[i] * y_1[i] * y_2[i] * y_3[i] * (w_1[i] * y_1[i] * y_2_3[i] + w_3[i] * y_3[i] * y_1_2[i] + w_2[i] * y_2[i] * (-y_1_3[i])) * lambda[i]) / (w_1_2[i]^4 * y_1[i]^2 * y_2[i]^2 - 2 * w_1_2[i]^2 * y_1[i] * y_2[i] * (w_1_3[i]^2 * y_1[i] + w_2_3[i]^2 * y_2[i]) * y_3[i] + (w_1_3[i]^2 * y_1[i] - w_2_3[i]^2 * y_2[i])^2 * y_3[i]^2)
        # If y and w are 0, then 0/0=NaN
        if (is.nan(A_result)) {
            A_result <- 0
        }
        A <- c(A, A_result)
    }
    spec$lc$A <- A
    spec$lc$lambda <- lambda
    spec$lc$w <- w
    spec$lc$w_delta <- w_delta
    spec
}

init_lorentz_curves_v10 <- function(spectrum_x, spectrum_y, filtered_peaks, filtered_left_position, filtered_right_position, save_scores) {

    w_1 <- c()
    w_2 <- c()
    w_3 <- c()
    y_1 <- c()
    y_2 <- c()
    y_3 <- c()
    w_1_2 <- c()
    w_1_3 <- c()
    w_2_3 <- c()
    y_1_2 <- c()
    y_1_3 <- c()
    y_2_3 <- c()
    w_delta <- c()
    w <- c()
    lambda <- c()
    A <- c()

    # Calculate parameters w, lambda and A for the initial lorentz curves
    for (i in seq_along(filtered_peaks)) {
        # Calculate position of peak triplets
        w_1 <- c(w_1, spectrum_x[filtered_left_position[i] + 1])
        w_2 <- c(w_2, spectrum_x[filtered_peaks[i] + 1])
        w_3 <- c(w_3, spectrum_x[filtered_right_position[i] + 1])

        # Calculate intensity of peak triplets
        y_1 <- c(y_1, spectrum_y[filtered_left_position[i] + 1])
        y_2 <- c(y_2, spectrum_y[filtered_peaks[i] + 1])
        y_3 <- c(y_3, spectrum_y[filtered_right_position[i] + 1])

        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if ((y_1[i] < y_2[i]) && (y_2[i] < y_3[i])) {
            w_3[i] <- 2 * w_2[i] - w_1[i]
            y_3[i] <- y_1[i]
        }
        # For descending shoulders
        if ((y_1[i] > y_2[i]) && (y_2[i] > y_3[i])) {
            w_1[i] <- 2 * w_2[i] - w_3[i]
            y_1[i] <- y_3[i]
        }

        # Move triplet to zero position
        w_delta[i] <- w_1[i]
        w_1[i] <- w_1[i] - w_delta[i]
        w_2[i] <- w_2[i] - w_delta[i]
        w_3[i] <- w_3[i] - w_delta[i]

        # Calculate difference of position of peak triplets
        w_1_2 <- c(w_1_2, w_1[i] - w_2[i])
        w_1_3 <- c(w_1_3, w_1[i] - w_3[i])
        w_2_3 <- c(w_2_3, w_2[i] - w_3[i])

        # Calculate difference of intensity values of peak triplets
        y_1_2 <- c(y_1_2, y_1[i] - y_2[i])
        y_1_3 <- c(y_1_3, y_1[i] - y_3[i])
        y_2_3 <- c(y_2_3, y_2[i] - y_3[i])

        # Calculate w for each peak triplet
        w_result <- (w_1[i]^2 * y_1[i] * y_2_3[i] + w_3[i]^2 * y_3[i] * y_1_2[i] + w_2[i]^2 * y_2[i] * (-y_1_3[i])) / (2 * w_1_2[i] * y_1[i] * y_2[i] - 2 * (w_1_3[i] * y_1[i] + (-w_2_3[i]) * y_2[i]) * y_3[i])
        w_result <- w_result + w_delta[i]
        w <- c(w, w_result)
        # Wenn y Werte nach der H?henanpassung 0 werden, so ist w_new[i] NaN
        if (is.nan(w[i])) {
            w[i] <- 0
        }

        # Calculate lambda for each peak triplet
        lambda_result <- -((sqrt(abs((-w_2[i]^4 * y_2[i]^2 * y_1_3[i]^2 - w_1[i]^4 * y_1[i]^2 * y_2_3[i]^2 - w_3[i]^4 * y_1_2[i]^2 * y_3[i]^2 + 4 * w_2[i] * w_3[i]^3 * y_2[i] * ((-y_1[i]) + y_2[i]) * y_3[i]^2 + 4 * w_2[i]^3 * w_3[i] * y_2[i]^2 * y_3[i] * ((-y_1[i]) + y_3[i]) + 4 * w_1[i]^3 * y_1[i]^2 * y_2_3[i] * (w_2[i] * y_2[i] - w_3[i] * y_3[i]) + 4 * w_1[i] * y_1[i] * (w_2[i]^3 * y_2[i]^2 * y_1_3[i] - w_2[i] * w_3[i]^2 * y_2[i] * (y_1[i] + y_2[i] - 2 * y_3[i]) * y_3[i] + w_3[i]^3 * y_1_2[i] * y_3[i]^2 - w_2[i]^2 * w_3[i] * y_2[i] * y_3[i] * (y_1[i] - 2 * y_2[i] + y_3[i])) + 2 * w_2[i]^2 * w_3[i]^2 * y_2[i] * y_3[i] * (y_1[i]^2 - 3 * y_2[i] * y_3[i] + y_1[i] * (y_2[i] + y_3[i])) + 2 * w_1[i]^2 * y_1[i] * (-2 * w_2[i] * w_3[i] * y_2[i] * y_3[i] * (-2 * y_1[i] + y_2[i] + y_3[i]) + w_3[i]^2 * y_3[i] * (y_1[i] * (y_2[i] - 3 * y_3[i]) + y_2[i] * (y_2[i] + y_3[i])) + w_2[i]^2 * y_2[i] * (y_1[i] * (-3 * y_2[i] + y_3[i]) + y_3[i] * (y_2[i] + y_3[i])))))))) / (2 * sqrt((w_1[i] * y_1[i] * y_2_3[i] + w_3[i] * y_1_2[i] * y_3[i] + w_2[i] * y_2[i] * ((-y_1[i]) + y_3[i]))^2))
        # If y and w are 0, then 0/0=NaN
        if (is.nan(lambda_result)) {
            lambda_result <- 0
        }
        lambda <- c(lambda, lambda_result)

        # Calculate scaling factor A for each peak triplet
        A_result <- (-4 * w_1_2[i] * w_1_3[i] * w_2_3[i] * y_1[i] * y_2[i] * y_3[i] * (w_1[i] * y_1[i] * y_2_3[i] + w_3[i] * y_3[i] * y_1_2[i] + w_2[i] * y_2[i] * (-y_1_3[i])) * lambda[i]) / (w_1_2[i]^4 * y_1[i]^2 * y_2[i]^2 - 2 * w_1_2[i]^2 * y_1[i] * y_2[i] * (w_1_3[i]^2 * y_1[i] + w_2_3[i]^2 * y_2[i]) * y_3[i] + (w_1_3[i]^2 * y_1[i] - w_2_3[i]^2 * y_2[i])^2 * y_3[i]^2)
        # If y and w are 0, then 0/0=NaN
        if (is.nan(A_result)) {
            A_result <- 0
        }
        A <- c(A, A_result)
    }
    list(A = A, lambda = lambda, w = w, w_delta = w_delta)
}
