deconvolute_spectrum <- function(filepath,
                                 name,
                                 file_format,
                                 same_parameter,
                                 processing_value,
                                 number_iterations,
                                 range_water_signal_ppm,
                                 signal_free_region,
                                 smoothing_param,
                                 delta,
                                 scale_factor,
                                 current_filenumber,
                                 number_of_files) {
    msgf("Start deconvolution of %s:", name)

    x <- deconvolution(
        filepath,
        name,
        file_format,
        same_parameter,
        processing_value,
        number_iterations,
        range_water_signal_ppm,
        signal_free_region,
        smoothing_param,
        delta,
        scale_factor,
        current_filenumber
    )
    y <- list(
        "number_of_files" = number_of_files, # [1] add entry
        "filename" = x$filename,
        "x_values" = x$spectrum_x, # [3] rename
        "x_values_ppm" = x$spectrum_x_ppm, # [4] rename
        "y_values" = x$spectrum_y, # [5] rename
        "spectrum_superposition" = x$spectrum_approx, # [6] rename
        "mse_normed" = x$mse_normed,
        "index_peak_triplets_middle" = x$index_peak_triplets_middle,
        "index_peak_triplets_left" = x$index_peak_triplets_left,
        "index_peak_triplets_right" = x$index_peak_triplets_right,
        "peak_triplets_middle" = x$peak_triplets_middle,
        "peak_triplets_left" = x$peak_triplets_left,
        "peak_triplets_right" = x$peak_triplets_right,
        "integrals" = x$integrals,
        "signal_free_region" = x$signal_free_region %||% signal_free_region,
        "range_water_signal_ppm" = x$range_water_signal_ppm %||% range_water_signal_ppm,
        "A" = x$A,
        "lambda" = x$lambda,
        "x_0" = x$w # [19] rename
    )
    return(y)
}


#' @title Deconvolute one single spectrum
#' @description Deconvolute one single spectrum
#' @param path Path to file or folder containing the spectra files.
#' @param type Format of the spectra files. Either `"bruker"` or `"jcampdx"`.
#' @param procno Processing value for the file. E.g. `"10"`. Called `procno` in the Bruker TopSpin Manual.
#' @param expno Spectroscopy value for the file. E.g. `"10"`. Called `expno` in the Bruker TopSpin Manual.
#' @param nfit Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param wshw Half width of the water artefact in ppm.
#' @param sfr Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smopts Row vector with two entries consisting of the number of smoothing repeats for the whole spectrum and the number of data points (uneven) for the mean calculation.
#' @param delta Threshold value to distinguish between signal and noise.
#' @param sf Row vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask Whether the function should ask the user to confirm the signal free region and the water signal. Must be TRUE if either sfr or wshw are not given.
#' @param filno Current file number. Only used for progress prints.
#' @param nfils Total number of files. Only used for progress prints.
#' @param bwc Use the old, slightly incorrect method for calculating the signal free region and water signal to maintain backwards compatibility with MetaboDecon1D results? For details see `Check: ...` issues in `TODOS.md`.
#' @return A list containing the deconvoluted spectrum data.
#' @examples \dontrun{
#' xds_path <- download_example_datasets()
#' path <- file.path(xds_path, "jcampdx/urine/urine_1")
#' type <- "bruker"
#' deconvolute_spectrum_v2(path, type)
#' }
#' @noRd
deconvolute_spectrum_v2 <- function(path = file.path(download_example_datasets(), "bruker/urine/urine_1"),
                                    type = "bruker",
                                    expno = 10,
                                    procno = 10,
                                    nfit = 10,
                                    wshw = 0.1527692,
                                    sfr = c(11.44494, -1.8828),
                                    smopts = c(2, 5),
                                    delta = 6.4,
                                    sf = c(1e3, 1e6),
                                    filno = 1,
                                    nfils = 1,
                                    ask = TRUE,
                                    bwc = TRUE) {
    # Parse arguments
    type <- match.arg(type, c("bruker", "jcampdx"))

    # Load and deconvolute spectrum
    spec <- read_spectrum(path, type, sf, expno, procno)
    spec <- determine_signal_free_region(spec, sfr, ask)
    spec <- determine_water_signal(spec, hwidth_ppm = wshw, bwc, ask)
    spec <- remove_water_signal(spec, bwc)
    spec <- remove_negative_signals(spec)
    spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2], bwc)
    spec <- select_inflection_points_v2(spec, bwc)

    if (!is.null(.GlobalEnv$debugenv)) {
        vcomp(n1 <- debugenv$spectrum_length, n <- n2 <- spec$n)
        vcomp(v1 <- debugenv$spectrum_x, v2 <- spec$sdp)
        vcomp(v1 <- debugenv$spectrum_x_ppm, v2 <- spec$ppm)
        vcomp(v1 <- debugenv$spectrum_y_raw, v2 <- spec$Y$raw)
        vcomp(v1 <- debugenv$spectrum_y_scaled, v2 <- spec$Y$scaled)
        vcomp(v1 <- debugenv$signal_free_region_left, v2 <- spec$sfr$left_sdp)
        vcomp(v1 <- debugenv$signal_free_region_right, v2 <- spec$sfr$right_sdp)
        vcomp(v1 <- debugenv$water_signal_left, v2 <- spec$ws$left_dp)
        vcomp(v1 <- debugenv$water_signal_right, v2 <- spec$ws$right_dp)
        vcomp(v1 <- debugenv$spectrum_y_no_ws, v2 <- spec$Y$nows)
        vcomp(v1 <- debugenv$spectrum_y_no_ws_no_neg_smoothed, v2 <- spec$Y$smooth)
        vcomp(v1 <- debugenv$second_derivative[1, ], v2 <- spec$sdp[2:(n-1)])
        vcomp(v1 <- debugenv$second_derivative[2, ], v2 <- spec$d)
        vcomp(v1 <- debugenv$peaks_index, v2 <- spec$ip)
        vcomp(v1 <- debugenv$peaks_x, v2 <- spec$ip)
    }

    # Peak selection procedure
    left_position <- find_left_positions_v0(ip) # 1)
    find_left_v1 <-
        right_position <- find_right_positions(ip)
    # 1) not quite sure if these are supposed to be left positions of inflection points or left positions of peaks


    # Find left positions of inflection points
    left_position <- matrix(nrow = 1, ncol = length(ip$sdp))
    for (i in 1:length(ip$sdp)) {
        # Save next left position of current local minima
        next_left <- ip$idx[i] + 1
        while ((ip$idx[i] < next_left) & (next_left < length(d2))) {
            if (d2[2, next_left - 1] < d2[2, next_left]) {
                if (((d2[2, next_left - 1] < d2[2, next_left]) & (d2[2, next_left + 1] <= d2[2, next_left])) | ((d2[2, next_left] < 0) & (d2[2, next_left + 1] >= 0))) {
                    left_position[i] <- next_left
                    break
                } else {
                    next_left <- next_left + 1
                }
            } else {
                next_left <- next_left + 1
            }
        }
    }

    # Find all right positions of all local minima of second derivative
    right_position <- matrix(nrow = 1, ncol = length(ip$sdp))
    for (i in 1:length(ip$sdp)) {
        # Save next right position of current local minima
        next_right <- ip$idx[i] - 1
        while ((next_right < ip$idx[i]) & (next_right >= 2)) {
            if (d2[2, next_right + 1] < d2[2, next_right]) {
                if (((d2[2, next_right + 1] < d2[2, next_right]) & (d2[2, next_right - 1] <= d2[2, next_right])) | ((d2[2, next_right] < 0) & (d2[2, next_right - 1] >= 0))) {
                    right_position[i] <- next_right
                    break
                } else {
                    next_right <- next_right - 1
                }
            } else {
                next_right <- next_right - 1
            }
        }
    }


    # Check borders of peak triplets
    # If NA values are available, remove corresponding peak triplet
    for (i in length(left_position):1) {
        if (is.na(left_position[i]) | (is.na(right_position[i]))) {
            ip$sdp <- ip$sdp[-i]
            ip$idx <- ip$idx[-i]
            left_position <- left_position[-i]
            right_position <- right_position[-i]
        }
    }

    # Calculate peak triplet score to distinguish between signal and noise
    scores <- matrix(nrow = 1, ncol = length(ip$sdp))
    scores_left <- matrix(nrow = 1, ncol = length(ip$sdp))
    scores_right <- matrix(nrow = 1, ncol = length(ip$sdp))
    for (i in 1:length(ip$sdp)) {
        # Calculate left score
        left_score <- 0
        for (j in ip$idx[i]:left_position[i]) {
            left_score <- sum(left_score, abs(d2[2, j]))
        }
        scores_left[i] <- left_score
        # Calculate right score
        right_score <- 0
        for (k in right_position[i]:ip$idx[i]) {
            right_score <- sum(right_score, abs(d2[2, k]))
        }
        scores_right[i] <- right_score
        # Save minimum score
        scores[i] <- min(left_score, right_score)
    }

    # Calculate mean of the score and standard deviation of the score of the signal free region R
    index_left <- which(spec$sdp[ip$idx + 1] >= sfrl_sdp)
    index_right <- which(spec$sdp[ip$idx + 1] <= sfrr_sdp)

    mean_score <- mean(c(scores[index_left], scores[index_right]))
    sd_score <- stats::sd(c(scores[index_left], scores[index_right]))

    # Filter peak triplets
    filtered_peaks <- c()
    filtered_left_position <- c()
    filtered_right_position <- c()
    save_scores <- c()
    for (i in 1:length(ip$sdp)) {
        if (scores[i] >= mean_score + delta * sd_score) {
            # Save peak position
            filtered_peaks <- c(filtered_peaks, ip$idx[i])
            # Save left position
            filtered_left_position <- c(filtered_left_position, left_position[i])
            # Save right position
            filtered_right_position <- c(filtered_right_position, right_position[i])
            # Save value of scores of filtered peaks
            save_scores <- c(save_scores, scores[i])
        }
    }


    # Parameter approximation method
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
    for (i in 1:length(filtered_peaks)) {
        # Calculate position of peak triplets
        w_1 <- c(w_1, spec$sdp[filtered_left_position[i] + 1])
        w_2 <- c(w_2, spec$sdp[filtered_peaks[i] + 1])
        w_3 <- c(w_3, spec$sdp[filtered_right_position[i] + 1])

        # Calculate intensity of peak triplets
        y_1 <- c(y_1, spec$y[filtered_left_position[i] + 1])
        y_2 <- c(y_2, spec$y[filtered_peaks[i] + 1])
        y_3 <- c(y_3, spec$y[filtered_right_position[i] + 1])

        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if ((y_1[i] < y_2[i]) & (y_2[i] < y_3[i])) {
            w_3[i] <- 2 * w_2[i] - w_1[i]
            y_3[i] <- y_1[i]
        }
        # For descending shoulders
        if ((y_1[i] > y_2[i]) & (y_2[i] > y_3[i])) {
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

    # Calculate all initial lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spec$sdp))
    for (i in 1:length(filtered_peaks)) {
        # If A = 0, then the lorentz curve is a zero line
        if (A[i] == 0) {
            lorentz_curves_initial[i, ] <- 0
        } else {
            lorentz_curves_initial[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (spec$sdp - w[i])^2)))
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

        for (i in 1:length(filtered_peaks)) {
            # Calculate the position of the peak triplets
            w_1_new <- c(w_1_new, spec$sdp[filtered_left_position[i] + 1])
            w_2_new <- c(w_2_new, spec$sdp[filtered_peaks[i] + 1])
            w_3_new <- c(w_3_new, spec$sdp[filtered_right_position[i] + 1])

            # Calculate the sum of all lorentz curves for each data point
            sum_left[i] <- sum(lorentz_curves_initial[1:length(filtered_left_position), filtered_left_position[i] + 1])
            sum_peaks[i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), filtered_peaks[i] + 1])
            sum_right[i] <- sum(lorentz_curves_initial[1:length(filtered_right_position), filtered_right_position[i] + 1])

            # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
            proportion_left[i] <- spec$y[filtered_left_position[i] + 1] / sum_left[i]
            proportion_peaks[i] <- spec$y[filtered_peaks[i] + 1] / sum_peaks[i]
            proportion_right[i] <- spec$y[filtered_right_position[i] + 1] / sum_right[i]

            # Calculate the new heights of the peak triplets
            y_1_new[i] <- lorentz_curves_initial[i, filtered_left_position[i] + 1] * proportion_left[i]
            y_2_new[i] <- lorentz_curves_initial[i, filtered_peaks[i] + 1] * proportion_peaks[i]
            y_3_new[i] <- lorentz_curves_initial[i, filtered_right_position[i] + 1] * proportion_right[i]

            # Calculate mirrored points if necesccary
            # For ascending shoulders
            if ((y_1_new[i] < y_2_new[i]) & (y_2_new[i] < y_3_new[i])) {
                w_3_new[i] <- 2 * w_2_new[i] - w_1_new[i]
                y_3_new[i] <- y_1_new[i]
            }
            # For descending shoulders
            if ((y_1_new[i] > y_2_new[i]) & (y_2_new[i] > y_3_new[i])) {
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
            if ((w_new[i] == 0) | (lambda_new[i] == 0) | (A_new[i] == 0)) {
                lorentz_curves_initial[i, ] <- 0
            } else {
                lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spec$sdp - w_new[i])^2)))
            }
        }

        # Calculate sum of lorentz curves
        spectrum_approx <- matrix(nrow = 1, ncol = length(spec$sdp))
        for (i in 1:length(spec$sdp)) {
            spectrum_approx[1, i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), i])
        }
        # ToSc: use vectorized functions, e.g.
        # spectrum_approx <- colSums(lorentz_curves_initial[1:length(filtered_peaks), ])

        # Standardize the spectra so that total area equals 1
        spectrum_y_normed <- spec$y / sum(spec$y)
        spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

        # Calculate the difference between normed original spectrum and normed approximated spectrum
        difference_normed <- c()
        for (i in 1:length(spec$sdp)) {
            difference_normed[i] <- (spectrum_y_normed[i] - spectrum_approx_normed[i])^2
        }
        mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
        message(paste("\nNormed MSE value of iteration", b, "is: "))
        print(mse_normed)
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(lambda_new))
    for (i in 1:length(lambda_new)) {
        integrals[1, i] <- A_new[i] * (atan((-w_new[i] + (spec$length / sfx)) / lambda_new[i]) - atan((-w_new[i]) / lambda_new[i]))
    }


    # Save index of peak triplets
    index_peak_triplets_middle <- c()
    index_peak_triplets_left <- c()
    index_peak_triplets_right <- c()
    for (i in 1:length(filtered_peaks)) {
        index_peak_triplets_middle[i] <- filtered_peaks[i] + 1
        index_peak_triplets_left[i] <- filtered_left_position[i] + 1
        index_peak_triplets_right[i] <- filtered_right_position[i] + 1
    }

    # Save ppm x position of peak triplets
    peak_triplets_middle <- c()
    peak_triplets_left <- c()
    peak_triplets_right <- c()
    for (i in 1:length(filtered_peaks)) {
        peak_triplets_middle[i] <- spec$x_ppm[index_peak_triplets_middle[i]]
        peak_triplets_left[i] <- spec$x_ppm[index_peak_triplets_left[i]]
        peak_triplets_right[i] <- spec$x_ppm[index_peak_triplets_right[i]]
    }

    # Save values A_new, lambda_new, w_new and noise_threshold to txt document
    noise_threshold <- replicate(length(w_new), 0)
    noise_threshold[1] <- mean_score + delta * sd_score
    spectrum_info <- data.frame(rbind(w_new, lambda_new, A_new, noise_threshold))
    spectrum_output <- data.frame(spectrum_approx)
    name_info_txt <- paste(name, "parameters.txt")
    name_output_txt <- paste(name, "approximated_spectrum.txt")

    message(paste("\nSaving parameters to txt documents..."))
    utils::write.table(spectrum_info, name_info_txt, sep = ",", col.names = FALSE, append = FALSE)
    utils::write.table(spectrum_output, name_output_txt, sep = ",", col.names = FALSE, append = FALSE)

    return_list <- list(
        "filename" = name,
        "spectrum_x" = spec$sdp,
        "spectrum_x_ppm" = spec$x_ppm,
        "spectrum_y" = spec$y,
        "lorentz_curves" = lorentz_curves_initial,
        "mse_normed" = mse_normed,
        "spectrum_approx" = spectrum_approx,
        "index_peak_triplets_middle" = index_peak_triplets_middle,
        "index_peak_triplets_left" = index_peak_triplets_left,
        "index_peak_triplets_right" = index_peak_triplets_right,
        "peak_triplets_middle" = peak_triplets_middle,
        "peak_triplets_left" = peak_triplets_left,
        "peak_triplets_right" = peak_triplets_right,
        "integrals" = integrals,
        "sfr" = c(sfrl_sdp, sfrr_sdp),
        "ws$hwidth_ppm" = ws$hwidth_ppm,
        "A" = A_new,
        "lambda" = lambda_new,
        "w" = w_new
    )
    return(return_list)
}

#' @title Determine Signal Free Region
#' @description This function determines the signal free region (SFR) of a given spectrum. It asks the user to confirm the left and right borders of the SFR, and allows them to adjust these borders if necessary. The function returns a list containing the left and right borders in both ppm and data points (dp), as well as the scaled data points (sdp).
#' @param spec A list representing the spectrum, which should include the minimum and maximum ppm (`$ppm_min` and `$ppm_max`), and the scaling factor (`$sfx`).
#' @param sfr Initial values for the left and right borders of the SFR in ppm. If not provided, the function will ask the user to select the borders.
#' @param ask Logical. If TRUE, the function will ask the user to confirm or adjust the borders of the SFR. Default is TRUE.
#' @param bwc Use the old, slightly incorrect method for conversion from ppm to data points to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: ppm to dp conversion` in TODOS.md
#' @return A list containing the left and right borders of the SFR in ppm (`$left_ppm` and `$right_ppm`), data points (`$left_dp` and `$right_dp`), and scaled data points (`$left_sdp` and `$right_sdp`).
#' @noRd
determine_signal_free_region <- function(spec, sfr = NULL, ask = TRUE, bwc = TRUE) {
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    if (is.null(sfr) && isFALSE(ask)) {
        stop("No signal free region (SFR) provided and `ask` is FALSE. Please provide the SFR or set `ask` to TRUE.")
    }
    if (ask) {
        plot_sfr(spec, left_ppm, right_ppm)
        sfr_ok <- get_yn_input("Signal free region borders correct selected? (Area left and right of the green lines)")
        while (!sfr_ok) {
            left_ppm <- get_num_input("Choose another left border: [e.g. 12]", min = spec$ppm_min, max = spec$ppm_max)
            right_ppm <- get_num_input("Choose another right border: [e.g. -2]", min = spec$ppm_min, max = spec$ppm_max)
            plot_sfr(spec, left_ppm, right_ppm)
            sfr_ok <- get_yn_input("Signal free region borders correct selected? (Area left and right of the green lines)")
        }
    }
    left_dp <- ppm_to_dp(left_ppm, spec, bwc)
    left_sdp <- left_dp / spec$sfx
    right_dp <- ppm_to_dp(right_ppm, spec, bwc)
    right_sdp <- right_dp / spec$sfx
    spec$sfr <- list(left_ppm = left_ppm, right_ppm = right_ppm, left_sdp = left_sdp, right_sdp = right_sdp, left_dp = left_dp, right_dp = right_dp)
    spec
}

#' @title Calculate water signal parameters
#' @description Calculates water signal parameters for a given spectrum.
#' @param spec A list representing the spectrum.
#' @param hwidth_ppm The half width in ppm. Default is wshw.
#' @param bwc Use the old, slightly incorrect methods for calculating water signal values to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: water signal calculation` in `TODOS.md`.
#' @return List of parameters including half width in dp and ppm, center line in dp and ppm and right and left borders in dp and ppm.
determine_water_signal <- function(spec, hwidth_ppm, bwc = TRUE, ask = TRUE) {
    if (ask) {
        plot_ws(spec, hwidth_ppm)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
        while (!ws_ok) {
            hwidth_ppm <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
            plot_ws(spec, hwidth_ppm)
            ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
        }
    }
    hwidth_dp <- if (bwc) hwidth_ppm / spec$ppm_nstep else hwidth_ppm / spec$ppm_step # half width in dp
    center_dp <- if (bwc) spec$n / 2 else (spec$n - 1) / 2 # center line in dp
    right_dp <- center_dp + hwidth_dp # right border in dp
    left_dp <- center_dp - hwidth_dp # left border in dp
    center_ppm <- if (bwc) spec$ppm[center_dp] else dp_to_ppm(center_dp, spec) # center in ppm
    right_ppm <- if (bwc) spec$ppm[right_dp] else dp_to_ppm(right_dp, spec) # right border in ppm
    left_ppm <- if (bwc) spec$ppm[left_dp] else dp_to_ppm(left_dp, spec) # left border in ppm
    spec$ws <- list(center_dp = center_dp, hwidth_dp = hwidth_dp, left_dp = left_dp, right_dp = right_dp, center_ppm = center_ppm, hwidth_ppm = hwidth_ppm, left_ppm = left_ppm, right_ppm = right_ppm)
    spec
}

remove_water_signal <- function(spec, bwc = TRUE) {
    if (is.null(spec$ws)) {
        stop("No water signal parameters found. Please call `determine_water_signal()`.")
    }
    y <- spec$Y$scaled
    if (bwc) {
        left <- spec$ws$left_dp
        right <- spec$ws$right_dp
        y[right:left] <- 0.01 / spec$sfy # (1)
        # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats in the backwards compatible case. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    } else {
        ppm <- spec$ppm
        ws_min <- spec$ws$center_ppm - spec$ws$hwidth_ppm
        ws_max <- spec$ws$center_ppm + spec$ws$hwidth_ppm
        ws_idx <- which(ppm >= ws_min & ppm <= ws_max)
        y[wsidx] <- 0
    }
    spec$Y$nows <- y
    spec
}

remove_negative_signals <- function(spec) {
    if (is.null(spec$Y$nows)) stop("Water signal not removed yet. Please call `remove_water_signal()` first.")
    spec$Y$pos <- abs(spec$Y$nows)
    spec
}

#' @inherit smooth_signals_v1
#' @param bwc Maintain backwards compatibility with MetaboDecon1D results by using the old and slow method for smoothing ([smooth_signals_v1()])? If FALSE, the new and fast method ([smooth_signals_v2()]) is used instead.
#' @noRd
smooth_signals <- function(spec, reps = 2, k = 5, bwc = TRUE) {
    if (bwc) {
        smooth_signals_v1(spec, reps, k)
    } else {
        smooth_signals_v2(spec, reps, k)
    }
}

#' @inherit smooth_signals_v1
#' @details New and fast version for smoothing of signals. Implements the same algorithm as `smooth_signal_v1` using different R functions (e.g. [stats::filter()]), causing a massive speedup but also numeric differences compared to the old version.
#' @noRd
smooth_signals_v2 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")

    Z <- vector("list", length = reps)
    y <- spec$Y$pos
    n <- length(y)

    for (i in 1:reps) {
        filter <- rep(1 / k, k)
        z <- stats::filter(y, filter, sides = 2) # (1)
        q <- (k - 1) / 2 # (2)
        for (j in seq_len(q)) {
            z[j] <- mean(y[1:(q + j)]) # (3)
            z[n - j + 1] <- mean(y[(n - q - j + 1):n]) # (4)
        }
        y <- Z[[i]] <- as.numeric(z)
        # Calling (1) gives NAs at both sides of vector, as there are not enough values for the moving average. The number of NAs at each side is given by (2). Example: if n==100 and k==5, then q==2, so z[1]==NA, z[2]==NA, z[99]==NA and z[100]==NA. To stay backwards compatible, these values must be filled with the mean of the values that are available. To do so, we iterate from 1:q, i.e. j==1 and j==2 and set
        # z[1]   <- mean(y[1:3])    # 3 == 2+1 == q+j            # (3)
        # z[2]   <- mean(y[1:4])    # 4 == 2+2 == q+j            # (3)
        # z[100] <- mean(y[98:100]) # 98 == 100-2-1+1 == n-q-j+1 # (4)
        # z[99]  <- mean(y[97:100]) # 97 == 100-2-2+1 == n-q-j+1 # (4)
        # Note: we could also think of leaving the NAs as they are, which would be more correct I think and even faster, but would break compatibility with the old version completely. So not even `all.equal(v1, v2)` would be TRUE anymore.
    }

    spec$Z <- Z
    spec$Y$smooth <- Z[[reps]]
    spec
}

#' @title Smooth signal intensities using a moving average
#' @description This function smooths signal intensities by applying a [moving average](https://en.wikipedia.org/wiki/Moving_average) filter with a window size of k.
#' @param spec A list representing the spectrum, which should include the scaled signal intensities, after removal of the water artefact and negative values (`spec$Y$pos`).
#' @param reps The number of times to apply the moving average.
#' @param k The number of points within the moving average window. Must be odd, so the smoothed point is in the middle of the window.
#' @return A numeric vector of the smoothed values.
#' @details Old and slow version producing the same results as the implementation within `deconvolution` from `MetaboDecon1D_deconvolution.R`.
#' @noRd
smooth_signals_v1 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")
    Z <- vector("list", length = reps)
    y <- spec$Y$pos
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
    spec$Y$smooth <- Z[[reps]]
    spec
}

#' @title Plot Water Signal
#' @description Draws the water signal as red vertical lines into the given spectrum.
#' @param spec A list representing the spec as returned by [load_jcampdx_spectrum()] or [load_bruker_spectrum()].
#' @param hwidth_ppm The half width of the water signal in ppm.
#' @return NULL. Called for side effect of plotting the water signal.
plot_ws <- function(spec, hwidth_ppm) {
    center_ppm <- (spec$ppm_max + spec$ppm_min) / 2
    plot(
        spec$ppm,
        spec$y,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(center_ppm + 2 * hwidth_ppm, center_ppm - 2 * hwidth_ppm)
    )
    graphics::abline(v = center_ppm + hwidth_ppm, col = "red")
    graphics::abline(v = center_ppm - hwidth_ppm, col = "red")
}

#' @title Plot Signal Free Region
#' @description Draws the signal free region as green vertical lines into the given spectrum.
#' @param spec A list representing the spectrum as returned by [load_jcampdx_spectrum()] or [load_bruker_spectrum()].
#' @param left_ppm The left border of the signal free region in ppm.
#' @param right_ppm The right border of the signal free region in ppm.
#' @return NULL. Called for side effect of plotting the signal free region.
plot_sfr <- function(spec, left_ppm, right_ppm) {
    plot(
        x = spec$ppm,
        y = spec$y,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(spec$ppm_max, spec$ppm_min)
    )
    graphics::abline(v = left_ppm, col = "green")
    graphics::abline(v = right_ppm, col = "green")
}

#' @title Select Inflection Points
#' @description This function calculates the second derivative of a spec and finds all local minima of the derivate, i.e. inflection points. It returns a list containing the second derivative, the x-values of the local minima, and their indices.
#' @param spec A list with elements `sdp` and `ssss` representing the x-values in scaled data points and the y values in "smoothed and scaled signal strength" respectively.
#' @param details this is the code copy pasted from the original `deconvolution` function of `MetaboDecon1D` with only the inputs renamed to `spec$sdp` and `spec$ssss`.
select_inflection_points_v0 <- function(spec) {
    # Calculate second derivative of spec
    second_derivative <- matrix(nrow = 2, ncol = length(spec$sdp) - 2)
    for (i in 2:length(spec$sdp) - 1) {
        second_derivative[1, i - 1] <- spec$sdp[i]
        second_derivative[2, i - 1] <- spec$ssss[i - 1] + spec$ssss[i + 1] - 2 * spec$ssss[i]
    }

    # Find all local minima of second derivative
    peaks_x <- c()
    peaks_index <- c()
    second_derivative_border <- ncol(second_derivative) - 1
    for (i in 2:second_derivative_border) {
        if (second_derivative[2, i] < 0) {
            if ((second_derivative[2, i] <= second_derivative[2, i - 1]) & (second_derivative[2, i] < second_derivative[2, i + 1])) {
                # Add local minima to peak list
                peaks_x <- c(peaks_x, second_derivative[1, i])
                peaks_index <- c(peaks_index, i)
            }
        }
    }
    return(list(
        d2 = second_derivative[2, ],
        sdp = peaks_x,
        idx = peaks_index
    ))
}

#' @inherit select_inflection_points_v0
#' @param details This is the same as `select_inflection_points_v0`, but using different variable names within the code to make it more readable.
select_inflection_points_v1 <- function(spec) {
    y <- spec$ssss
    x <- spec$sdp
    n <- spec$n
    # Calculate second derivative of smoothed and scaled data points
    d <- matrix(nrow = 2, ncol = n - 2)
    for (i in 2:n - 1) { # Hack: this should be 2:(n-1), because we want to iterate from 2 to (n-1), but the actual code iterates from 1 to (n-1), so we have `d[1, 0] <- x[1]` in the first iteration. It's not too bad, because R silently ignores this error, so the code works. But again, it's very hard to reason about the code.
        d[1, i - 1] <- x[i]
        d[2, i - 1] <- y[i - 1] + y[i + 1] - 2 * y[i]
    }
    # Find all inflection points, i.e. local minima of second derivative
    ipx <- c()
    ipi <- c()
    for (i in 2:(ncol(d) - 1)) {
        if (d[2, i] < 0) {
            if ((d[2, i] <= d[2, i - 1]) & (d[2, i] < d[2, i + 1])) {
                # Add local minima to inflection point list
                ipx <- c(ipx, d[1, i])
                ipi <- c(ipi, i)
            }
        }
    }
    return(list(d2 = d[2, ], sdp = ipx, idx = ipi))
}


select_inflection_points_v2 <- function(spec, bwc = TRUE) {
    y <- spec$Y$smooth
    d <- calculate_second_derivative(y, bwc)
    m <- length(d) # Calculating second derivative of y removes points 1 and n as they dont have neighbours that can be used for the calculation
    d1 <- d[1:(m - 2)]
    d2 <- d[2:(m - 1)]
    d3 <- d[3:(m - 0)]
    spec$ip <- as.integer(which(d2 < 0 & d2 <= d1 & d2 < d3) + 1) # TODO: why <= instad of < once?
    spec$d <- d
    return(spec)
}

find_left_positions_v0 <- function(ip) {
    left_position <- matrix(nrow = 1, ncol = length(ip$sdp))
    d2 <- ip$d2
    idx <- ip$idx
    for (i in 1:length(ip$sdp)) {
        # Save next left position of current local minima
        next_left <- idx[i] + 1
        while ((idx[i] < next_left) & (next_left < length(d2))) {
            if (d2[next_left - 1] < d2[next_left]) {
                if (((d2[next_left - 1] < d2[next_left]) & (d2[next_left + 1] <= d2[next_left])) |
                    ((d2[next_left] < 0) & (d2[next_left + 1] >= 0)))
                {
                    left_position[i] <- next_left
                    break
                } else {
                    next_left <- next_left + 1
                }
            } else {
                next_left <- next_left + 1
            }
        }
    }
    left_position
}

find_left_positions_v1 <- function(ip) {
    d2 <- ip$d2
    sdp <- ip$sdp
    idx <- ip$idx
    left <- rep(NA, length(idx))
    for (i in seq_along(idx)) {
        candidates <- idx[i]:(idx[i + 1] - 1)
        for (j in candidates) { # j == <index of current inflection point>
            c0 <- d2[j - 1] # curvature of data point before current inflection point
            c1 <- d2[j + 0] # curvature of current inflection point
            c2 <- d2[j + 1] # curvature of data point after current inflection point
            cond1 <- c0 < c1 & c1 <= c2
            cond2 <- c0 < c1 & c1 < 0 & c2 >= 0
            if (cond1 | cond2) {
                left[i] <- j
                break
            }
        }
    }
    left
}

calculate_second_derivative <- function(y, bwc) {
    n <- length(y)
    if (bwc) {
        y1 <- y[1:(n - 2)]
        y2 <- y[2:(n - 1)]
        y3 <- y[3:(n - 0)]
        d <- y1 + y3 - 2 * y2
    } else {
        # Using diff is almost equivalent to the above, but due to numeric instabilities, it sometimes gives slightly different results (e.g. -5.51600000000001e-05 instead of -5.51599999999998e-05). Since we do `<=` and `<` comparisons further below this is not backwards compatible. Nonetheless, we keep it here, because it would make the code slightly more readable (runtime is almost the same).
        d <- diff(y, differences = 2)
    }
}
