# Private API Functions #####

#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spetra and generates a Lorentz curve for each detected signal within a spectra.
#' @param data_path Either the path to an existing directory containing measured NMR spectra or a dataframe with columns `ppm` (parts per million) and `si` (signal intensity) or a list of such dataframes.
#' @param file_format Format of the spectra files. Either `"bruker"` or `"jcampdx"`. Only relevant if `data_path` is a directory.
#' @param make_rds Store results as a rds file on disk? Should be set to TRUE if many spectra are evaluated to decrease computation time.
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
#' path <- file.path(xds_path, "jcampdx/urine/urine_1")
#' type <- "bruker"
#' deconvolution_v12(path, type)
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
                                        debug = FALSE) {
    spectra <- get_spectra(data_path, file_format, expno, procno, ask, sf)
    spectra <- determine_wsr_and_sfr(spectra, sfr, wshw, ask)
    nfils <- length(spectra)
    for (i in seq_along(spectra)) {
        filno <- i
        spec <- spectra[[i]]
        spec <- remove_water_signal_v12(spec)
        spec <- remove_negative_signals_v12(spec)
        spec <- smooth_signals_v12(spec, reps = smopts[1], k = smopts[2])
        spec <- find_peaks_v12(spec)
        spec <- rm_peaks_with_low_scores_v12(spec, delta)
        spec <- init_lorentz_curves_v12(spec)
        spec <- refine_lorentz_curves_v12(spec, nfit)
        spec <- create_return_list_v12(spec, args)
        # check_spec(spec, compare_against = "deconvolution_urine1_spF_ni10_cf1_nf2.rds")
        spectra[[i]] <- spec
    }
    return(spectra)
}

# Private helpers #####

get_spectra <- function(data_path, file_format, expno, procno, ask, sf) {
    if (!file_format %in% c("bruker", "jcampdx")) {
        stop("Argument `file_format` should be either 'bruker' or 'jcampdx'")
    }
    p <- data_path
    spectra <- if (is.character(p) && length(p) == 1 && dir.exists(p)) {
        read_spectra(path = p, type = file_format, expno, procno, ask, sf)
    } else if (is.data.frame(p) && all(c("ppm", "si") %in% colnames(p))) {
        list(p)
    } else if (is.list(p) && all(sapply(p, function(p) all(c("ppm", "si") %in% colnames(p))))) {
        p
    } else {
        stop("Argument `data_path` should be either:\n- a character vector of length 1, denoting an existing directory containing measured NMR spectra or\n- a dataframe with columns `ppm` (parts per million) and `si` (signal intensity) or\n- a list of dataframes with columns `ppm` and `si`")
    }
}

determine_wsr_and_sfr <- function(spectra, sfr, wshw, ask) {
    # length(spectra): 1, n
    # min(length(sfr), length(wshw)): 1, n
    # ask: F, T
    n <- length(spectra)
    sfr_info <- sprintf("Argument `sfr` should be either\n- a numeric vector of length 2, giving the left and right boundaries of the signal free region in ppm or\n- a list of length %d of such vectors (one for each spectrum)", n)
    wshw_info <- sprintf("Argument `wshw` should be either\n- a numeric value, giving the half width of the water artefact in ppm or\n- a vector of length %d of such values (one for each spectrum)", n)
    sfr_defaults <- if (is_list_of_nums(sfr, n, 2)) sfr else if (is_num(sfr, 2)) rep(list(sfr), n) else stop(sfr_info)
    wshw_defaults <- if (is_list_of_nums(wshw, n, 1)) wshw else if (is_num(wshw, 1)) rep(list(wshw), n) else stop(wshw_info)
    if (isFALSE(ask)) {
        for (i in seq_along(spectra)) {
            spectra[[i]]$sfr <- sfr_defaults[[i]]
            spectra[[i]]$wshw <- wshw_defaults[[i]]
        }
        return(spectra)
    }
    same_parameter <- get_yn_input("Use same parameters (sfr, wshw) for all spectra?")
    if (same_parameter) {
        msgf("[%d] %s\n", seq_along(spectra), names(spectra), appendLF = FALSE)
        main_spectrum <- get_num_input("Number of spectrum for adjusting parameters (e.g. 1):", min = 1, max = length(files), int = TRUE)
        # CONTINUE HERE
    }
    return(spectra)
}

create_return_list_v12 <- function(spec, number_of_files, debug = debug) {
    # args: path, type, expno, procno, nfit, wshw, sfr, smopts, delta, sf, filno, nfils, ask
    # spec: Y, n, sfx, sfy, dp, sdp, ppm, ppm_min, ppm_max, ppm_range, ppm_step, ppm_nstep, sfr, ws, Z, d, peak, lc, lcr
    ret <- list(
        number_of_files = number_of_files,
        filename = basename(args$path),
        x_values = spec$sdp,
        x_values_ppm = spec$x_ppm,
        y_values = spec$Y$smooth,
        spectrum_superposition = spec$lcr$spectrum_approx,
        mse_normed = spec$lcr$mse_normed,
        index_peak_triplets_middle = as.integer(spec$peak$center[spec$peak$high] - 1),
        index_peak_triplets_left = spec$peak$right[spec$peak$high] - 1,
        index_peak_triplets_right = spec$peak$left[spec$peak$high] - 1,
        peak_triplets_middle = spec$peak$peak_triplets_middle,
        peak_triplets_left = spec$peak$peak_triplets_left,
        peak_triplets_right = spec$peak$peak_triplets_right,
        integrals = spec$lcr$integrals,
        sfr = c(spec$sfr$left_sdp, spec$sfr$right_sdp),
        wshw = spec$ws$hwidth_ppm,
        A = spec$lcr$A_new,
        lambda = spec$lcr$lambda_new,
        x_0 = spec$lcr$w_new
    )
    spec$ret <- ret
}

#' @title Determine Signal Free Region
#' @description This function determines the signal free region (SFR) of a given spectrum. It asks the user to confirm the left and right borders of the SFR, and allows them to adjust these borders if necessary. The function returns a list containing the left and right borders in both ppm and data points (dp), as well as the scaled data points (sdp).
#' @param spec A list representing the spectrum, which should include the minimum and maximum ppm (`$ppm_min` and `$ppm_max`), and the scaling factor (`$sfx`).
#' @param sfr Initial values for the left and right borders of the SFR in ppm. If not provided, the function will ask the user to select the borders.
#' @param ask Logical. If TRUE, the function will ask the user to confirm or adjust the borders of the SFR. Default is TRUE.
#' @return A list containing the left and right borders of the SFR in ppm (`$left_ppm` and `$right_ppm`), data points (`$left_dp` and `$right_dp`), and scaled data points (`$left_sdp` and `$right_sdp`).
#' @noRd
determine_signal_free_region_v12 <- function(spec, sfr = NULL, ask = TRUE) {
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
    spec$sfr <- within(list(), {
        left_ppm <- left_ppm
        right_ppm <- right_ppm
        left_dp <- (spec$n + 1) - (spec$ppm_max - left_ppm) / spec$ppm_nstep
        left_sdp <- left_dp / spec$sfx
        right_dp <- (spec$n + 1) - (spec$ppm_max - right_ppm) / spec$ppm_nstep
        right_sdp <- right_dp / spec$sfx
    })
    spec
}

#' @title Calculate water signal parameters
#' @description Calculates water signal parameters for a given spectrum.
#' @param spec A list representing the spectrum.
#' @param hwidth_ppm The half width in ppm. Default is wshw.
#' @return List of parameters including half width in dp and ppm, center line in dp and ppm and right and left borders in dp and ppm.
determine_water_signal_v12 <- function(spec, hwidth_ppm, ask = TRUE) {
    if (ask) {
        plot_ws(spec, hwidth_ppm)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
        while (!ws_ok) {
            hwidth_ppm <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
            plot_ws(spec, hwidth_ppm)
            ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
        }
    }
    spec$ws <- within(list(), {
        hwidth_ppm <- hwidth_ppm
        hwidth_dp <- hwidth_ppm / spec$ppm_nstep # half width in dp
        center_dp <- spec$n / 2 # center line in dp
        right_dp <- center_dp + hwidth_dp # right border in dp
        left_dp <- center_dp - hwidth_dp # left border in dp
        center_ppm <- spec$ppm[center_dp] # center in ppm
        right_ppm <- spec$ppm[right_dp] # right border in ppm
        left_ppm <- spec$ppm[left_dp] # left border in ppm
    })
    spec
}

remove_water_signal_v12 <- function(spec) {
    y <- spec$Y$scaled
    left <- spec$ws$left_dp
    right <- spec$ws$right_dp
    y[right:left] <- 0.01 / spec$sfy # Order, i.e. `right:left` instead of `left:right`, is important here, because `right` and `left` are floats. Example: `right <- 3.3; left <- 1.4` ==> `right:left == c(3.3, 2.3)` and `left:right == c(1.4, 2.4)`.
    spec$Y$nows <- y
    spec
}

remove_negative_signals_v12 <- function(spec) {
    if (is.null(spec$Y$nows)) stop("Water signal not removed yet. Please call `remove_water_signal_v12()` first.")
    spec$Y$pos <- abs(spec$Y$nows)
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
smooth_signals_v12 <- function(spec, reps = 2, k = 5) {
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

#' @inherit find_peaks_v12
#' @param details Successor of `select_peaks_v0`, `find_left_positions_v0` and `find_right_positions_v0`. The new function `find_peaks_v12()` is a combination of these three functions with a corrected naming convention: what was incorrectly referred to as "left" is now correctly called "right" and vice versa.
#' noRd
find_peaks_v12 <- function(spec) {
    d <- spec$d <- calc_second_derivative_v12(y = spec$Y$smooth)
    a <- abs(d)
    m <- length(d)
    dl <- c(NA, d[-m]) # dl[i] == d[i-1]
    dr <- c(d[-1], NA) # dr[i] == d[i+1]
    center <- which(d < 0 & d <= dl & d < dr)
    spec$peak <- data.frame(left = NA, center = center, right = NA, score = NA)
    for (i in seq_along(center)) {
        j <- center[i]
        l <- spec$peak$left[i] <- get_left_border_v12(j, d, m)
        r <- spec$peak$right[i] <- get_right_border_v12(j, d, m)
        s <- spec$peak$score[i] <- get_peak_score_v12(j, l, r, a)
    }
    return(spec)
}

rm_peaks_with_low_scores_v12 <- function(spec, delta = 6.4) {
    score <- spec$peak$score
    l <- which(spec$sdp[spec$peak$center] >= spec$sfr$left_sdp)
    r <- which(spec$sdp[spec$peak$center] <= spec$sfr$right_sdp)
    mu <- mean(score[c(l, r)])
    sigma <- sd(c(score[l], score[r]))
    spec$peak$high <- score >= mu + delta * sigma
    spec$peak$region <- "norm"
    spec$peak$region[l] <- "sfrl"
    spec$peak$region[r] <- "sfrr"
    spec
}

init_lorentz_curves_v12 <- function(spec) {
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
    spectrum_y <- spec$Y$smooth
    filtered_peaks <- as.integer(spec$peak$center[spec$peak$high] - 1)
    filtered_left_position <- spec$peak$right[spec$peak$high] - 1
    filtered_right_position <- spec$peak$left[spec$peak$high] - 1
    save_scores <- spec$peak$score[spec$peak$high]

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

    spec$lc$A <- A
    spec$lc$lambda <- lambda
    spec$lc$w <- w
    spec$lc$w_delta <- w_delta
    spec
}

refine_lorentz_curves_v12 <- function(spec, nfit) {
    spectrum_x <- spec$sdp
    spectrum_y <- spec$Y$smooth
    filtered_peaks <- as.integer(spec$peak$center[spec$peak$high] - 1)
    filtered_left_position <- spec$peak$right[spec$peak$high] - 1
    filtered_right_position <- spec$peak$left[spec$peak$high] - 1
    A <- spec$lc$A
    lambda <- spec$lc$lambda
    w <- spec$lc$w
    nfit <- nfit

    # Calculate all initial lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spectrum_x))
    for (i in 1:length(filtered_peaks)) {
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

        for (i in 1:length(filtered_peaks)) {
            # Calculate the position of the peak triplets
            w_1_new <- c(w_1_new, spectrum_x[filtered_left_position[i] + 1])
            w_2_new <- c(w_2_new, spectrum_x[filtered_peaks[i] + 1])
            w_3_new <- c(w_3_new, spectrum_x[filtered_right_position[i] + 1])

            # Calculate the sum of all lorentz curves for each data point
            sum_left[i] <- sum(lorentz_curves_initial[1:length(filtered_left_position), filtered_left_position[i] + 1])
            sum_peaks[i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), filtered_peaks[i] + 1])
            sum_right[i] <- sum(lorentz_curves_initial[1:length(filtered_right_position), filtered_right_position[i] + 1])

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
                lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum_x - w_new[i])^2)))
            }
        }

        # Calculate sum of lorentz curves
        spectrum_approx <- matrix(nrow = 1, ncol = length(spectrum_x))
        for (i in 1:length(spectrum_x)) {
            spectrum_approx[1, i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), i])
        }

        # Standardization of spectra so that total area equals 1
        spectrum_y_normed <- c()
        spectrum_approx_normed <- c()

        # Standardize the spectra
        spectrum_y_normed <- spectrum_y / sum(spectrum_y)
        spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

        # Calculate the difference between normed original spectrum and normed approximated spectrum
        difference_normed <- c()
        for (i in 1:length(spectrum_x)) {
            difference_normed[i] <- (spectrum_y_normed[i] - spectrum_approx_normed[i])^2
        }
        mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
        message(sprintf("Normed MSE after iteration %d: %.22f", b, mse_normed))
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(lambda_new))
    for (i in 1:length(lambda_new)) {
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
        integrals = integrals,
        lorentz_curves_initial = lorentz_curves_initial
    )
    spec
}

# Private helpers of `find_peaks_v12` #####

calc_second_derivative_v12 <- function(y) {
    n <- length(y)
    x <- c(NA, y[-n]) # x[i] == y[i-1]
    z <- c(y[-1], NA) # z[i] == y[i+1]
    d <- x + z - 2 * y
}

get_right_border_v12 <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (is_right_border) {
            return(r)
        }
        r <- r + 1
    }
    return(NA)
}

get_left_border_v12 <- function(j, d, m) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (is_left_border) {
            return(l)
        }
        l <- l - 1
    }
}

get_peak_score_v12 <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        NA
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}
