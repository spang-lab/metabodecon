# Private API Functions #####

#' @title Generate Lorentz Curves from NMR Spectra
#' @description Deconvolutes NMR spetra and generates a Lorentz curve for each detected signal within a spectra.
#' @param data_path (string) Path to the folder where the original spectra are stored. After deconvolution this folder contains two additional .txt files for each spectrum which contain the spectrum approximated from all deconvoluted signals and a parameter file that contains all numerical values of the deconvolution.
#' @param file_format (string) Format of the spectra files.
#' @param make_rds (bool) Store results as a rds file on disk? Should be set to TRUE if many spectra are evaluated to decrease computation time.
#' @param number_iterations (int) Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param range_water_signal_ppm (float) Half width of the water artefact in ppm.
#' @param signal_free_region (float) Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smoothing_param (int) Row vector with two entries consisting of the number of smoothing repeats for the whole spectrum and the number of data points (uneven) for the mean calculation.
#' @param delta (float) Threshold value to distinguish between signal and noise.
#' @param scale_factor (int) Row vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param ask (bool) Whether to ask for user input during the deconvolution process. If set to FALSE, the provided default values will be used.
#' @details First, an automated curvature based signal selection is performed. Each signal is represented by 3 data points to allow the determination of initial Lorentz curves. These Lorentz curves are then iteratively adjusted to optimally approximate the measured spectrum. For each spectrum two text files will be created in the parent folder i.e. the folder given in data path. The spectrum approximated from all deconvoluted signals and a parameter file that contains all numerical values of the deconvolution. Furthermore, the numerical values of the deconvolution are also stored in a data_frame.
#' @details Shall replace `generate_lorentz_curves` as soon as implementation is finished.
#' @noRd
generate_lorentz_curves_v20 <- function(data_path,
                                       file_format = c("bruker", "jcampdx"),
                                       make_rds = FALSE,
                                       number_iterations = 10,
                                       range_water_signal_ppm = 0.1527692,
                                       signal_free_region = c(11.44494, -1.8828),
                                       smoothing_param = c(2, 5),
                                       delta = 6.4,
                                       scale_factor = c(1000, 1000000),
                                       ask = TRUE) {
    # Check arguments
    file_format <- match.arg(file_format)

    # Switch to data directory
    data_path <- normalizePath(data_path)
    owd <- getwd()
    setwd(data_path)
    on.exit(setwd(owd))

    # Get input files
    if (file_format == "jcampdx") {
        files <- dir(data_path, pattern = "\\.dx$") # `.dx` files inside `data_path`
        spectroscopy_value <- NULL
        processing_value <- NULL
    } else if (file_format == "bruker") {
        files <- list.dirs(data_path, recursive = FALSE, full.names = FALSE) # folders inside `data_path`
        spectroscopy_value <- readline(prompt = "What is the name of the subfolder of your filepath: (e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10) ")
        processing_value <- readline(prompt = "What is the name of the subsubsubfolder of your filepath: (e.g. 10 for C:/Users/Username/Desktop/spectra_folder/spectrum_name/10/pdata/10) ")
    }

    # Reorder files in case user wants to use one spectrum to determine parameters for all spectra
    same_parameter <- get_yn_input("Do you want to use the same parameters (signal_free_region, range_water_signal_ppm) for all spectra?")
    if (same_parameter) {
        print(files)
        number <- get_num_input("Choose number of file which is used to adjust all parameters: [e.g. 1] ", min = 1, max = length(files), int = TRUE)
        message(paste("The selected file to adjust all parameters for all spectra is: ", files[number]))
        files <- c(files[number], files[-number])
    }

    # Do actual deconvolution
    spectrum_data <- list()
    for (i in seq_along(files)) {
        name <- files[i] # bruker: `urine_2`, jcampdx: `urine_2.dx`
        filepath <- switch(file_format, # see [FAQ](../vignettes/FAQ.Rmd#file-structure) for example file structures
            "bruker" = paste(data_path, name, spectroscopy_value, sep = "/"),
            "jcampdx" = data_path,
            stop("Invalid file format")
        )
        x <- deconvolution_v20(filepath, name, file_format, same_parameter, processing_value, number_iterations, range_water_signal_ppm, signal_free_region, smoothing_param, delta, scale_factor, current_filenumber = i, number_of_files = length(files))
        spectrum_data[[name]] <- x
        # Save `range_water_signal` and `signal_free_region` for next loop passage as those might have been adjusted interactively by the user
        range_water_signal_ppm <- x$range_water_signal_ppm
        signal_free_region <- x$signal_free_region
    }

    # Save results
    if (make_rds) {
        saveRDS(object = spectrum_data, file = file.path(data_path, "spectrum_data.rds"))
    }
    return(spectrum_data)
}

#' @title Deconvolute one single spectrum
#' @description Deconvolutes one single spectrum
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
deconvolute_spectrum_v20 <- function(path = file.path(download_example_datasets(), "bruker/urine/urine_1"),
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

    toscutil::stub(deconvolute_spectrum_v2, ask = FALSE)
    type <- match.arg(type, c("bruker", "jcampdx"))

    # Implemented
    spec <- read_spectrum(path, type, sf, expno, procno)
    spec <- determine_signal_free_region_v20(spec, sfr, ask)
    spec <- determine_water_signal_v20(spec, hwidth_ppm = wshw, bwc, ask)
    spec <- remove_water_signal_v20(spec, bwc)
    spec <- remove_negative_signals_v20(spec)
    spec <- smooth_signals(spec, reps = smopts[1], k = smopts[2], bwc)
    spec <- find_peaks_v20(spec)
    spec <- rm_peaks_with_low_scores(spec, delta)
    spec <- init_lorentz_curves_v1(spec)

    # In progress
    spec <- refine_lorentz_curves_v1(spec, nfit)

    # To be done
    spec <- calculate_lorentz_curve_integrals(spec)

    compare_spectra(spec)
    plot_peaks(spec)

    ret <- create_return_list(spec)
}

# Private Helpers #####

#' @title Determine Signal Free Region
#' @description Determines the signal free region (SFR) of a given spectrum. It asks the user to confirm the left and right borders of the SFR, and allows them to adjust these borders if necessary. The function returns a list containing the left and right borders in both ppm and data points (dp), as well as the scaled data points (sdp).
#' @param spec A list representing the spectrum, which should include the minimum and maximum ppm (`$ppm_min` and `$ppm_max`), and the scaling factor (`$sfx`).
#' @param sfr Initial values for the left and right borders of the SFR in ppm. If not provided, the function will ask the user to select the borders.
#' @param ask Logical. If TRUE, the function will ask the user to confirm or adjust the borders of the SFR. Default is TRUE.
#' @param bwc Use the old, slightly incorrect method for conversion from ppm to data points to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: ppm to dp conversion` in TODOS.md
#' @return A list containing the left and right borders of the SFR in ppm (`$left_ppm` and `$right_ppm`), data points (`$left_dp` and `$right_dp`), and scaled data points (`$left_sdp` and `$right_sdp`).
#' @noRd
determine_signal_free_region_v20 <- function(spec, sfr = NULL, ask = TRUE, bwc = TRUE) {
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
determine_water_signal_v20 <- function(spec, hwidth_ppm, bwc = TRUE, ask = TRUE) {
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

remove_water_signal_v20 <- function(spec, bwc = TRUE) {
    if (is.null(spec$ws)) {
        stop("No water signal parameters found. Please call `determine_water_signal_v20()`.")
    }
    y <- spec$y_scaled
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
    spec$y_nows <- y
    spec
}

remove_negative_signals_v20 <- function(spec) {
    if (is.null(spec$y_nows)) stop("Water signal not removed yet. Please call `remove_water_signal_v20()` first.")
    spec$y_pos <- abs(spec$y_nows)
    spec
}

#' @inherit find_peaks_v20
#' @param details Successor of `select_peaks_v0`, `find_left_positions_v0` and `find_right_positions_v0`. The new function `find_peaks_v20()` is a combination of these three functions with a corrected naming convention: what was incorrectly referred to as "left" is now correctly called "right" and vice versa.
#' @noRd
find_peaks_v20 <- function(spec) {
    d <- spec$d <- calc_second_derivative_v20(y = spec$y_smooth, bwc = bwc)
    a <- abs(d)
    m <- length(d)
    dl <- c(NA, d[-m]) # dl[i] == d[i-1]
    dr <- c(d[-1], NA) # dr[i] == d[i+1]
    center <- which(d < 0 & d <= dl & d < dr)
    spec$peak <- data.frame(left = NA, center = center, right = NA, score = NA)
    for (i in seq_along(center)) {
        j <- center[i]
        r <- spec$peak$right[i] <- get_right_border_v20(j, d, m)
        l <- spec$peak$left[i] <- get_left_border_v20(j, d, m)
        s <- spec$peak$score[i] <- get_peak_score_v20(j, l, r, a)
    }
    return(spec)
}

#' @inherit smooth_signals_v1
#' @details New and fast version for smoothing of signals. Implements the same algorithm as `smooth_signal_v1` using different R functions (e.g. [stats::filter()]), causing a massive speedup but also numeric differences compared to the old version.
#' @noRd
smooth_signals_v20 <- function(spec, reps = 2, k = 5) {
    if (k %% 2 == 0) stop("k must be odd")

    Z <- vector("list", length = reps)
    y <- spec$y_pos
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
    spec$y_smooth <- Z[[reps]]
    spec
}

init_lorentz_curves_v20 <- function(spec) {
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
    save_scores <- spec$peak$score[spec$peak$high]

    # Calculate parameters w, lambda and A for the initial lorentz curves
    for (i in 1:length(filtered_peaks)) {
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

refine_lorentz_curves_v20 <- function(spec, nfit) {

    spectrum_x <- spec$sdp
    spectrum_y <- spec$y_smooth
    filtered_peaks <- as.integer(spec$peak$center[spec$peak$high] - 1)
    filtered_left_position <- spec$peak$right[spec$peak$high] - 1
    filtered_right_position <- spec$peak$left[spec$peak$high] - 1
    save_scores <- spec$peak$score[spec$peak$high]

    w <- spec$lc$w
    w_delta <- spec$lc$w_delta
    lambda <- spec$lc$lambda
    A <- spec$lc$A
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spectrum_x))
    for (i in 1:length(filtered_peaks)) {
        lorentz_curves_initial[i, ] <- if (A[i] == 0) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (spectrum_x - w[i])^2)))
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

    spec$lc_refined$A <- A_new
    spec$lc_refined$lambda <- lambda_new
    spec$lc_refined$w <- w_new
    spec$lc_refined$w_delta <- w_delta_new
    spec$lc_refined$spectrum_approx <- spectrum_approx
    spec$lc_refined$spectrum_y_normed <- spectrum_y_normed
    spec$lc_refined$spectrum_approx_normed <- spectrum_approx_normed
    spec$lc_refined$difference_normed <- difference_normed
    spec$lc_refined$mse_normed <- mse_normed
    spec
}

calculate_lorentz_curve_integrals_v20 <- function() {
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
}

create_return_list_v20 <- function(spec) {
    return_list <- list(
        filename = name,
        spectrum_x = spec$sdp,
        spectrum_x_ppm = spec$x_ppm,
        spectrum_y = spec$y,
        lorentz_curves = lorentz_curves_initial,
        mse_normed = mse_normed,
        spectrum_approx = spectrum_approx,
        index_peak_triplets_middle = index_peak_triplets_middle,
        index_peak_triplets_left = index_peak_triplets_left,
        index_peak_triplets_right = index_peak_triplets_right,
        peak_triplets_middle = peak_triplets_middle,
        peak_triplets_left = peak_triplets_left,
        peak_triplets_right = peak_triplets_right,
        integrals = integrals,
        sfr = c(sfrl_sdp, sfrr_sdp),
        wshwidth_ppm = ws$hwidth_ppm,
        A = A_new,
        lambda = lambda_new,
        w = w_new
    )
}

# Private Helpers of find_peaks_v20 #####

calc_second_derivative_v20 <- function(y, bwc) {
    n <- length(y)
    if (bwc) {
        x <- c(NA, y[-n]) # x[i] == y[i-1]
        z <- c(y[-1], NA) # z[i] == y[i+1]
        d <- x + z - 2 * y
    } else {
        # Using diff is almost equivalent to the above, but due to numeric instabilities, it sometimes gives slightly different results (e.g. -5.51600000000001e-05 instead of -5.51599999999998e-05). Since we do `<=` and `<` comparisons further below this is not backwards compatible. Nonetheless, we keep it here, because it would make the code slightly more readable (runtime is almost the same).
        d <- c(NA, diff(y, differences = 2), NA)
    }
}

get_right_border_v20 <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (is_right_border) return(r)
        r <- r + 1
    }
    return(NA)
}

get_left_border_v20 <- function(j, d, m) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (is_left_border) return(l)
        l <- l - 1
    }
}

get_peak_score_v20 <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        NA
     } else {
        min(sum(a[l:j]), sum(a[j:r]))
     }
}

rescale <- function(x, range = c(0, 1)) {
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    x_len <- x_max - x_min
    y_len <- max(range) - min(range)
    x_percent <- (x - x_min) / x_len
    y <- min(range) + x_percent * y_len
    y
}