#' @title Deconvolute one single spectrum
#' @description Deconvolute one single spectrum
#' @param path Path to file or folder containing the spectra files.
#' @param type Format of the spectra files. Either `"bruker"` or `"jcampdx"`.
#' @param processing_value Processing value for the file. E.g. `"10"`. Called `procno` in the Bruker TopSpin Manual.
#' @param spectroscopy_value Spectroscopy value for the file. E.g. `"10"`. Called `specno` in the Bruker TopSpin Manual.
#' @param number_iterations Number of iterations for the approximation of the parameters for the Lorentz curves.
#' @param ws$hwidth_ppm Half width of the water artefact in ppm.
#' @param signal_free_region Row vector with two entries consisting of the ppm positions for the left and right border of the signal free region of the spectrum.
#' @param smoothing_param Row vector with two entries consisting of the number of smoothing repeats for the whole spectrum and the number of data points (uneven) for the mean calculation.
#' @param delta Threshold value to distinguish between signal and noise.
#' @param scale_factor Row vector with two entries consisting of the factor to scale the x-axis and the factor to scale the y-axis.
#' @param same_parameter TODO
#' @param current_filenumber TODO
#' @param number_of_files TODO
#' @param bwc Use the old, slightly incorrect methods for calculating signal free region and water signal to maintain backwards compatibility with MetaboDecon1D results? For details see `Check: ...` issues in `TODOS.md`.
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
                                    same_parameter = FALSE,
                                    spectroscopy_value = 10,
                                    processing_value = 10,
                                    number_iterations = 10,
                                    range_water_signal_ppm = 0.1527692,
                                    signal_free_region = c(11.44494, -1.8828),
                                    smoothing_param = c(2, 5),
                                    delta = 6.4,
                                    scale_factor = c(1e3, 1e6),
                                    current_filenumber = 1,
                                    number_of_files = 1,
                                    bwc = TRUE) {

    # Parse arguments
    type <- match.arg(type, c("bruker", "jcampdx"))
    sfx <- scale_factor[1]
    sfy <- scale_factor[2]
    expno <- spectroscopy_value
    procno <- processing_value

    # Load spectrum
    spectrum <- switch(type, "bruker" = load_bruker_spectrum_v2(path, sfx, sfy, expno, procno), "jcampdx" = load_jcampdx_spectrum_v2(path, sfx, sfy))
    ws <- calculate_water_signal(spectrum, hwidth_ppm = range_water_signal_ppm, bwc = bwc)
    sfr <- calculate_signal_free_region(spectrum, left_ppm = signal_free_region[1], right_ppm = signal_free_region[2], sfx = sfx, bwc = bwc)

    # Ask the user to confirm/adjust the signal free region and/or water signal
    if (isFALSE(same_parameter) || current_filenumber == 1) {

        # Confirm/adjust the signal free region
        plot_sfr(spectrum, sfr)
        sfr_ok <- get_yn_input("Signal free region borders correct selected? (Area left and right of the green lines)")
        while (!sfr_ok) {
            left_ppm <- get_num_input("Choose another left border: [e.g. 12] ", min = spectrum$ppm_min, max = spectrum$ppm_max)
            right_ppm <- get_num_input("Choose another right border: [e.g. -2] ", min = spectrum$ppm_min, max = spectrum$ppm_max)
            sfr <- calculate_signal_free_region(spectrum, left_ppm, right_ppm, sf = sfx, bwc = bwc)
            plot_sfr(spectrum, sfr)
            sfr_ok <- get_yn_input("Signal free region borders correct selected? (Area left and right of the green lines)")
        }

        # Confirm/adjust water signal
        plot_ws(spectrum, ws)

        # Check for correct range of water artefact
        check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

        # Set parameter to TRUE or FALSE
        if (check_range_water_signal == "y" | check_range_water_signal == "n") {
            correct_input <- TRUE
        } else {
            correct_input <- FALSE
        }

        # Check if User input is correct or not
        while (correct_input == FALSE) {
            # Ask User if he want to use same parameters for all spectra of the folder
            message("Error. Please type only y or n.")
            check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

            if (check_range_water_signal == "y" | check_range_water_signal == "n") {
                correct_input <- TRUE
            } else {
                correct_input <- FALSE
            }
        }


        while (check_range_water_signal == "n") {
            ws$hwidth_ppm <- readline(prompt = "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] ")

            # Check if input is a digit
            digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", ws$hwidth_ppm)

            while (digit_true != TRUE) {
                # Ask User which of the files should be used to adjust the parameters
                message("Error. Please only type a digit.")
                ws$hwidth_ppm <- readline(prompt = "Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154] ")
                # Check if input is a digit
                digit_true <- grepl("[+-]?([0-9]*[.])?[0-9]+", ws$hwidth_ppm)
            }
            # Save as numeric
            ws$hwidth_ppm <- as.numeric(ws$hwidth_ppm)



            # Remove water signal
            ws$center_dp <- length(spectrum$x) / 2
            ws$center_ppm <- spectrum$x_ppm[length(spectrum$x_ppm) / 2]
            # Recalculate ppm into data points
            ws$hwidth_dp <- ws$hwidth_ppm / (spectrum$ppm_nstep)
            ws$left_dp <- ws$center_dp - ws$hwidth_dp
            ws$right_dp <- ws$center_dp + ws$hwidth_dp

            plot(spectrum$x_ppm, spectrum$y, type = "l", xlab = "[ppm]", ylab = "Intensity [a.u.]", xlim = rev(range((ws$center_ppm - 2 * ws$hwidth_ppm), (ws$center_ppm + 2 * ws$hwidth_ppm))))
            graphics::abline(v = spectrum$x_ppm[ws$left_dp], col = "red")
            graphics::abline(v = spectrum$x_ppm[ws$right_dp], col = "red")

            check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

            # Set parameter to TRUE or FALSE
            if (check_range_water_signal == "y" | check_range_water_signal == "n") {
                correct_input <- TRUE
            } else {
                correct_input <- FALSE
            }

            # Check if User input is correct or not
            while (correct_input == FALSE) {
                # Ask User if he want to use same parameters for all spectra of the folder
                message("Error. Please type only y or n.")
                check_range_water_signal <- readline(prompt = "Water artefact fully inside red vertical lines? (y/n): ")

                if (check_range_water_signal == "y" | check_range_water_signal == "n") {
                    correct_input <- TRUE
                } else {
                    correct_input <- FALSE
                }
            }
        }

    # Save adjusted ws$hwidth_dp
    ws$hwidth_ppm <- ws$hwidth_ppm
    }


    # Remove water signal
    for (i in ws$right_dp:ws$left_dp) {
        spectrum$y[i] <- 0.00000001
    }

    # Remove negative values of spectrum by Saving the absolut values
    for (i in 1:length(spectrum$y)) {
        spectrum$y[i] <- abs(spectrum$y[i])
    }

    # Variable Mean Filter
    smoothing_iteration <- smoothing_param[1]
    smoothing_pairs <- smoothing_param[2]

    # Check if number of smoothing pairs is uneven
    while (smoothing_pairs %% 2 == "0") {
        smoothing_pairs <- as.numeric(readline(prompt = "Number of smoothing pairs is even. Please choose uneven number: "))
    }

    for (j in 1:smoothing_iteration) {
        smoothed_spectrum_y <- c()
        for (i in 1:(length(spectrum$x))) {
            # Calculate borders
            left_border <- i - floor(smoothing_pairs / 2)
            right_border <- i + floor(smoothing_pairs / 2)

            # Calculate smoothed spectrum for borders
            if (left_border <= 0) {
                left_border <- 1
                smoothed_spectrum_y[i] <- (1 / right_border) * sum(spectrum$y[left_border:right_border])
            } else if (right_border >= length(spectrum$x)) {
                right_border <- length(spectrum$x)
                smoothed_spectrum_y[i] <- (1 / (right_border - left_border + 1)) * sum(spectrum$y[left_border:right_border])
            } else {
                # Calculate smoothed spectrum
                smoothed_spectrum_y[i] <- (1 / smoothing_pairs) * sum(spectrum$y[left_border:right_border])
            }
        }
        # Save smoothed spectrum
        spectrum$y <- smoothed_spectrum_y
    }



    # Peak selection procedure

    # Calculate second derivative of spectrum
    second_derivative <- matrix(nrow = 2, ncol = length(spectrum$x) - 2)
    for (i in 2:length(spectrum$x) - 1) {
        second_derivative[1, i - 1] <- spectrum$x[i]
        second_derivative[2, i - 1] <- spectrum$y[i - 1] + spectrum$y[i + 1] - 2 * spectrum$y[i]
    }

    # Find all local minima of second derivative
    peaks_x <- c()
    peaks_index <- c()
    second_derivative_border <- ncol(second_derivative) - 1
    for (i in 2:second_derivative_border) {
        if (second_derivative[2, i] < 0) {
            if ((second_derivative[2, i] <= second_derivative[2, i - 1]) & (second_derivative[2, i] < second_derivative[2, i + 1])) {
                # if(((spectrum$y[i+1] >= spectrum$y[i]) & (spectrum$y[i+1] > spectrum$y[i+2])) | ((spectrum$y[i+1] > spectrum$y[i]) & (spectrum$y[i+1] >= spectrum$y[i+2]))){
                # Add local minima to peak list
                peaks_x <- c(peaks_x, second_derivative[1, i])
                peaks_index <- c(peaks_index, i)
            }
        }
    }



    # Find all left positions of all local minima of second derivative
    left_position <- matrix(nrow = 1, ncol = length(peaks_x))
    for (i in 1:length(peaks_x)) {
        # Save next left position of current local minima
        next_left <- peaks_index[i] + 1
        while ((peaks_index[i] < next_left) & (next_left < ncol(second_derivative))) {
            if (second_derivative[2, next_left - 1] < second_derivative[2, next_left]) {
                if (((second_derivative[2, next_left - 1] < second_derivative[2, next_left]) & (second_derivative[2, next_left + 1] <= second_derivative[2, next_left])) | ((second_derivative[2, next_left] < 0) & (second_derivative[2, next_left + 1] >= 0))) {
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
    right_position <- matrix(nrow = 1, ncol = length(peaks_x))
    for (i in 1:length(peaks_x)) {
        # Save next right position of current local minima
        next_right <- peaks_index[i] - 1
        while ((next_right < peaks_index[i]) & (next_right >= 2)) {
            if (second_derivative[2, next_right + 1] < second_derivative[2, next_right]) {
                if (((second_derivative[2, next_right + 1] < second_derivative[2, next_right]) & (second_derivative[2, next_right - 1] <= second_derivative[2, next_right])) | ((second_derivative[2, next_right] < 0) & (second_derivative[2, next_right - 1] >= 0))) {
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
            peaks_x <- peaks_x[-i]
            peaks_index <- peaks_index[-i]
            left_position <- left_position[-i]
            right_position <- right_position[-i]
        }
    }

    # Calculate peak triplet score to distinguish between signal and noise
    scores <- matrix(nrow = 1, ncol = length(peaks_x))
    scores_left <- matrix(nrow = 1, ncol = length(peaks_x))
    scores_right <- matrix(nrow = 1, ncol = length(peaks_x))
    for (i in 1:length(peaks_x)) {
        # Calculate left score
        left_score <- 0
        for (j in peaks_index[i]:left_position[i]) {
            left_score <- sum(left_score, abs(second_derivative[2, j]))
        }
        scores_left[i] <- left_score
        # Calculate right score
        right_score <- 0
        for (k in right_position[i]:peaks_index[i]) {
            right_score <- sum(right_score, abs(second_derivative[2, k]))
        }
        scores_right[i] <- right_score
        # Save minimum score
        scores[i] <- min(left_score, right_score)
    }

    # Calculate mean of the score and standard deviation of the score of the signal free region R
    index_left <- which(spectrum$x[peaks_index + 1] >= sfrl_sdp)
    index_right <- which(spectrum$x[peaks_index + 1] <= sfrr_sdp)

    mean_score <- mean(c(scores[index_left], scores[index_right]))
    sd_score <- stats::sd(c(scores[index_left], scores[index_right]))

    # Filter peak triplets
    filtered_peaks <- c()
    filtered_left_position <- c()
    filtered_right_position <- c()
    save_scores <- c()
    for (i in 1:length(peaks_x)) {
        if (scores[i] >= mean_score + delta * sd_score) {
            # Save peak position
            filtered_peaks <- c(filtered_peaks, peaks_index[i])
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
        w_1 <- c(w_1, spectrum$x[filtered_left_position[i] + 1])
        w_2 <- c(w_2, spectrum$x[filtered_peaks[i] + 1])
        w_3 <- c(w_3, spectrum$x[filtered_right_position[i] + 1])

        # Calculate intensity of peak triplets
        y_1 <- c(y_1, spectrum$y[filtered_left_position[i] + 1])
        y_2 <- c(y_2, spectrum$y[filtered_peaks[i] + 1])
        y_3 <- c(y_3, spectrum$y[filtered_right_position[i] + 1])

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
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spectrum$x))
    for (i in 1:length(filtered_peaks)) {
        # If A = 0, then the lorentz curve is a zero line
        if (A[i] == 0) {
            lorentz_curves_initial[i, ] <- 0
        } else {
            lorentz_curves_initial[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (spectrum$x - w[i])^2)))
        }
    }


    # Approximation of lorentz curves
    for (b in 1:number_iterations) {
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
            w_1_new <- c(w_1_new, spectrum$x[filtered_left_position[i] + 1])
            w_2_new <- c(w_2_new, spectrum$x[filtered_peaks[i] + 1])
            w_3_new <- c(w_3_new, spectrum$x[filtered_right_position[i] + 1])

            # Calculate the sum of all lorentz curves for each data point
            sum_left[i] <- sum(lorentz_curves_initial[1:length(filtered_left_position), filtered_left_position[i] + 1])
            sum_peaks[i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), filtered_peaks[i] + 1])
            sum_right[i] <- sum(lorentz_curves_initial[1:length(filtered_right_position), filtered_right_position[i] + 1])

            # Calculate the proprotion between original spectrum an the sum of the lorentz curves for each peak triplets position
            proportion_left[i] <- spectrum$y[filtered_left_position[i] + 1] / sum_left[i]
            proportion_peaks[i] <- spectrum$y[filtered_peaks[i] + 1] / sum_peaks[i]
            proportion_right[i] <- spectrum$y[filtered_right_position[i] + 1] / sum_right[i]

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
                lorentz_curves_initial[i, ] <- abs(A_new[i] * (lambda_new[i] / (lambda_new[i]^2 + (spectrum$x - w_new[i])^2)))
            }
        }

        # Calculate sum of lorentz curves
        spectrum_approx <- matrix(nrow = 1, ncol = length(spectrum$x))
        for (i in 1:length(spectrum$x)) {
            spectrum_approx[1, i] <- sum(lorentz_curves_initial[1:length(filtered_peaks), i])
        }
        # ToSc: use vectorized functions, e.g.
        # spectrum_approx <- colSums(lorentz_curves_initial[1:length(filtered_peaks), ])

        # Standardize the spectra so that total area equals 1
        spectrum_y_normed <- spectrum$y / sum(spectrum$y)
        spectrum_approx_normed <- spectrum_approx / sum(spectrum_approx)

        # Calculate the difference between normed original spectrum and normed approximated spectrum
        difference_normed <- c()
        for (i in 1:length(spectrum$x)) {
            difference_normed[i] <- (spectrum_y_normed[i] - spectrum_approx_normed[i])^2
        }
        mse_normed <- (1 / length(difference_normed)) * sum(difference_normed)
        message(paste("\nNormed MSE value of iteration", b, "is: "))
        print(mse_normed)
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(lambda_new))
    for (i in 1:length(lambda_new)) {
        integrals[1, i] <- A_new[i] * (atan((-w_new[i] + (spectrum$length / sfx)) / lambda_new[i]) - atan((-w_new[i]) / lambda_new[i]))
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
        peak_triplets_middle[i] <- spectrum$x_ppm[index_peak_triplets_middle[i]]
        peak_triplets_left[i] <- spectrum$x_ppm[index_peak_triplets_left[i]]
        peak_triplets_right[i] <- spectrum$x_ppm[index_peak_triplets_right[i]]
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

    return_list <- list("filename" = name, "spectrum_x" = spectrum$x, "spectrum_x_ppm" = spectrum$x_ppm, "spectrum_y" = spectrum$y, "lorentz_curves" = lorentz_curves_initial, "mse_normed" = mse_normed, "spectrum_approx" = spectrum_approx, "index_peak_triplets_middle" = index_peak_triplets_middle, "index_peak_triplets_left" = index_peak_triplets_left, "index_peak_triplets_right" = index_peak_triplets_right, "peak_triplets_middle" = peak_triplets_middle, "peak_triplets_left" = peak_triplets_left, "peak_triplets_right" = peak_triplets_right, "integrals" = integrals, "signal_free_region" = c(sfrl_sdp, sfrr_sdp), "ws$hwidth_ppm" = ws$hwidth_ppm, "A" = A_new, "lambda" = lambda_new, "w" = w_new)
    return(return_list)
}


#' @title Calculate water signal parameters
#' @description Calculates water signal parameters for a given spectrum.
#' @param spectrum A list representing the spectrum.
#' @param hwidth_ppm The half width in ppm. Default is range_water_signal_ppm.
#' @param bwc Use the old, slightly incorrect methods for calculating water signal values to maintain backwards compatibility with MetaboDecon1D results? For details see issue `Check: water signal calculation` in `TODOS.md`.
#' @return List of parameters including half width in dp and ppm, center line in dp and ppm and right and left borders in dp and ppm.
calculate_water_signal <- function(spectrum, hwidth_ppm, bwc = TRUE) {
    ws <- list()
    if (bwc) {
        ws$hwidth_dp <- hwidth_ppm / (spectrum$ppm_nstep) # half width in dp
        ws$hwidth_ppm <- hwidth_ppm # half width in ppm
        ws$center_dp <- spectrum$n / 2 # center line in dp
        ws$center_ppm <- spectrum$x_ppm[ws$center_dp] # center in ppm
        ws$right_dp <- ws$center_dp + ws$hwidth_dp # right border in dp
        ws$right_ppm <- spectrum$x_ppm[ws$right_dp] # right border in ppm
        ws$left_dp <- ws$center_dp - ws$hwidth_dp # left border in dp
        ws$left_ppm <- spectrum$x_ppm[ws$left_dp] # left border in ppm
    } else {
        ws$hwidth_dp <- hwidth_ppm / (spectrum$ppm_step) # ppm_step instead of ppm_nstep
        ws$hwidth_ppm <- hwidth_ppm
        ws$center_dp <- (spectrum$n - 1) / 2
        ws$center_ppm <- dp_to_ppm(ws$center_dp, spectrum)
        ws$right_dp <- ws$center_dp + ws$hwidth_dp
        ws$right_ppm <- dp_to_ppm(ws$right_dp, spectrum)
        ws$left_dp <- ws$center_dp - ws$hwidth_dp
        ws$left_ppm <- dp_to_ppm(ws$left_dp, spectrum)
    }
    ws
}

calculate_signal_free_region <- function(spectrum, left_ppm, right_ppm, sf = 1000, bwc = TRUE) {
    sfr <- list()
    sfr$left_ppm <- left_ppm
    sfr$left_dp <- ppm_to_dp(left_ppm, spectrum, bwc)
    sfr$left_sdp <- sfr$left_dp / sf
    sfr$right_ppm <- right_ppm
    sfr$right_dp <- ppm_to_dp(right_ppm, spectrum, bwc)
    sfr$right_sdp <- sfr$right_dp / sf
    sfr
}

#' @title Plot Signal Free Region
#' @description Draws the signal free region as green vertical lines into the given spectrum.
#' @param spectrum A list representing the spectrum as returned by [load_jcampdx_spectrum()] or [load_bruker_spectrum()].
#' @param sfr A list representing the signal free region as returned by [calculate_signal_free_region()].
#' @return NULL. Called for side effect of plotting the signal free region.
plot_sfr <- function(spectrum, sfr) {
    plot(
        x = spectrum$x_ppm,
        y = spectrum$y,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = c(spectrum$ppm_max, spectrum$ppm_min)
    )
    graphics::abline(v = sfr$left_ppm, col = "green")
    graphics::abline(v = sfr$right_ppm, col = "green")
}

#' @title Plot Water Signal
#' @description Draws the water signal as red vertical lines into the given spectrum.
#' @param spectrum A list representing the spectrum as returned by [load_jcampdx_spectrum()] or [load_bruker_spectrum()].
#' @param ws A list representing the water signal as returned by [calculate_water_signal()].
#' @return NULL. Called for side effect of plotting the water signal.
plot_ws <- function(spectrum, ws) {
    plot(
        spectrum$x_ppm,
        spectrum$y,
        type = "l",
        xlab = "[ppm]",
        ylab = "Intensity [a.u.]",
        xlim = rev(range(
            ws$center_ppm - 2 * ws$hwidth_ppm,
            ws$center_ppm + 2 * ws$hwidth_ppm
        ))
    )
    graphics::abline(v = spectrum$x_ppm[ws$left_dp], col = "red")
    graphics::abline(v = spectrum$x_ppm[ws$right_dp], col = "red")
}