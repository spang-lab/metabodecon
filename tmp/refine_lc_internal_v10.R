refine_lc_internal_v10 <- function(x, y, pl, pc, pr, Yt) {

    cat3("Refining Lorentz Curves")

    spectrum_x <- x
    spectrum_y <- y
    filtered_left_position <- pl - 1
    filtered_peaks <- pc - 1
    filtered_right_position <- pr - 1
    lorentz_curves_initial <- Yt

    # Calculate new heights of peak triplets
    w_1_new <- c(); w_2_new <- c(); w_3_new <- c()
    y_1_new <- c(); y_2_new <- c(); y_3_new <- c()
    w_1_2_new <- c(); w_1_3_new <- c(); w_2_3_new <- c()
    y_1_2_new <- c(); y_1_3_new <- c(); y_2_3_new <- c()
    w_delta_new <- c()
    w_new <- c()
    lambda_new <- c()
    A_new <- c()
    sum_left <- c(); sum_peaks <- c(); sum_right <- c()
    proportion_left <- c(); proportion_peaks <- c(); proportion_right <- c()
    w_1_bak1 <- c(); w_2_bak1 <- c(); w_3_bak1 <- c() # DEBUG
    w_1_bak2 <- c(); w_2_bak2 <- c(); w_3_bak2 <- c() # DEBUG
    w_1_bak3 <- c(); w_2_bak3 <- c(); w_3_bak3 <- c() # DEBUG
    y_1_bak1 <- c(); y_2_bak1 <- c(); y_3_bak1 <- c() # DEBUG
    y_1_bak2 <- c(); y_2_bak2 <- c(); y_3_bak2 <- c() # DEBUG
    y_1_bak3 <- c(); y_2_bak3 <- c(); y_3_bak3 <- c() # DEBUG
    y_1_bak4 <- c(); y_2_bak4 <- c(); y_3_bak4 <- c() # DEBUG

    for (i in seq_along(filtered_peaks)) {
        # Calculate the position of the peak triplets

        w_1_new <- c(w_1_new, spectrum_x[filtered_left_position[i] + 1])
        w_2_new <- c(w_2_new, spectrum_x[filtered_peaks[i] + 1])
        w_3_new <- c(w_3_new, spectrum_x[filtered_right_position[i] + 1])
        w_1_bak1[i] <- c(w_1_new[i]) # DEBUG
        w_2_bak1[i] <- c(w_2_new[i]) # DEBUG
        w_3_bak1[i] <- c(w_3_new[i]) # DEBUG

        # Calculate the sum of all lorentz curves for each data point
        sum_left[i] <- sum(lorentz_curves_initial[seq_along(filtered_left_position), filtered_left_position[i] + 1])
        sum_peaks[i] <- sum(lorentz_curves_initial[seq_along(filtered_peaks), filtered_peaks[i] + 1])
        sum_right[i] <- sum(lorentz_curves_initial[seq_along(filtered_right_position), filtered_right_position[i] + 1])

        # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
        proportion_left[i] <- spectrum_y[filtered_left_position[i] + 1] / sum_left[i]
        proportion_peaks[i] <- spectrum_y[filtered_peaks[i] + 1] / sum_peaks[i]
        proportion_right[i] <- spectrum_y[filtered_right_position[i] + 1] / sum_right[i]
        y_1_bak1[i] <- spectrum_y[filtered_left_position[i] + 1] # DEBUG
        y_2_bak1[i] <- spectrum_y[filtered_peaks[i] + 1] # DEBUG
        y_3_bak1[i] <- spectrum_y[filtered_right_position[i] + 1] # DEBUG

        # Calculate the new heights of the peak triplets
        y_1_new[i] <- lorentz_curves_initial[i, filtered_left_position[i] + 1] * proportion_left[i]
        y_2_new[i] <- lorentz_curves_initial[i, filtered_peaks[i] + 1] * proportion_peaks[i]
        y_3_new[i] <- lorentz_curves_initial[i, filtered_right_position[i] + 1] * proportion_right[i]
        y_1_bak2[i] <- y_1_new[i] # DEBUG
        y_2_bak2[i] <- y_2_new[i] # DEBUG
        y_3_bak2[i] <- y_3_new[i] # DEBUG

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
        w_1_bak2[i] <- w_1_new[i] # DEBUG
        w_2_bak2[i] <- w_2_new[i] # DEBUG
        w_3_bak2[i] <- w_3_new[i] # DEBUG
        y_1_bak3[i] <- y_1_new[i] # DEBUG
        y_2_bak3[i] <- y_2_new[i] # DEBUG
        y_3_bak3[i] <- y_3_new[i] # DEBUG

        # Move triplet to zero position
        w_delta_new[i] <- w_1_new[i]
        w_1_new[i] <- w_1_new[i] - w_delta_new[i]
        w_2_new[i] <- w_2_new[i] - w_delta_new[i]
        w_3_new[i] <- w_3_new[i] - w_delta_new[i]
        w_1_bak3[i] <- w_1_new[i] # DEBUG
        w_2_bak3[i] <- w_2_new[i] # DEBUG
        w_3_bak3[i] <- w_3_new[i] # DEBUG

        # Calculate difference of peak triplet positions
        w_1_2_new <- c(w_1_2_new, w_1_new[i] - w_2_new[i])
        w_1_3_new <- c(w_1_3_new, w_1_new[i] - w_3_new[i])
        w_2_3_new <- c(w_2_3_new, w_2_new[i] - w_3_new[i])

        # Calculate difference of new intensity values of peak triplets
        y_1_2_new <- c(y_1_2_new, y_1_new[i] - y_2_new[i])
        y_1_3_new <- c(y_1_3_new, y_1_new[i] - y_3_new[i])
        y_2_3_new <- c(y_2_3_new, y_2_new[i] - y_3_new[i])
        y_1_bak4[i] <- y_1_new[i] # DEBUG
        y_2_bak4[i] <- y_2_new[i] # DEBUG
        y_3_bak4[i] <- y_3_new[i] # DEBUG

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
    cat3(sprintf("Normed MSE: %.22f", mse_normed))

    return(list(
        w_1_bak1 = w_1_bak1, w_2_bak1 = w_2_bak1, w_3_bak1 = w_3_bak1, # DEBUG
        w_1_bak2 = w_1_bak2, w_2_bak2 = w_2_bak2, w_3_bak2 = w_3_bak2, # DEBUG
        w_1_bak3 = w_1_bak3, w_2_bak3 = w_2_bak3, w_3_bak3 = w_3_bak3, # DEBUG
        y_1_bak1 = y_1_bak1, y_2_bak1 = y_2_bak1, y_3_bak1 = y_3_bak1, # DEBUG
        y_1_bak2 = y_1_bak2, y_2_bak2 = y_2_bak2, y_3_bak2 = y_3_bak2, # DEBUG
        y_1_bak3 = y_1_bak3, y_2_bak3 = y_2_bak3, y_3_bak3 = y_3_bak3, # DEBUG
        y_1_bak4 = y_1_bak4, y_2_bak4 = y_2_bak4, y_3_bak4 = y_3_bak4, # DEBUG
        w_1_new = w_1_new, w_2_new = w_2_new, w_3_new = w_3_new,
        y_1_new = y_1_new, y_2_new = y_2_new, y_3_new = y_3_new,
        w_1_2_new = w_1_2_new, w_1_3_new = w_1_3_new, w_2_3_new = w_2_3_new,
        y_1_2_new = y_1_2_new, y_1_3_new = y_1_3_new, y_2_3_new = y_2_3_new,
        w_delta_new = w_delta_new,
        w_new = w_new,
        lambda_new = lambda_new,
        A_new = A_new,
        sum_left = sum_left, sum_peaks = sum_peaks, sum_right = sum_right,
        proportion_left = proportion_left, proportion_peaks = proportion_peaks, proportion_right = proportion_right,
        lorentz_curves_initial = lorentz_curves_initial,
        spectrum_approx = spectrum_approx,
        mse_normed = mse_normed
    ))
}
