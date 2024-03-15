init_lorentz_curves_v0 <- function(spec) {

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
}
