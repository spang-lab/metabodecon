init_lorentz_curves_v2 <- function(spec) {

    h <- which(spec$peak$high)
    l <- spec$peak$left[h]
    z <- spec$peak$center[h] # use letter z instead c for center, because c is an important function in R
    r <- spec$peak$right[h]
    x <- spec$sdp
    y <- spec$Y$smooth

    # Calculate parameters w, lambda and A for the initial lorentz curves
    xl <- x[l] # called w_1 in the [algorithm paper](https://doi.org/10.1016/j.jmr.2009.09.003)
    xz <- x[z] # called w_2 ...
    xr <- x[r] # called w_3 ...
    yl <- y[l]
    yz <- y[z]
    yr <- y[r]

    # Calculate mirrored points for ascending shoulders if necessary
    i <- which((yl < yz) & (yz < yr))
    xr[i] <- 2 * xz[i] - xl[i]
    yr[i] <- yl[i]

    for (i in 1:length(z)) {
        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if ((yl[i] < yz[i]) & (yz[i] < yr[i])) {
            xr[i] <- 2 * xz[i] - xl[i]
            yr[i] <- yl[i]
        }
    }
    for (i in 1:length(z)) {
        # Calculate mirrored points if necesccary
        # For descending shoulders
        if ((yl[i] > yz[i]) & (yz[i] > yr[i])) {
            xl[i] <- 2 * xz[i] - xr[i]
            yl[i] <- yr[i]
        }
    }
    for (i in 1:length(z)) {

        # Move triplet to zero position
        w_delta[i] <- xl[i]
        xl[i] <- xl[i] - w_delta[i]
        xz[i] <- xz[i] - w_delta[i]
        xr[i] <- xr[i] - w_delta[i]

        # Calculate difference of position of peak triplets
        w_1_2 <- c(w_1_2, xl[i] - xz[i])
        w_1_3 <- c(w_1_3, xl[i] - xr[i])
        w_2_3 <- c(w_2_3, xz[i] - xr[i])

        # Calculate difference of intensity values of peak triplets
        y_1_2 <- c(y_1_2, yl[i] - yz[i])
        y_1_3 <- c(y_1_3, yl[i] - yr[i])
        y_2_3 <- c(y_2_3, yz[i] - yr[i])

        # Calculate w for each peak triplet
        w_result <- (xl[i]^2 * yl[i] * y_2_3[i] +
            xr[i]^2 * yr[i] * y_1_2[i] +
            xz[i]^2 * yz[i] * (-y_1_3[i])
        ) / (
            2 * w_1_2[i] * yl[i] * yz[i] -
                2 * (w_1_3[i] * yl[i] + (-w_2_3[i]) * yz[i]) * yr[i]
        )
        w_result <- w_result + w_delta[i]
        w <- c(w, w_result)
        # Wenn y Werte nach der Hoehenanpassung 0 werden, so ist w_new[i] NaN
        if (is.nan(w[i])) {
            w[i] <- 0
        }

        # Calculate lambda for each peak triplet
        lambda_result <- -((sqrt(abs((-xz[i]^4 * yz[i]^2 * y_1_3[i]^2 - xl[i]^4 * yl[i]^2 * y_2_3[i]^2 - xr[i]^4 * y_1_2[i]^2 * yr[i]^2 + 4 * xz[i] * xr[i]^3 * yz[i] * ((-yl[i]) + yz[i]) * yr[i]^2 + 4 * xz[i]^3 * xr[i] * yz[i]^2 * yr[i] * ((-yl[i]) + yr[i]) + 4 * xl[i]^3 * yl[i]^2 * y_2_3[i] * (xz[i] * yz[i] - xr[i] * yr[i]) + 4 * xl[i] * yl[i] * (xz[i]^3 * yz[i]^2 * y_1_3[i] - xz[i] * xr[i]^2 * yz[i] * (yl[i] + yz[i] - 2 * yr[i]) * yr[i] + xr[i]^3 * y_1_2[i] * yr[i]^2 - xz[i]^2 * xr[i] * yz[i] * yr[i] * (yl[i] - 2 * yz[i] + yr[i])) + 2 * xz[i]^2 * xr[i]^2 * yz[i] * yr[i] * (yl[i]^2 - 3 * yz[i] * yr[i] + yl[i] * (yz[i] + yr[i])) + 2 * xl[i]^2 * yl[i] * (-2 * xz[i] * xr[i] * yz[i] * yr[i] * (-2 * yl[i] + yz[i] + yr[i]) + xr[i]^2 * yr[i] * (yl[i] * (yz[i] - 3 * yr[i]) + yz[i] * (yz[i] + yr[i])) + xz[i]^2 * yz[i] * (yl[i] * (-3 * yz[i] + yr[i]) + yr[i] * (yz[i] + yr[i])))))))) / (2 * sqrt((xl[i] * yl[i] * y_2_3[i] + xr[i] * y_1_2[i] * yr[i] + xz[i] * yz[i] * ((-yl[i]) + yr[i]))^2))
        # If y and w are 0, then 0/0=NaN
        if (is.nan(lambda_result)) {
            lambda_result <- 0
        }
        lambda <- c(lambda, lambda_result)

        # Calculate scaling factor A for each peak triplet
        A_result <- (-4 * w_1_2[i] * w_1_3[i] * w_2_3[i] * yl[i] * yz[i] * yr[i] * (xl[i] * yl[i] * y_2_3[i] + xr[i] * yr[i] * y_1_2[i] + xz[i] * yz[i] * (-y_1_3[i])) * lambda[i]) / (w_1_2[i]^4 * yl[i]^2 * yz[i]^2 - 2 * w_1_2[i]^2 * yl[i] * yz[i] * (w_1_3[i]^2 * yl[i] + w_2_3[i]^2 * yz[i]) * yr[i] + (w_1_3[i]^2 * yl[i] - w_2_3[i]^2 * yz[i])^2 * yr[i]^2)
        # If y and w are 0, then 0/0=NaN
        if (is.nan(A_result)) {
            A_result <- 0
        }
        A <- c(A, A_result)
    }

    # Calculate all initial lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(z), ncol = length(x))
    for (i in 1:length(z)) {
        # If A = 0, then the lorentz curve is a zero line
        if (A[i] == 0) {
            lorentz_curves_initial[i, ] <- 0
        } else {
            lorentz_curves_initial[i, ] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (x - w[i])^2)))
        }
    }
}