init_lorentz_curves_v12 <- function(x, y, pc, pl, pr) {
    # x <- spec$x
    # y <- spec$y_smooth
    # pc <- as.integer(spec$peak$center[spec$peak$high] - 1)
    # pl <- spec$peak$right[spec$peak$high] - 1
    # pr <- spec$peak$left[spec$peak$high] - 1
    msg("Initializing Lorentz curves")

    xl <- x[pl + 1]
    xc <- x[pc + 1]
    xr <- x[pr + 1]

    yl <- y[pl + 1]
    yc <- y[pc + 1]
    yr <- y[pr + 1]


    # Calculate mirrored points for asecnding/descending shoulders
    i <- which((yl < yc) & (yc < yr)) # ascending shoulders
    xr[i] <- 2 * xc[i] - xl[i]
    yr[i] <- yl[i]
    j <- which((yl > yc) & (yc > yr)) # descending shoulders
    xl[j] <- 2 * xc[j] - xr[j]
    yl[j] <- yr[j]

    

    # Calculate parameters w, lambda and A for the initial lorentz curves
    for (i in seq_along(pc)) {

        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if ((yl[i] < yc[i]) && (yc[i] < yr[i])) {
            xr[i] <- 2 * xc[i] - xl[i]
            yr[i] <- yl[i]
        }
        # For descending shoulders
        if ((yl[i] > yc[i]) && (yc[i] > yr[i])) {
            xl[i] <- 2 * xc[i] - xr[i]
            yl[i] <- yr[i]
        }
    }


    xl_2 <- c()
    xl_3 <- c()
    xc_3 <- c()

    yl_2 <- c()
    yl_3 <- c()
    yc_3 <- c()

    w_delta <- c()
    w <- c()
    lambda <- c()
    A <- c()

    # Calculate parameters w, lambda and A for the initial lorentz curves
    for (i in seq_along(pc)) {

        # Calculate mirrored points if necesccary
        # For ascending shoulders
        if ((yl[i] < yc[i]) && (yc[i] < yr[i])) {
            xr[i] <- 2 * xc[i] - xl[i]
            yr[i] <- yl[i]
        }
        # For descending shoulders
        if ((yl[i] > yc[i]) && (yc[i] > yr[i])) {
            xl[i] <- 2 * xc[i] - xr[i]
            yl[i] <- yr[i]
        }

        # Move triplet to zero position
        w_delta[i] <- xl[i]
        xl[i] <- xl[i] - w_delta[i]
        xc[i] <- xc[i] - w_delta[i]
        xr[i] <- xr[i] - w_delta[i]

        # Calculate difference of position of peak triplets
        xl_2 <- c(xl_2, xl[i] - xc[i])
        xl_3 <- c(xl_3, xl[i] - xr[i])
        xc_3 <- c(xc_3, xc[i] - xr[i])

        # Calculate difference of intensity values of peak triplets
        yl_2 <- c(yl_2, yl[i] - yc[i])
        yl_3 <- c(yl_3, yl[i] - yr[i])
        yc_3 <- c(yc_3, yc[i] - yr[i])

        # Calculate w for each peak triplet
        w_result <- (xl[i]^2 * yl[i] * yc_3[i] + xr[i]^2 * yr[i] * yl_2[i] + xc[i]^2 * yc[i] * (-yl_3[i])) / (2 * xl_2[i] * yl[i] * yc[i] - 2 * (xl_3[i] * yl[i] + (-xc_3[i]) * yc[i]) * yr[i])
        w_result <- w_result + w_delta[i]
        w <- c(w, w_result)
        # Wenn y Werte nach der H?henanpassung 0 werden, so ist w_new[i] NaN
        if (is.nan(w[i])) {
            w[i] <- 0
        }

        # Calculate lambda for each peak triplet
        lambda_result <- -((sqrt(abs((-xc[i]^4 * yc[i]^2 * yl_3[i]^2 - xl[i]^4 * yl[i]^2 * yc_3[i]^2 - xr[i]^4 * yl_2[i]^2 * yr[i]^2 + 4 * xc[i] * xr[i]^3 * yc[i] * ((-yl[i]) + yc[i]) * yr[i]^2 + 4 * xc[i]^3 * xr[i] * yc[i]^2 * yr[i] * ((-yl[i]) + yr[i]) + 4 * xl[i]^3 * yl[i]^2 * yc_3[i] * (xc[i] * yc[i] - xr[i] * yr[i]) + 4 * xl[i] * yl[i] * (xc[i]^3 * yc[i]^2 * yl_3[i] - xc[i] * xr[i]^2 * yc[i] * (yl[i] + yc[i] - 2 * yr[i]) * yr[i] + xr[i]^3 * yl_2[i] * yr[i]^2 - xc[i]^2 * xr[i] * yc[i] * yr[i] * (yl[i] - 2 * yc[i] + yr[i])) + 2 * xc[i]^2 * xr[i]^2 * yc[i] * yr[i] * (yl[i]^2 - 3 * yc[i] * yr[i] + yl[i] * (yc[i] + yr[i])) + 2 * xl[i]^2 * yl[i] * (-2 * xc[i] * xr[i] * yc[i] * yr[i] * (-2 * yl[i] + yc[i] + yr[i]) + xr[i]^2 * yr[i] * (yl[i] * (yc[i] - 3 * yr[i]) + yc[i] * (yc[i] + yr[i])) + xc[i]^2 * yc[i] * (yl[i] * (-3 * yc[i] + yr[i]) + yr[i] * (yc[i] + yr[i])))))))) / (2 * sqrt((xl[i] * yl[i] * yc_3[i] + xr[i] * yl_2[i] * yr[i] + xc[i] * yc[i] * ((-yl[i]) + yr[i]))^2))
        # If y and w are 0, then 0/0=NaN
        if (is.nan(lambda_result)) {
            lambda_result <- 0
        }
        lambda <- c(lambda, lambda_result)

        # Calculate scaling factor A for each peak triplet
        A_result <- (-4 * xl_2[i] * xl_3[i] * xc_3[i] * yl[i] * yc[i] * yr[i] * (xl[i] * yl[i] * yc_3[i] + xr[i] * yr[i] * yl_2[i] + xc[i] * yc[i] * (-yl_3[i])) * lambda[i]) / (xl_2[i]^4 * yl[i]^2 * yc[i]^2 - 2 * xl_2[i]^2 * yl[i] * yc[i] * (xl_3[i]^2 * yl[i] + xc_3[i]^2 * yc[i]) * yr[i] + (xl_3[i]^2 * yl[i] - xc_3[i]^2 * yc[i])^2 * yr[i]^2)
        # If y and w are 0, then 0/0=NaN
        if (is.nan(A_result)) {
            A_result <- 0
        }
        A <- c(A, A_result)
    }
    lc <- list(w = w, lambda = lambda, A = A, w_delta = w_delta)
    spec$lc$A <- A
    spec$lc$lambda <- lambda
    spec$lc$w <- w
    spec$lc$w_delta <- w_delta
    spec
}