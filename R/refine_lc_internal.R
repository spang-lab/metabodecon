refine_lc_internal_example <- function() {

    # Prepare inputs
    spectra <- read_spectra()
    spectra <- add_sfrs(spectra, sfr = c(11.44494, -1.8828), ask = FALSE, adjno = 1)
    spectra <- add_wsrs(spectra, wshw = 0.1527692, ask = FALSE, adjno = 1)
    spec <- spectra[[1]]
    spec <- rm_water_signal_v12(spec)
    spec <- rm_negative_signals_v12(spec)
    spec <- smooth_signals_v12(spec, reps = 2, k = 5)
    spec <- find_peaks_v12(spec)
    spec <- rm_peaks_with_low_scores_v12(spec, delta = 6.4)
    s <- spec; p <- spec$peak; h <- p$high; # styler: off
    x <- s$sdp; y <- s$y_smooth; pc <- p$center[h]; pl <- p$right[h]; pr <- p$left[h] # styler: off

    # Init lorentz curves
    specv12 <- init_lorentz_curves_v12(spec)
    lciv12 <- specv12$lc
    lci <- init_lorentz_curves_v13(x, y, pc, pl, pr)
    cat2("all.equal(lci, lcv12):", all.equal(lci, lciv12))

    # Refine lorentz curves
    p <- length(pc); n <- length(x) # styler: off
    A <- lci$A; lambda <- lci$lambda; w <- lci$w # styler: off
    Y <- matrix(nrow = n, ncol = p) # 131072 x 1227 for urine_1
    for (j in seq_along(pc)) {
        Y[, j] <- if (A[j] == 0) 0 else abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))
    }
    lcr1 <- refine_lc_internal(x, y, pc, pl, pr, Y)
    lcr2 <- refine_lc_internal(x, y, pc, pl, pr, Y)
    lcr3 <- refine_lc_internal(x, y, pc, pl, pr, Y)
    lcr <- list(lcr1, lcr2, lcr3)

    # Compare against old result
    lcrv12 <- refine_lorentz_curves_v12(specv12, nfit = 3)
    cat2("all.equal(lcr, lcrv12)", all.equal(lcr, lcrv12))

}

# x <- spec$sdp
# y <- spec$y_smooth
# pc <- spec$peak$center[spec$peak$high]
# pl <- spec$peak$right[spec$peak$high]
# pr <- spec$peak$left[spec$peak$high]
refine_lc_internal <- function(x, y, pc, pl, pr, Y) {

    # Old
    system.time({
        x_left <- x_left_1 <- x_left_2 <- c()
        x_center <- x_center_1 <- x_center_2 <- c()
        x_right <- x_right_1 <- x_right_2 <- c()
        y_left <- y_left_1 <- y_left_2 <- c()
        y_center <- y_center_1 <- y_center_2 <- c()
        y_right <- y_right_1 <- y_right_2 <- c()
        sum_left <- c()
        sum_peaks <- c()
        sum_right <- c()
        proportion_left <- c()
        proportion_peaks <- c()
        proportion_right <- c()
        ascending_shoulders <- c()
        descending_shoulders <- c()
        for (i in seq_along(pc)) {
            # Calculate the position of the peak triplets
            x_left[i] <- x[pl[i]]; x_center[i] <- x[pc[i]]; x_right[i] <- x[pr[i]]
            y_left[i] <- y[pl[i]]; y_center[i] <- y[pc[i]]; y_right[i] <- y[pr[i]]
            y_left_1[i] <- y_left[i]; y_center_1[i] <- y_center[i]; y_right_1[i] <- y_right[i] # backup
            x_left_1[i] <- x_left[i]; x_center_1[i] <- x_center[i]; x_right_1[i] <- x_right[i] # backup
            # Calculate the sum of all lorentz curves for each peak
            sum_left[i] <- sum(Y[pl[i], ])
            sum_peaks[i] <- sum(Y[pc[i], ])
            sum_right[i] <- sum(Y[pr[i], ])
            # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
            proportion_left[i] <- y[pl[i]] / sum_left[i]
            proportion_peaks[i] <- y[pc[i]] / sum_peaks[i]
            proportion_right[i] <- y[pr[i]] / sum_right[i]
            # Calculate the new heights of the peak triplets
            y_left[i] <- Y[pl[i], i] * proportion_left[i]
            y_right[i] <- Y[pr[i], i] * proportion_right[i]
            y_center[i] <- Y[pc[i], i] * proportion_peaks[i]
            y_left_2[i] <- y_left[i]; y_center_2[i] <- y_center[i]; y_right_2[i] <- y_right[i] # backup
            x_left_2[i] <- x_left[i]; x_center_2[i] <- x_center[i]; x_right_2[i] <- x_right[i] # backup
            # Calculate mirrored points if necesccary
            if ((y_left[i] < y_center[i]) && (y_center[i] < y_right[i])) { # For ascending shoulders
                x_right[i] <- 2 * x_center[i] - x_left[i]
                y_right[i] <- y_left[i]
                ascending_shoulders <- c(ascending_shoulders, i)
            }
            if ((y_left[i] > y_center[i]) && (y_center[i] > y_right[i])) { # For descending shoulders
                x_left[i] <- 2 * x_center[i] - x_right[i]
                y_left[i] <- y_right[i]
                descending_shoulders <- c(descending_shoulders, i)
            }
        }
    })

    # New
    # styler: off
    system.time({
        # Calculate the position and height of the peak triplets
        xl <- x[pl]; xc <- x[pc]; xr <- x[pr]
        yl <- y[pl]; yc <- y[pc]; yr <- y[pr]
        xl1 <- xl; xc1 <- xc; xr1 <- xr # backup
        yl1 <- yl; yc1 <- yc; yr1 <- yr # backup
        # Calculate the sum of all lorentz curves for each peak
        sl <- rowSums(Y[pl, ]); sc <- rowSums(Y[pc, ]); sr <- rowSums(Y[pr, ]) # sum_left, sum_peaks, sum_right
        # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
        ql <- yl / sl; qc <- yc / sc; qr <- yr / sr # proportion_left, proportion_peaks, proportion_right
        # Calculate new heights of peak triplets
        yl <- yl * ql; yc <- yc * qc; yr <- yr * qr
        xl2 <- xl; xc2 <- xc; xr2 <- xr # backup
        yl2 <- yl; yc2 <- yc; yr2 <- yr # backup
        # Calculate mirrored points if necesccary
        ia <- which((yl < yc) & (yc < yr)) # ascending_shoulders
        xr[ia] <- 2 * xc[ia] - xl[ia] # example: xl==7, xc==10, xr==12 ==> xr==2*10-7==13
        yr[ia] <- yl[ia]
        id <- which((yl > yc) & (yc > yr)) # descending_shoulders
        xl[id] <- 2 * xc[id] - xr[id] # example: xl==7, xc==10, xr==12 ==> xl==2*10-12==8
        yl[id] <- yr[id]
    })
    # styler: on

    # Check
    all.equal(xl1, x_left_1); all.equal(xc1, x_center_1); all.equal(xr, x_right_1)
    all.equal(yl1, y_left_1); all.equal(yc1, y_center_1); all.equal(yr, y_right_1)
    all.equal(sl, sum_left); all.equal(sc, sum_peaks); all.equal(sr, sum_right)
    all.equal(ql, proportion_left); all.equal(qc, proportion_peaks); all.equal(qr, proportion_right)
    all.equal(xl2, x_left_2); all.equal(xc2, x_center_2); all.equal(xr2, x_right_2)
    all.equal(yl2, y_left_2); all.equal(yc2, y_center_2); all.equal(yr2, y_right_2)
    all.equal(ia, ascending_shoulders)
    all.equal(id, descending_shoulders)
    all.equal(xl, x_left)
    all.equal(xr, x_right)
    all.equal(yl, y_left)

    for (i in seq_along(pc)) {
        # Move triplet to zero position
        wd[i] <- xl[i]
        xl[i] <- xl[i] - wd[i]
        xc[i] <- xc[i] - wd[i]
        xr[i] <- xr[i] - wd[i]

        # Calculate difference of peak triplet positions
        xlc <- c(xlc, xl[i] - xc[i])
        xlr <- c(xlr, xl[i] - xr[i])
        xcr <- c(xcr, xc[i] - xr[i])

        # Calculate difference of new intensity values of peak triplets
        Y <- c(Y, yl[i] - yc[i])
        ylr <- c(ylr, yl[i] - yr[i])
        ycr <- c(ycr, yc[i] - yr[i])

        # Calculate w for each peak triplet
        w_result <- (xl[i]^2 * yl[i] * ycr[i] + xr[i]^2 * yr[i] * Y[i] + xc[i]^2 * yc[i] * (-ylr[i])) / (2 * xlc[i] * yl[i] * yc[i] - 2 * (xlr[i] * yl[i] + (-xcr[i]) * yc[i]) * yr[i])
        w_result <- w_result + wd[i]
        w <- c(w, w_result)

        # If y values are getting 0 after height adjustment, then w[i]=NaN
        if (is.nan(w[i])) {
            w[i] <- 0
        }

        # Calculate lambda for each peak triplet
        lambda_result <- -((sqrt(abs(((-xc[i]^4 * yc[i]^2 * ylr[i]^2 - xl[i]^4 * yl[i]^2 * ycr[i]^2 - xr[i]^4 * Y[i]^2 * yr[i]^2 + 4 * xc[i] * xr[i]^3 * yc[i] * ((-yl[i]) + yc[i]) * yr[i]^2 + 4 * xc[i]^3 * xr[i] * yc[i]^2 * yr[i] * ((-yl[i]) + yr[i]) + 4 * xl[i]^3 * yl[i]^2 * ycr[i] * (xc[i] * yc[i] - xr[i] * yr[i]) + 4 * xl[i] * yl[i] * (xc[i]^3 * yc[i]^2 * ylr[i] - xc[i] * xr[i]^2 * yc[i] * (yl[i] + yc[i] - 2 * yr[i]) * yr[i] + xr[i]^3 * Y[i] * yr[i]^2 - xc[i]^2 * xr[i] * yc[i] * yr[i] * (yl[i] - 2 * yc[i] + yr[i])) + 2 * xc[i]^2 * xr[i]^2 * yc[i] * yr[i] * (yl[i]^2 - 3 * yc[i] * yr[i] + yl[i] * (yc[i] + yr[i])) + 2 * xl[i]^2 * yl[i] * (-2 * xc[i] * xr[i] * yc[i] * yr[i] * (-2 * yl[i] + yc[i] + yr[i]) + xr[i]^2 * yr[i] * (yl[i] * (yc[i] - 3 * yr[i]) + yc[i] * (yc[i] + yr[i])) + xc[i]^2 * yc[i] * (yl[i] * (-3 * yc[i] + yr[i]) + yr[i] * (yc[i] + yr[i]))))))))) / (2 * sqrt((xl[i] * yl[i] * ycr[i] + xr[i] * Y[i] * yr[i] + xc[i] * yc[i] * ((-yl[i]) + yr[i]))^2))

        # If y and w are 0, then 0/0=NaN
        if (is.nan(lambda_result)) {
            lambda_result <- 0
        }
        lambda <- c(lambda, lambda_result)

        # Calculate scaling factor A for each peak triplet
        A_result <- (-4 * xlc[i] * xlr[i] * xcr[i] * yl[i] * yc[i] * yr[i] * (xl[i] * yl[i] * ycr[i] + xr[i] * yr[i] * Y[i] + xc[i] * yc[i] * (-ylr[i])) * lambda[i]) / (xlc[i]^4 * yl[i]^2 * yc[i]^2 - 2 * xlc[i]^2 * yl[i] * yc[i] * (xlr[i]^2 * yl[i] + xcr[i]^2 * yc[i]) * yr[i] + (xlr[i]^2 * yl[i] - xcr[i]^2 * yc[i])^2 * yr[i]^2)

        # If y and w are 0, then 0/0=NaN
        if (is.nan(A_result)) {
            A_result <- 0
        }
        A <- c(A, A_result)

        # Calculate new lorentz curves
        # If y values are zero, then lorentz curves should also be zero
        for (j in seq_along(pc)) {
            Y[, j] <- if ((w[j] == 0) || (lambda[j] == 0) || (A[j] == 0)) 0 else abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))
        }
    }
}


refine_lc_internal_old <- function(spectrum_x, spectrum_y, filtered_left_position, filtered_peaks, filtered_right_position) {

    spectrum_x <- x
    spectrum_y <- y
    filtered_left_position <- pl - 1
    filtered_peaks <- pc - 1
    filtered_right_position <- pr - 1

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
    message(paste("\nNormed MSE value of iteration", b, "is: "))
    print(mse_normed)

    return(list(lorentz_curves_initial, spectrum_approx, mse_normed))
}