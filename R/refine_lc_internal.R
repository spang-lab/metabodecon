refine_lc_internal_example <- function() {

    # Prepare Inputs
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

    # Init Lorentz Curves
    specv12 <- init_lorentz_curves_v12(spec)
    lciv12 <- specv12$lc
    lci <- init_lorentz_curves_v13(x, y, pc, pl, pr)
    cat2("all.equal(lci, lcv12):", all.equal(lci, lciv12))

    # Refine Lorentz Curves
    p <- length(pc); n <- length(x) # styler: off
    A <- lci$A; lambda <- lci$lambda; w <- lci$w # styler: off
    Y <- matrix(nrow = n, ncol = p) # 131072 x 1227 for urine_1
    for (j in seq_along(pc)) {
        Y[, j] <- if (A[j] == 0) 0 else abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))
    }
    Yt <- t(Y)

    # Refine Lorentz Curves Internal
    system.time(old2 <- refine_lc_internal_old(x, y, pl, pc, pr, Yt))
    system.time(new <- refine_lc_internal(x, y, pc, pl, pr, Y))

    # Check correctness
    all.equal(new$xl1, old$w_1_bak1);     all.equal(new$xc1, old$w_2_bak1);         all.equal(new$xr1, old$w_3_bak1) # x1 == first backup x
    all.equal(new$yl1, old$y_1_bak1);     all.equal(new$yc1, old$y_2_bak1);         all.equal(new$yr1, old$y_3_bak1) # y1 == first backup y

    all.equal(new$sl,  old$sum_left);     all.equal(new$sc,  old$sum_peaks);        all.equal(new$sr,  old$sum_right) # s == sum of lorentz curves aeptp
    all.equal(new$ql,  old$proportion_1); all.equal(new$qc,  old$proportion_peaks); all.equal(new$qr,  old$proportion_3) # q == proportion y/s
    all.equal(new$yl2, y_1_2);            all.equal(new$yc2, old$y_2_2);            all.equal(new$yr2, old$y_3_2) # y2 == second backup y
    all.equal(new$ia, ascending_shoulders) # ia == index of ascending shoulders
    all.equal(new$id, descending_shoulders) # id == index of descending shoulders
    all.equal(new$xl2, x_1_2); all.equal(new$xc2, x_2_2); all.equal(new$xr2, x_3_2) # x2 == second backup x
    all.equal(new$yl3, y_1_3); all.equal(new$yc3, y_2_3); all.equal(new$yr3, y_3_3) # y3 == third backup y
    all.equal(new$xd, w_delta)
    all.equal(new$xl2, x_1_2); all.equal(new$xc2, x_2_2); all.equal(new$xr2, x_3_2) # x3 == third backup x

}

# x <- spec$sdp
# y <- spec$y_smooth
# pc <- spec$peak$center[spec$peak$high]
# pl <- spec$peak$right[spec$peak$high]
# pr <- spec$peak$left[spec$peak$high]
refine_lc_internal <- function(x, y, pc, pl, pr, Y) {

    xl <- x[pl]; xc <- x[pc]; xr <- x[pr] # x == coords at each peak triplet position (aeptp)
    yl <- y[pl]; yc <- y[pc]; yr <- y[pr] # y == height aeptp
    xl1 <- xl; xc1 <- xc; xr1 <- xr # x1 == first backup x
    yl1 <- yl; yc1 <- yc; yr1 <- yr # y1 == first backup y

    sl <- rowSums(Y[pl, ]); sc <- rowSums(Y[pc, ]); sr <- rowSums(Y[pr, ]) # s == sum of lorentz curves aeptp
    ql <- yl / sl; qc <- yc / sc; qr <- yr / sr # q == proportion y/s
    yl <- diag(Y[pl, ]) * ql; yc <- diag(Y[pc, ]) * qc; yr <- diag(Y[pr, ]) * qr # y == new height aeptp
    yl2 <- yl; yc2 <- yc; yr2 <- yr # y2 == second backup y

    ia <- which((yl < yc) & (yc < yr)) # ia == index of ascending shoulders
    xr[ia] <- 2 * xc[ia] - xl[ia] # x == x with ascending shoulders replaced by mirrored points (1)
    yr[ia] <- yl[ia] # same for y
    id <- which((yl > yc) & (yc > yr)) # id == index of descending shoulders
    xl[id] <- 2 * xc[id] - xr[id] # x == x with descending shoulders replaced by mirrored points (2)
    yl[id] <- yr[id] # same for y
    xl2 <- xl; xc2 <- xc; xr2 <- xr # x2 == second backup x
    yl3 <- yl; yc3 <- yc; yr3 <- yr # y3 == third backup y
    # (1) Example: xl == 7, xc == 10, xr == 12 ==> xr == 2 * 10 - 7  == 13
    # (2) Example: xl == 7, xc == 10, xr == 12 ==> xl == 2 * 10 - 12 == 8

    xd <- xl; xl <- xl - xd; xc <- xc - xd; xr <- xr - xd # x == distance to left border
    xl3 <- xl; xc3 <- xc; xr3 <- xr # x3 == third backup x

    return(as.list(environment()))
}


refine_lc_internal_old <- function(x, y, pl, pc, pr, Yt) {

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

        # Calculate the proprotion between original spectrum an the sum of the lorentz curves for each peak triplets position
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
    message(paste("\nNormed MSE value is: "))
    print(mse_normed)

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