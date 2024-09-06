# Private Main #####

fit_lorentz_curves <- function(spec, nfit = 3) {
    logf("Initializing Lorentz curves")
    spec$lci <- lc <- init_lc_v14(spec)
    spec$lca <- vector("list", length = nfit)
    logf("Refining Lorentz Curves")
    for (i in 1:nfit) {
        spec$lca[[i]] <- lc <- refine_lc_v14(spec, lc$Z)
    }
    A <- lc$A; lambda <- lc$lambda; w <- lc$w
    integrals <- A * (atan((-w + (spec$n / spec$sf[1])) / lambda) - atan((-w) / lambda))
    spec$lcr <- list(A = A, lambda = lambda, w = w, integrals = integrals)
    spec
}

# Private Helpers #####

#' @title Initialize Lorentz Curve parameters from a spectrum with selected peak triplets
#' @param spec List with elements: `x`, `y`, `peak` where `peak` is a list with elements `center`, `left`, `right` and `high`.
#' @noRd
init_lc_v14 <- function(spec, verbose = TRUE) {
    # Init values
    p <- spec$peak
    ir <- p$right[p$high]; ic <- p$center[p$high]; il <- p$left[p$high] # index of each peak triplet position (PTP)
    lmr <- sort(unique(c(il, ic, ir))) # combined PTP indices
    rr <- match(ir, lmr);  rc <- match(ic, lmr);   rl <- match(il, lmr) # rank  of each PTP
    x <- spec$sdp; y <- spec$y_smooth # x and y value for each data point
    xlmr <- x[lmr]; # x value for each PTP
    yr <- y[ir]; yc <- y[ic]; yl <- y[il]; # intensity of each PTP
    xr <- x[ir]; xc <- x[ic]; xl <- x[il]; # position of each PTP

    # Replace shoulders
    as <- (yr > yc) & (yc > yl) # ascending shoulders (AS)
    ds <- (yr < yc) & (yc < yl) # descending shoulders (DS)
    xl[ds] <- 2 * xc[ds] - xr[ds] # replace DS with mirrored points (MP)
    xr[as] <- 2 * xc[as] - xl[as] # replace AS with MP
    yl[ds] <- yr[ds] # replace DS with MP
    yr[as] <- yl[as] # replace AS with MP

    # Calculate distances
    wr  <- xr - xr; wc  <- xc - xr; wl  <- xl - xr # express positions wr/wc/wl as "distance to right border"
    wrc <- wr - wc; wrl <- wr - wl; wcl <- wc - wl # x-distance between PTPs
    yrc <- yr - yc; yrl <- yr - yl; ycl <- yc - yl # y-distance between PTPs

    # Estimate parameters
    w <- calc_w_v14(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr)
    lambda <- calc_lambda_v14(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr)
    A <- calc_A_v14(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, lambda, w, xr)

    # Calculate contribution of each lorentz curve to each PTP data point
    Z <- matrix(nrow = length(lmr), ncol = length(ic)) # 3614 x 1227 for urine_1
    for (j in seq_along(ic)) Z[, j] <- if (A[j] == 0) 0 else abs(A[j] * (lambda[j] / (lambda[j]^2 + (xlmr - w[j])^2)))

    # Print MSE
    mse <- mean((y[lmr] - rowSums(Z))^2)
    if (verbose) logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, as, ds)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

refine_lc_v14 <- function(spec, Z) {

    # Init x and y values
    x <- spec$sdp; y <- spec$y_smooth # x and y value for each data point

    # Init peak related variables
    p <- spec$peak
    ir <- p$right[p$high]; ic <- p$center[p$high]; il <- p$left[p$high] # index of each peak triplet position (PTP)
    lmr <- sort(unique(c(il, ic, ir))) # combined PTP indices
    rr <- match(ir, lmr);  rc <- match(ic, lmr);   rl <- match(il, lmr) # rank  of each PTP
    xlmr <- x[lmr]; # x value for each PTP
    yr <- y[ir]; yc <- y[ic]; yl <- y[il]; # intensity of each PTP
    xr <- x[ir]; xc <- x[ic]; xl <- x[il]; # position of each PTP
    sl <- sc <- sr <- numeric(length(ic));
    ql <- qc <- qr <- numeric(length(ic));

    # Init distance related variables
    wr  <- wc  <- wl  <- numeric(length(ic));
    wrc <- wrl <- wcl <- numeric(length(ic));
    yrc <- yrl <- ycl <- numeric(length(ic));

    # Init lorentz curves parameters and matrices
    A <- lambda <- w <- numeric(length(ic))

    for (i in seq_along(il)) {

        # Calculate the sum of all lorentz curves for each data point
        sl[i] <- sum(Z[rl[i], ])
        sc[i] <- sum(Z[rc[i], ])
        sr[i] <- sum(Z[rr[i], ])

        # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
        ql[i] <- yl[i] / sl[i]
        qc[i] <- yc[i] / sc[i]
        qr[i] <- yr[i] / sr[i]

        # Calculate the new heights of the peak triplets
        yl[i] <- Z[rl[i], i] * ql[i]
        yc[i] <- Z[rc[i], i] * qc[i]
        yr[i] <- Z[rr[i], i] * qr[i]

        # Calculate mirrored points if necesccary
        if ((yl[i] < yc[i]) && (yc[i] < yr[i])) { # Ascending shoulder
            xr[i] <- 2 * xc[i] - xl[i]
            yr[i] <- yl[i]
        }
        if ((yl[i] > yc[i]) && (yc[i] > yr[i])) { # Descending shoulder
            xl[i] <- 2 * xc[i] - xr[i]
            yl[i] <- yr[i]
        }

        # Calculate distances between peak triplet positions and intensities
        wr[i] <- xr[i] - xr[i]
        wc[i] <- xc[i] - xr[i]
        wl[i] <- xl[i] - xr[i]
        wrc[i] <- wr[i] - wc[i]
        wrl[i] <- wr[i] - wl[i]
        wcl[i] <- wc[i] - wl[i]
        yrc[i] <- yr[i] - yc[i]
        yrl[i] <- yr[i] - yl[i]
        ycl[i] <- yc[i] - yl[i]

        # Estimate parameters
        w[i] <- calc_w_v14(wr[i], wc[i], wl[i], yr[i], yc[i], yl[i], wrc[i], wrl[i], wcl[i], yrc[i], yrl[i], ycl[i], xr[i])
        lambda[i] <- calc_lambda_v14(wr[i], wc[i], wl[i], yr[i], yc[i], yl[i], wrc[i], wrl[i], wcl[i], yrc[i], yrl[i], ycl[i], xr[i])
        A[i] <- calc_A_v14(wr[i], wc[i], wl[i], yr[i], yc[i], yl[i], wrc[i], wrl[i], wcl[i], yrc[i], yrl[i], ycl[i], lambda[i], w[i], xr[i])

        # Calculate contribution of each lorentz curve to each data point
        cond <- (w[i] == 0) || (lambda[i] == 0) || (A[i] == 0)
        Z[, i] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (xlmr - w[i])^2)))
    }

    # Print MSE
    mse <- mean((y[lmr] - rowSums(Z))^2)
    logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, sl, sc, sr, ql, qc, qr)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

calc_w_v14 <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    t1 <- wr^2 * yr * ycl
    t2 <- wl^2 * yl * yrc
    t3 <- wc^2 * yc * yrl
    t4 <- 2 * wrc * yr * yc
    t5 <- 2 * wcl * yc * yl
    t6 <- 2 * wrl * yr * yl
    w <- (t1 + t2 - t3) / (t4 + t5 - t6) + xr
    w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.
    w
}

calc_lambda_v14 <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    lambda <- -((sqrt(abs((-wc^4 * yc^2 * yrl^2 - wr^4 * yr^2 * ycl^2 - wl^4 * yrc^2 * yl^2 + 4 * wc * wl^3 * yc * ((-yr) + yc) * yl^2 + 4 * wc^3 * wl * yc^2 * yl * ((-yr) + yl) + 4 * wr^3 * yr^2 * ycl * (wc * yc - wl * yl) + 4 * wr * yr * (wc^3 * yc^2 * yrl - wc * wl^2 * yc * (yr + yc - 2 * yl) * yl + wl^3 * yrc * yl^2 - wc^2 * wl * yc * yl * (yr - 2 * yc + yl)) + 2 * wc^2 * wl^2 * yc * yl * (yr^2 - 3 * yc * yl + yr * (yc + yl)) + 2 * wr^2 * yr * (-2 * wc * wl * yc * yl * (-2 * yr + yc + yl) + wl^2 * yl * (yr * (yc - 3 * yl) + yc * (yc + yl)) + wc^2 * yc * (yr * (-3 * yc + yl) + yl * (yc + yl)))))))) / (2 * sqrt((wr * yr * ycl + wl * yrc * yl + wc * yc * ((-yr) + yl))^2))
    lambda[is.nan(lambda)] <- 0
    lambda
}

calc_A_v14 <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, lambda, w, xr) {
    A <- (-4 * wrc * wrl * wcl * yr * yc * yl * (wr * yr * ycl + wl * yl * yrc + wc * yc * (-yrl)) * lambda) / (wrc^4 * yr^2 * yc^2 - 2 * wrc^2 * yr * yc * (wrl^2 * yr + wcl^2 * yc) * yl + (wrl^2 * yr - wcl^2 * yc)^2 * yl^2)
    A[is.nan(A)] <- 0
    A
}

# Deprecated #####

init_lorentz_curves_v10 <- function(spectrum_x, spectrum_y, filtered_peaks, filtered_left_position, filtered_right_position, save_scores) {

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
        if ((y_1[i] < y_2[i]) && (y_2[i] < y_3[i])) {
            w_3[i] <- 2 * w_2[i] - w_1[i]
            y_3[i] <- y_1[i]
        }
        # For descending shoulders
        if ((y_1[i] > y_2[i]) && (y_2[i] > y_3[i])) {
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
    list(A = A, lambda = lambda, w = w, w_delta = w_delta)
}

init_lorentz_curves_v12 <- function(spec) {
    logf("Initializing Lorentz curves")
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
        if (is.na(((y_1[i] < y_2[i]) & (y_2[i] < y_3[i])))) {
            1
        }
        if ((y_1[i] < y_2[i]) && (y_2[i] < y_3[i])) {
            w_3[i] <- 2 * w_2[i] - w_1[i]
            y_3[i] <- y_1[i]
        }
        # For descending shoulders
        if ((y_1[i] > y_2[i]) && (y_2[i] > y_3[i])) {
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

# x <- spec$sdp
# y <- spec$y_smooth
# pc <- spec$peak$center[spec$peak$high]
# pl <- spec$peak$right[spec$peak$high]
# pr <- spec$peak$left[spec$peak$high]
init_lorentz_curves_v13 <- function(x, y, pc, pl, pr) {

    logf("Initializing Lorentz curves")

    xl <- x[pl] # position of peak triplets
    xc <- x[pc]
    xr <- x[pr]

    yl <- y[pl] # intensity of peak triplets
    yc <- y[pc]
    yr <- y[pr]

    # Calculate mirrored points for ascending/descending shoulders
    i <- which((yl < yc) & (yc < yr)) # ascending shoulders
    j <- which((yl > yc) & (yc > yr)) # descending shoulders
    xr[i] <- 2 * xc[i] - xl[i]
    xl[j] <- 2 * xc[j] - xr[j]
    yr[i] <- yl[i]
    yl[j] <- yr[j]

    # Move triplet to zero position
    xd <- xl
    xl <- xl - xd
    xc <- xc - xd
    xr <- xr - xd

    # Calculate difference of position of peak triplets
    xlc <- xl - xc
    xlr <- xl - xr
    xcr <- xc - xr

    # Calculate difference of intensity values of peak triplets
    ylc <- yl - yc
    ylr <- yl - yr
    ycr <- yc - yr

    t1 <- xl^2 * yl * ycr
    t2 <- xr^2 * yr * ylc
    t3 <- xc^2 * yc * ylr
    t4 <- 2 * xlc * yl * yc
    t5 <- 2 * xcr * yc * yr
    t6 <- 2 * xlr * yl * yr
    w <- (t1 + t2 - t3) / (t4 + t5 - t6) + xd
    w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.

    lambda <- -((sqrt(abs((-xc^4 * yc^2 * ylr^2 - xl^4 * yl^2 * ycr^2 - xr^4 * ylc^2 * yr^2 + 4 * xc * xr^3 * yc * ((-yl) + yc) * yr^2 + 4 * xc^3 * xr * yc^2 * yr * ((-yl) + yr) + 4 * xl^3 * yl^2 * ycr * (xc * yc - xr * yr) + 4 * xl * yl * (xc^3 * yc^2 * ylr - xc * xr^2 * yc * (yl + yc - 2 * yr) * yr + xr^3 * ylc * yr^2 - xc^2 * xr * yc * yr * (yl - 2 * yc + yr)) + 2 * xc^2 * xr^2 * yc * yr * (yl^2 - 3 * yc * yr + yl * (yc + yr)) + 2 * xl^2 * yl * (-2 * xc * xr * yc * yr * (-2 * yl + yc + yr) + xr^2 * yr * (yl * (yc - 3 * yr) + yc * (yc + yr)) + xc^2 * yc * (yl * (-3 * yc + yr) + yr * (yc + yr)))))))) / (2 * sqrt((xl * yl * ycr + xr * ylc * yr + xc * yc * ((-yl) + yr))^2))
    lambda[is.nan(lambda)] <- 0

    A <- (-4 * xlc * xlr * xcr * yl * yc * yr * (xl * yl * ycr + xr * yr * ylc + xc * yc * (-ylr)) * lambda) / (xlc^4 * yl^2 * yc^2 - 2 * xlc^2 * yl * yc * (xlr^2 * yl + xcr^2 * yc) * yr + (xlr^2 * yl - xcr^2 * yc)^2 * yr^2)
    A[is.nan(A)] <- 0

    lc <- list(A = A, lambda = lambda, w = w, w_delta = xd)
    lc
}

#' @noRd
refine_lorentz_curves_v12 <- function(spec, nfit) {
    logf("Refining Lorentz curves")

    spectrum_x <- spec$sdp
    spectrum_y <- spec$y_smooth
    filtered_peaks <- as.integer(spec$peak$center[spec$peak$high] - 1)
    filtered_left_position <- spec$peak$right[spec$peak$high] - 1
    filtered_right_position <- spec$peak$left[spec$peak$high] - 1
    A <- spec$lc$A
    lambda <- spec$lc$lambda
    w <- spec$lc$w
    nfit <- nfit

    # Calculate all initial lorentz curves
    lorentz_curves_initial <- matrix(nrow = length(filtered_peaks), ncol = length(spectrum_x))
    for (i in seq_along(filtered_peaks)) {
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
            if ((y_1_new[i] < y_2_new[i]) && (y_2_new[i] < y_3_new[i])) {
                w_3_new[i] <- 2 * w_2_new[i] - w_1_new[i]
                y_3_new[i] <- y_1_new[i]
            }
            # For descending shoulders
            if ((y_1_new[i] > y_2_new[i]) && (y_2_new[i] > y_3_new[i])) {
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
        logf("Normed MSE after iteration %d: %.22f", b, mse_normed)
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(lambda_new))
    for (i in seq_along(lambda_new)) {
        integrals[1, i] <- A_new[i] * (atan((-w_new[i] + (spec$n / spec$sf[1])) / lambda_new[i]) - atan((-w_new[i]) / lambda_new[i]))
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
        integrals = integrals
    )
    spec
}
