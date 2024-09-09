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
