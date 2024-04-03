refine_lc_internal_v14 <- function(spec, Y, Z) {

    cat3("Refining Lorentz Curves")

    # Init x and y values
    x <- spec$sdp; y <- spec$y_smooth # x and y value for each data point

    # Init peak related variables
    # TODO: include this info directly inside spec$peak, so we dont have to calculate over and over again
    p <- spec$peak
    ir <- p$right[p$high]; ic <- p$center[p$high]; il <- p$left[p$high] # index of each peak triplet position (PTP)
    rr <- match(ir, lmr);  rc <- match(ic, lmr);   rl <- match(il, lmr) # rank  of each PTP
    lmr <- sort(unique(c(il, ic, ir))) # combined PTP indices
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
        # Y[, i] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (x - w[i])^2)))
        Z[, i] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (xlmr - w[i])^2)))
    }

    # y_normed <- y / sum(y)
    # yhat <- rowSums(Y)
    # yhat_normed <- yhat / sum(yhat)
    # se_normed <- (y_normed - yhat_normed)^2
    # mse_normed <- mean(se_normed)
    # cat3(sprintf("Normed MSE: %.22f", mse_normed))

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, sl, sc, sr, ql, qc, qr)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    cat3("Done")
    named_list(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}
