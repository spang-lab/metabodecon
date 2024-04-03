# x1 x1 x3 y1 y2 y3 i1 i2 i3 j1 j2 j3 s

refine_lc_internal_v11 <- function(x, y, pl, pc, pr, Yt) {

    cat3("Refining Lorentz Curves")

    # Debug variables
    wl1 <- c(); wc1 <- c(); wr1 <- c() # DEBUG
    wl2 <- c(); wc2 <- c(); wr2 <- c() # DEBUG
    wl3 <- c(); wc3 <- c(); wr3 <- c() # DEBUG
    yl1 <- c(); yc1 <- c(); yr1 <- c() # DEBUG
    yl2 <- c(); yc2 <- c(); yr2 <- c() # DEBUG
    yl3 <- c(); yc3 <- c(); yr3 <- c() # DEBUG
    yl4 <- c(); yc4 <- c(); yr4 <- c() # DEBUG

    # Calculate the position of the peak triplets
    wl <- x[pl]; wc <- x[pc]; wr <- x[pr]
    wl1 <- wl; wc1 <- wc; wr1 <- wr

    # Calculate new heights of peak triplets
    yl <- c(); yc <- c(); yr <- c()
    wlc <- c(); wlr <- c(); wcr <- c()
    ylc <- c(); ylr <- c(); ycr <- c()
    wd <- c(); w <- c(); lambda <- c(); A <- c()
    sl <- c(); sc <- c(); sr <- c()
    ql <- c(); qc <- c(); qr <- c()

    for (i in seq_along(pc)) {

        # Calculate the sum of all lorentz curves for each data point
        sl[i] <- sum(Yt[, pl[i]])
        sc[i] <- sum(Yt[, pc[i]])
        sr[i] <- sum(Yt[, pr[i]])

        # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
        ql[i] <- y[pl[i]] / sl[i]
        qc[i] <- y[pc[i]] / sc[i]
        qr[i] <- y[pr[i]] / sr[i]
        yl1[i] <- y[pl[i]] # DEBUG
        yc1[i] <- y[pc[i]] # DEBUG
        yr1[i] <- y[pr[i]] # DEBUG

        # Calculate the new heights of the peak triplets
        yl[i] <- Yt[i, pl[i]] * ql[i]
        yc[i] <- Yt[i, pc[i]] * qc[i]
        yr[i] <- Yt[i, pr[i]] * qr[i]
        yl2[i] <- yl[i] # DEBUG
        yc2[i] <- yc[i] # DEBUG
        yr2[i] <- yr[i] # DEBUG

        # Calculate mirrored points if necesccary
        if ((yl[i] < yc[i]) && (yc[i] < yr[i])) { # For ascending shoulders
            wr[i] <- 2 * wc[i] - wl[i]
            yr[i] <- yl[i]
        }
        if ((yl[i] > yc[i]) && (yc[i] > yr[i])) { # For descending shoulders
            wl[i] <- 2 * wc[i] - wr[i]
            yl[i] <- yr[i]
        }
        wl2[i] <- wl[i] # DEBUG
        wc2[i] <- wc[i] # DEBUG
        wr2[i] <- wr[i] # DEBUG
        yl3[i] <- yl[i] # DEBUG
        yc3[i] <- yc[i] # DEBUG
        yr3[i] <- yr[i] # DEBUG

        # Move triplet to zero position
        wd[i] <- wl[i]
        wl[i] <- wl[i] - wd[i]
        wc[i] <- wc[i] - wd[i]
        wr[i] <- wr[i] - wd[i]
        wl3[i] <- wl[i] # DEBUG
        wc3[i] <- wc[i] # DEBUG
        wr3[i] <- wr[i] # DEBUG

        # Calculate difference of peak triplet positions
        wlc <- c(wlc, wl[i] - wc[i])
        wlr <- c(wlr, wl[i] - wr[i])
        wcr <- c(wcr, wc[i] - wr[i])

        # Calculate difference of new intensity values of peak triplets
        ylc <- c(ylc, yl[i] - yc[i])
        ylr <- c(ylr, yl[i] - yr[i])
        ycr <- c(ycr, yc[i] - yr[i])
        yl4[i] <- yl[i] # DEBUG
        yc4[i] <- yc[i] # DEBUG
        yr4[i] <- yr[i] # DEBUG

        # Calculate new position w for each peak triplet
        w_new <- (wl[i]^2 * yl[i] * ycr[i] + wr[i]^2 * yr[i] * ylc[i] + wc[i]^2 * yc[i] * (-ylr[i])) / (2 * wlc[i] * yl[i] * yc[i] - 2 * (wlr[i] * yl[i] + (-wcr[i]) * yc[i]) * yr[i]) + wd[i]
        w[i] <- if (is.nan(w_new)) 0 else w_new # If y values are getting 0 after height adjustment, then w[i]=NaN

        # Calculate lambda for each peak triplet
        lambda_new <- -((sqrt(abs(((-wc[i]^4 * yc[i]^2 * ylr[i]^2 - wl[i]^4 * yl[i]^2 * ycr[i]^2 - wr[i]^4 * ylc[i]^2 * yr[i]^2 + 4 * wc[i] * wr[i]^3 * yc[i] * ((-yl[i]) + yc[i]) * yr[i]^2 + 4 * wc[i]^3 * wr[i] * yc[i]^2 * yr[i] * ((-yl[i]) + yr[i]) + 4 * wl[i]^3 * yl[i]^2 * ycr[i] * (wc[i] * yc[i] - wr[i] * yr[i]) + 4 * wl[i] * yl[i] * (wc[i]^3 * yc[i]^2 * ylr[i] - wc[i] * wr[i]^2 * yc[i] * (yl[i] + yc[i] - 2 * yr[i]) * yr[i] + wr[i]^3 * ylc[i] * yr[i]^2 - wc[i]^2 * wr[i] * yc[i] * yr[i] * (yl[i] - 2 * yc[i] + yr[i])) + 2 * wc[i]^2 * wr[i]^2 * yc[i] * yr[i] * (yl[i]^2 - 3 * yc[i] * yr[i] + yl[i] * (yc[i] + yr[i])) + 2 * wl[i]^2 * yl[i] * (-2 * wc[i] * wr[i] * yc[i] * yr[i] * (-2 * yl[i] + yc[i] + yr[i]) + wr[i]^2 * yr[i] * (yl[i] * (yc[i] - 3 * yr[i]) + yc[i] * (yc[i] + yr[i])) + wc[i]^2 * yc[i] * (yl[i] * (-3 * yc[i] + yr[i]) + yr[i] * (yc[i] + yr[i]))))))))) / (2 * sqrt((wl[i] * yl[i] * ycr[i] + wr[i] * ylc[i] * yr[i] + wc[i] * yc[i] * ((-yl[i]) + yr[i]))^2))
        lambda[i] <- if (is.nan(lambda_new)) 0 else lambda_new # If y and w are 0, then 0/0=NaN

        # Calculate scaling factor A for each peak triplet
        A_new <- (-4 * wlc[i] * wlr[i] * wcr[i] * yl[i] * yc[i] * yr[i] * (wl[i] * yl[i] * ycr[i] + wr[i] * yr[i] * ylc[i] + wc[i] * yc[i] * (-ylr[i])) * lambda[i]) / (wlc[i]^4 * yl[i]^2 * yc[i]^2 - 2 * wlc[i]^2 * yl[i] * yc[i] * (wlr[i]^2 * yl[i] + wcr[i]^2 * yc[i]) * yr[i] + (wlr[i]^2 * yl[i] - wcr[i]^2 * yc[i])^2 * yr[i]^2)
        A[i] <- if (is.nan(A_new)) 0 else A_new # If y and w are 0, then 0/0=NaN

        cond <- (w[i] == 0) || (lambda[i] == 0) || (A[i] == 0)
        Yt[i, ] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (x - w[i])^2)))
    }

    y_normed <- y / sum(y)
    yhat <- colSums(Yt)
    yhat_normed <- yhat / sum(yhat)
    se_normed <- (y_normed - yhat_normed)^2
    mse_normed <- (1 / length(se_normed)) * sum(se_normed)
    cat3(sprintf("Normed MSE: %.22f", mse_normed))

    return(list(
        wl1 = wl1, wc1 = wc1, wr1 = wr1, # DEBUG
        wl2 = wl2, wc2 = wc2, wr2 = wr2, # DEBUG
        wl3 = wl3, wc3 = wc3, wr3 = wr3, # DEBUG
        yl1 = yl1, yc1 = yc1, yr1 = yr1, # DEBUG
        yl2 = yl2, yc2 = yc2, yr2 = yr2, # DEBUG
        yl3 = yl3, yc3 = yc3, yr3 = yr3, # DEBUG
        yl4 = yl4, yc4 = yc4, yr4 = yr4, # DEBUG
        wl = wl, wc = wc, wr = wr,
        yl = yl, yc = yc, yr = yr,
        wlc = wlc, wlr = wlr, wcr = wcr,
        ylc = ylc, ylr = ylr, ycr = ycr,
        wd = wd,
        w = w,
        lambda = lambda,
        A = A,
        sl = sl, sc = sc, sr = sr,
        ql = ql, qc = qc, qr = qr,
        Yt = Yt,
        spectrum_approx = yhat,
        mse_normed = mse_normed
    ))
}