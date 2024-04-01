fit_lorentz_curves <- function(x, y, pc, pl, pr, nfit = 3) {

    msg("Initializing Lorentz curves")

    n <- length(x)
    p <- length(pc)
    wl <- x[pl]; wc <- x[pc]; wr <- x[pr] # position of peak triplets
    yl <- y[pl]; yc <- y[pc]; yr <- y[pr] # intensity of peak triplets

    # # Calculate mirrored points for ascending/descending shoulders
    # i <- which((yl < yc) & (yc < yr)) # ascending shoulders
    # j <- which((yl > yc) & (yc > yr)) # descending shoulders
    # xr[i] <- 2 * xc[i] - xl[i]
    # xl[j] <- 2 * xc[j] - xr[j]
    # yr[i] <- yl[i]
    # yl[j] <- yr[j]

    # # Move triplet to zero position
    # xd <- xl
    # xl <- xl - xd
    # xc <- xc - xd
    # xr <- xr - xd

    # # Calculate difference of position of peak triplets
    # xlc <- xl - xc
    # xlr <- xl - xr
    # xcr <- xc - xr

    # # Calculate difference of intensity values of peak triplets
    # ylc <- yl - yc
    # ylr <- yl - yr
    # ycr <- yc - yr

    # t1 <- xl^2 * yl * ycr
    # t2 <- xr^2 * yr * ylc
    # t3 <- xc^2 * yc * ylr
    # t4 <- 2 * xlc * yl * yc
    # t5 <- 2 * xcr * yc * yr
    # t6 <- 2 * xlr * yl * yr
    # w <- (t1 + t2 - t3) / (t4 + t5 - t6) + xd
    # w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.

    # lambda <- -((sqrt(abs((-xc^4 * yc^2 * ylr^2 - xl^4 * yl^2 * ycr^2 - xr^4 * ylc^2 * yr^2 + 4 * xc * xr^3 * yc * ((-yl) + yc) * yr^2 + 4 * xc^3 * xr * yc^2 * yr * ((-yl) + yr) + 4 * xl^3 * yl^2 * ycr * (xc * yc - xr * yr) + 4 * xl * yl * (xc^3 * yc^2 * ylr - xc * xr^2 * yc * (yl + yc - 2 * yr) * yr + xr^3 * ylc * yr^2 - xc^2 * xr * yc * yr * (yl - 2 * yc + yr)) + 2 * xc^2 * xr^2 * yc * yr * (yl^2 - 3 * yc * yr + yl * (yc + yr)) + 2 * xl^2 * yl * (-2 * xc * xr * yc * yr * (-2 * yl + yc + yr) + xr^2 * yr * (yl * (yc - 3 * yr) + yc * (yc + yr)) + xc^2 * yc * (yl * (-3 * yc + yr) + yr * (yc + yr)))))))) / (2 * sqrt((xl * yl * ycr + xr * ylc * yr + xc * yc * ((-yl) + yr))^2))
    # lambda[is.nan(lambda)] <- 0 

    # A <- (-4 * xlc * xlr * xcr * yl * yc * yr * (xl * yl * ycr + xr * yr * ylc + xc * yc * (-ylr)) * lambda) / (xlc^4 * yl^2 * yc^2 - 2 * xlc^2 * yl * yc * (xlr^2 * yl + xcr^2 * yc) * yr + (xlr^2 * yl - xcr^2 * yc)^2 * yr^2)
    # A[is.nan(A)] <- 0


    # Y <- matrix(nrow = n, ncol = p) # 131072 x 1227 for urine_1
    # for (j in seq_along(pc)) {
    #     Y[, j] <- if (A[j] == 0) 0 else Y[, j] <- abs(A[j] * (lambda[j] / (lambda[j]^2 + (x - w[j])^2)))
    # }

    # Approximation of lorentz curves
    for (b in 1:nfit) {
        cat3("Refining Lorentz Curves")

        # Debug variables
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
            sl[i] <- sum(Y[pl[i], ])
            sc[i] <- sum(Y[pc[i], ])
            sr[i] <- sum(Y[pr[i], ])

            # Calculate the proportion between original spectrum an the sum of the lorentz curves for each peak triplets position
            ql[i] <- y[pl[i]] / sl[i]
            qc[i] <- y[pc[i]] / sc[i]
            qr[i] <- y[pr[i]] / sr[i]
            yl1[i] <- y[pl[i]] # DEBUG
            yc1[i] <- y[pc[i]] # DEBUG
            yr1[i] <- y[pr[i]] # DEBUG

            # Calculate the new heights of the peak triplets
            yl[i] <- Y[pl[i], i] * ql[i]
            yc[i] <- Y[pc[i], i] * qc[i]
            yr[i] <- Y[pr[i], i] * qr[i]
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
            Y[, i] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (x - w[i])^2)))
        }

        y_normed <- y / sum(y)
        yhat <- rowSums(Y)
        yhat_normed <- yhat / sum(yhat)
        se_normed <- (y_normed - yhat_normed)^2
        mse_normed <- (1 / length(se_normed)) * sum(se_normed)
        cat3(sprintf("Normed MSE: %.22f", mse_normed))
    }

    # Calculate the integrals for each lorentz curve
    integrals <- matrix(nrow = 1, ncol = length(l))
    for (i in seq_along(l)) {
        integrals[1, i] <- A_new[i] * (atan((-w[i] + (spec$n / spec$sfx)) / l[i]) - atan((-w[i]) / l[i]))
    }

    spec$lcr <- list(
        A_new = A_new,
        l = l,
        w = w,
        spectrum_approx = spectrum_approx,
        y_normed = y_normed,
        spectrum_approx_normed = spectrum_approx_normed,
        difference_normed = difference_normed,
        mse_normed = mse_normed,
        integrals = integrals
    )
    spec
}
