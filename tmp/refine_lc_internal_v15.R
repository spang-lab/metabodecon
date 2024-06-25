# x1 x1 x3 y1 y2 y3 i1 i2 i3 j1 j2 j3 s

refine_lc_internal_v15 <- function(spec, lc) {

    # Phases:
    # 1. Increase/decrease y[i] at PTP[i] if SLC[i] is smaller/larger than y[i]
    # 2. Replace shoulders with mirrored points
    # 3. Move PTPs to zero position
    # 4. Estimate new LC parameters based on y

    cat3("Refining Lorentz Curves")
    # Info from last run
    x <- spec$sdp # Postion of data points
    y <- spec$y_smooth # Intensity of data points
    n <- length(x) # Number of data points
    m <- length(lc$uid) # Number of PT data points
    p <- nrow(lc$id) # Number of peaks == Number of Lorenz curves
    id <- lc$id # PT indices (p x 3)
    rnk <- lc$rnk # PT ranks (p x 3)

    w0 <- lc$w0  # PT positions unchanged (p x 3)
    y0 <- lc$y0  # PT intensities unchanged (p x 3)

    Y <- lc$Y # Intensity of each lorenz curve at each data points (n x p)
    Z <- lc$Z # Intensity of each lorenz curve at each PT data point (m x p)
    q <- s <- y1 <- y2 <- w2 <- matrix(nrow = p, ncol = 3)
    dw <- dy <- matrix(nrow = p, ncol = 6) # distances between PT positions/intensities


    for (i in seq_len(p)) {

        # Phase 1: update peak intensities to reduce difference between original spectrum and sum of lorentz curves
        s[i, ] <- rowSums(Z[rnk[i], ]) # Sum of lorentz curves (cols) at current peak positions (rows)
        q[i, ] <- y0[i, ] / s[i, ] # Proportion between original spectrum and sum of lorenz curves
        y1[i, ] <- Z[rnk[i], i] * q[i, ] # New intensity of the peak triplet
        w1[i, ] <- w0[i, ] # Positions stay the same in this phase

        # Phase 2: replace shoulders with mirrored points
        if (y1[i, 1] < y1[i, 2] && y1[i, 2] < y1[i, 3]) { # ascending shoulder
            w2[i, 3] <- 2 * w0[i, 2] - w0[i, 1]
            y2[i, 3] <- y1[i, 1]
        } else if (y1[i, 1] > y1[i, 2] && y1[i, 2] > y1[i, 3]) { # descending shoulder
            w2[i, 1] <- 2 * w0[i, 2] - w0[i, 3]
            y2[i, 1] <- y1[i, 3]
        }

        # Phase 3: Update parameters of lorentz curves based on updated peak intensities
        dw[i] <- calc_distances_in_x(w2)
        dy[i] <- calc_distances_in_y(y2)
        w[i] <- calc_w_v15(dw, dy)
        lambda[i] <- calc_lambda_v15(dw[i], dy[i], w[i])
        A[i] <- calc_A_v15(dw[i], dy[i], lambda[i], w[i])
    }

    y_normed <- y / sum(y)
    yhat <- rowSums(Y)
    yhat_normed <- yhat / sum(yhat)
    se_normed <- (y_normed - yhat_normed)^2
    mse_normed <- (1 / length(se_normed)) * sum(se_normed)
    cat3(sprintf("Normed MSE: %.22f", mse_normed))

    return(list())
}

calc_distances_in_x <- function(...) {
    # TODO
}

calc_distances_in_y <- function(...) {
    # TODO
}

calc_w_v15 <- function(...) {
    # TODO
}

calc_lambda_v15 <- function(...) {
    # TODO
}

calc_A_v15 <- function(...) {
    # TODO
}
