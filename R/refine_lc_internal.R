test_refine_lc_internal <- function() {

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
    lci2 <- init_lorentz_curves_v13(x, y, pc, pl, pr)

    # Refine Lorentz Curves

    Yt <- t(Y)

    # Refine Lorentz Curves Internal
    system.time(v10 <- refine_lc_internal_v10(x, y, pl, pc, pr, Yt))
    system.time(v11 <- refine_lc_internal_v11(x, y, pl, pc, pr, Yt))
    system.time(v11 <- refine_lc_internal_v12(x, y, pl, pc, pr, Y))

    # Check correctness
    checks <- c()
    checks[1] <- isTRUE(all.equal(v11$wl1, v10$w_1_bak1))
    checks[2] <- isTRUE(all.equal(v11$wc1, v10$w_2_bak1))
    checks[3] <- isTRUE(all.equal(v11$wr1, v10$w_3_bak1))
    checks[4] <- isTRUE(all.equal(v11$yl1, v10$y_1_bak1))
    checks[5] <- isTRUE(all.equal(v11$yc1, v10$y_2_bak1))
    checks[6] <- isTRUE(all.equal(v11$yr1, v10$y_3_bak1))
    checks[7] <- isTRUE(all.equal(v11$sl,  v10$sum_left))
    checks[8] <- isTRUE(all.equal(v11$sc,  v10$sum_peaks))
    checks[9] <- isTRUE(all.equal(v11$sr,  v10$sum_right))
    checks[10] <- isTRUE(all.equal(v11$ql,  v10$proportion_left))
    checks[11] <- isTRUE(all.equal(v11$qc,  v10$proportion_peaks))
    checks[12] <- isTRUE(all.equal(v11$qr,  v10$proportion_right))
    checks[13] <- isTRUE(all.equal(v11$yl2, v10$y_1_bak2))
    checks[14] <- isTRUE(all.equal(v11$yc2, v10$y_2_bak2))
    checks[15] <- isTRUE(all.equal(v11$yr2, v10$y_3_bak2))
    checks[16] <- isTRUE(all.equal(v11$ia, v10$ascending_shoulders))
    checks[17] <- isTRUE(all.equal(v11$id, v10$descending_shoulders))
    checks[18] <- isTRUE(all.equal(v11$wl2, v10$w_1_bak2))
    checks[19] <- isTRUE(all.equal(v11$wc2, v10$w_2_bak2))
    checks[20] <- isTRUE(all.equal(v11$wr2, v10$w_3_bak2))
    checks[21] <- isTRUE(all.equal(v11$yl3, v10$y_1_bak3))
    checks[22] <- isTRUE(all.equal(v11$yc3, v10$y_2_bak3))
    checks[23] <- isTRUE(all.equal(v11$yr3, v10$y_3_bak3))
    checks[24] <- isTRUE(all.equal(v11$yl3, v10$y_1_bak3))
    checks[25] <- isTRUE(all.equal(v11$yc3, v10$y_2_bak3))
    checks[26] <- isTRUE(all.equal(v11$yr3, v10$y_3_bak3))
    checks[27] <- isTRUE(all.equal(v11$yl4, v10$y_1_bak4))
    checks[28] <- isTRUE(all.equal(v11$yc4, v10$y_2_bak4))
    checks[29] <- isTRUE(all.equal(v11$yr4, v10$y_3_bak4))
    checks[30] <- isTRUE(all.equal(v11$wd,  v10$w_delta_new))
    checks[31] <- isTRUE(all.equal(v11$w, v10$w_new))
    checks[32] <- isTRUE(all.equal(v11$lambda, v10$lambda_new))
    checks[33] <- isTRUE(all.equal(v11$A, v10$A_new))
    checks[34] <- isTRUE(all.equal(v11$spectrum_approx, v10$spectrum_approx))
    checks[35] <- isTRUE(all.equal(v11$mse_normed, v10$mse_normed))
    if (FALSE) isTRUE(all.equal(v11$Yt, v10$lorentz_curves_initial))
    return(checks)
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
