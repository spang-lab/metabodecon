# styler: off
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

    lci_v13 <- init_lorentz_curves_v13(spec)
    lci_v14 <- init_lorentz_curves_v14(spec)
    compare_lci_v13_with_lci_v14(lci_v13, lci_v14)

    x <- spec$sdp; y <- spec$y_smooth; p <- spec$peak
    pl <- p$right[p$high]; pc <- p$center[p$high]; pr <- p$left[p$high];
    Y <- lci_v14$Y
    Yt <- t(Y)
    Z <- lci_v14$Z

    lc1_v10 <- refine_lc_internal_v10(x, y, pl, pc, pr, Yt)
    lc1_v12 <- refine_lc_internal_v12(x, y, pl, pc, pr, Y)
    lc1_v14 <- refine_lc_internal_v14(spec, Y, Z)
    compare_lc1_v10_with_lc1_v12(lc1_v10, lc1_v12)
    compare_lc1_v12_with_lc1_v14(lc1_v12, lc1_v14)

    Y <- lc1_v12$Y
    Yt <- t(Y)
    Z <- lc1_v14$Z

    lc2_v10 <- refine_lc_internal_v10(x, y, pl, pc, pr, Yt)
    lc2_v12 <- refine_lc_internal_v12(x, y, pl, pc, pr, Y)
    lc2_v14 <- refine_lc_internal_v14(spec, Y, Z)
    compare_lc1_v10_with_lc1_v12(lc2_v10, lc2_v12)
    compare_lc1_v12_with_lc1_v14(lc2_v12, lc2_v14)
}

compare_lci_v13_with_lci_v14 <- function(lci_v13, lci_v14) {

    vcomp(lci_v14[["P"]][["il"]],  lci_v13[["P1"]][["pr"]] )
    vcomp(lci_v14[["P"]][["ic"]],  lci_v13[["P1"]][["pc"]] )
    vcomp(lci_v14[["P"]][["ir"]],  lci_v13[["P1"]][["pl"]] )
    vcomp(lci_v14[["P"]][["xl"]],  lci_v13[["P2"]][["wr"]] )
    vcomp(lci_v14[["P"]][["xc"]],  lci_v13[["P2"]][["wc"]] )
    vcomp(lci_v14[["P"]][["xr"]],  lci_v13[["P2"]][["wl"]] )
    vcomp(lci_v14[["P"]][["yl"]],  lci_v13[["P2"]][["yr"]] )
    vcomp(lci_v14[["P"]][["yc"]],  lci_v13[["P2"]][["yc"]] )
    vcomp(lci_v14[["P"]][["yr"]],  lci_v13[["P2"]][["yl"]] )
    vcomp(lci_v14[["P"]][["ds"]],  lci_v13[["P2"]][["i"]]  )
    vcomp(lci_v14[["P"]][["as"]],  lci_v13[["P2"]][["j"]]  )
    vcomp(lci_v14[["D"]][["wl"]],  lci_v13[["P3"]][["wr"]] )
    vcomp(lci_v14[["D"]][["wc"]],  lci_v13[["P3"]][["wc"]] )
    vcomp(lci_v14[["D"]][["wr"]],  lci_v13[["P3"]][["wl"]] )
    vcomp(lci_v14[["D"]][["ycl"]], lci_v13[["P3"]][["ycr"]])
    vcomp(lci_v14[["D"]][["yrl"]], lci_v13[["P3"]][["ylr"]])
    vcomp(lci_v14[["D"]][["yrc"]], lci_v13[["P3"]][["ylc"]])
    vcomp(lci_v14[["D"]][["wcl"]], lci_v13[["P3"]][["xcr"]])
    vcomp(lci_v14[["D"]][["wrl"]], lci_v13[["P3"]][["xlr"]])
    vcomp(lci_v14[["D"]][["wrc"]], lci_v13[["P3"]][["xlc"]])
    vcomp(lci_v14[["P"]][["xr"]],  lci_v13[["P3"]][["wd"]] )
    vcomp(lci_v14[["A"]],          lci_v13[["A"]]          )
    vcomp(lci_v14[["lambda"]],     lci_v13[["lambda"]]     )
    vcomp(lci_v14[["w"]],          lci_v13[["w"]]          )
}


compare_lc1_v10_with_lc1_v12 <- function(lc1_v10, lc1_v12) {
    n <- nrow(lc1_v12$Y)
    p <- ncol(lc1_v12$Y)
    vcomp(lc1_v12$wl1,             lc1_v10$w_1_bak1                   )
    vcomp(lc1_v12$wc1,             lc1_v10$w_2_bak1                   )
    vcomp(lc1_v12$wr1,             lc1_v10$w_3_bak1                   )
    vcomp(lc1_v12$yl1,             lc1_v10$y_1_bak1                   )
    vcomp(lc1_v12$yc1,             lc1_v10$y_2_bak1                   )
    vcomp(lc1_v12$yr1,             lc1_v10$y_3_bak1                   )
    vcomp(lc1_v12$sl,              lc1_v10$sum_left                   )
    vcomp(lc1_v12$sc,              lc1_v10$sum_peaks                  )
    vcomp(lc1_v12$sr,              lc1_v10$sum_right                  )
    vcomp(lc1_v12$ql,              lc1_v10$proportion_left            )
    vcomp(lc1_v12$qc,              lc1_v10$proportion_peaks           )
    vcomp(lc1_v12$qr,              lc1_v10$proportion_right           )
    vcomp(lc1_v12$yl2,             lc1_v10$y_1_bak2                   )
    vcomp(lc1_v12$yc2,             lc1_v10$y_2_bak2                   )
    vcomp(lc1_v12$yr2,             lc1_v10$y_3_bak2                   )
    vcomp(lc1_v12$ia,              lc1_v10$ascending_shoulders        )
    vcomp(lc1_v12$id,              lc1_v10$descending_shoulders       )
    vcomp(lc1_v12$wl2,             lc1_v10$w_1_bak2                   )
    vcomp(lc1_v12$wc2,             lc1_v10$w_2_bak2                   )
    vcomp(lc1_v12$wr2,             lc1_v10$w_3_bak2                   )
    vcomp(lc1_v12$yl3,             lc1_v10$y_1_bak3                   )
    vcomp(lc1_v12$yc3,             lc1_v10$y_2_bak3                   )
    vcomp(lc1_v12$yr3,             lc1_v10$y_3_bak3                   )
    vcomp(lc1_v12$yl3,             lc1_v10$y_1_bak3                   )
    vcomp(lc1_v12$yc3,             lc1_v10$y_2_bak3                   )
    vcomp(lc1_v12$yr3,             lc1_v10$y_3_bak3                   )
    vcomp(lc1_v12$yl4,             lc1_v10$y_1_bak4                   )
    vcomp(lc1_v12$yc4,             lc1_v10$y_2_bak4                   )
    vcomp(lc1_v12$yr4,             lc1_v10$y_3_bak4                   )
    vcomp(lc1_v12$wd,              lc1_v10$w_delta_new                )
    vcomp(lc1_v12$w,               lc1_v10$w_new                      )
    vcomp(lc1_v12$lambda,          lc1_v10$lambda_new                 )
    vcomp(lc1_v12$A,               lc1_v10$A_new                      )
    vcomp(lc1_v12$spectrum_approx, lc1_v10$spectrum_approx[1, ]       )
    vcomp(lc1_v12$mse_normed,      lc1_v10$mse_normed                 )
    vcomp(lc1_v12$Y[1, ],          lc1_v10$lorentz_curves_initial[, 1])
    vcomp(lc1_v12$Y[n, ],          lc1_v10$lorentz_curves_initial[, n])
    vcomp(lc1_v12$Y[, 1],          lc1_v10$lorentz_curves_initial[1, ])
    vcomp(lc1_v12$Y[, p],          lc1_v10$lorentz_curves_initial[p, ])
}

compare_lc1_v12_with_lc1_v14 <- function(lc1_v12, lc1_v14) {
    vcomp(lc1_v12$w,               lc1_v14$w     )
    vcomp(lc1_v12$lambda,          lc1_v14$lambda)
    vcomp(lc1_v12$A,               lc1_v14$A     )
}
