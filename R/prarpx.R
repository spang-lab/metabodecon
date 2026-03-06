#' @noRd
#' @title Compute mdrb MSE for parameters
#' @description
#' Helper to evaluate MSE of mdrb deconvolution for a single spectrum.
mdrb_mse_for_params <- function(spectrum,
                                sfr,
                                smopts,
                                delta,
                                nfit) {
    check_mdrb(stop_on_fail = TRUE)
    spec <- mdrb::Spectrum$new(spectrum$cs, spectrum$si, sfr)
    dec <- mdrb::Deconvoluter$new()
    dec$set_moving_average_smoother(smopts[1], smopts[2])
    dec$set_noise_score_selector(delta)
    dec$set_analytical_fitter(nfit)
    decon <- dec$deconvolute_spectrum(spec)
    decon$mse()
}

#' @noRd
#' @title Evaluate PRARPX and MSE across a parameter grid
#' @description
#' Runs a grid of metabodecon parameters on a simulated spectrum and
#' returns PRARPX and MSE values for each combination.
prarpx_grid <- function(spectrum,
                        smopts_grid,
                        delta_grid,
                        nfit_grid,
                        sfr = NULL,
                        verbose = TRUE) {
    sfr <- aki_default_sfr(spectrum, sfr)
    truepar <- spectrum$meta$simpar
    n_grid <- length(smopts_grid) * length(delta_grid) * length(nfit_grid)
    rows <- vector("list", n_grid)
    idx <- 1
    for (sm in smopts_grid) {
        for (delta in delta_grid) {
            for (nfit in nfit_grid) {
                decon <- deconvolute(
                    x = spectrum,
                    nfit = nfit,
                    smopts = sm,
                    delta = delta,
                    sfr = sfr,
                    ask = FALSE,
                    verbose = FALSE,
                    use_rust = TRUE
                )
                prx <- calc_prarpx(decon, truepar = truepar)
                mse_val <- decon$mse$raw
                if (is.null(mse_val) || is.na(mse_val)) {
                    resid <- decon$si - decon$sit$sup
                    mse_val <- mean(resid * resid)
                }
                rows[[idx]] <- data.frame(
                    sm_iter = sm[1],
                    sm_win = sm[2],
                    delta = delta,
                    nfit = nfit,
                    n_peaks = length(decon$lcpar$x0),
                    prarpx = prx,
                    mse = mse_val
                )
                if (isTRUE(verbose)) {
                    logf("Grid %d/%d", idx, n_grid)
                }
                idx <- idx + 1
            }
        }
    }
    do.call(rbind, rows)
}

#' @noRd
#' @title Summarize best PRARPX and MSE rows
prarpx_best_rows <- function(df) {
    best_mse_idx <- which.min(df$mse)
    best_prx_idx <- which.max(df$prarpx)
    list(best_mse = df[best_mse_idx, , drop = FALSE],
         best_prarpx = df[best_prx_idx, , drop = FALSE])
}

#' @noRd
#' @title Plot PRARPX vs MSE
prarpx_plot <- function(df, path) {
    grDevices::png(path, width = 900, height = 420)
    opar <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    on.exit(graphics::par(opar), add = TRUE)
    graphics::plot(
        df$prarpx, df$mse,
        xlab = "PRARPX",
        ylab = "MSE",
        pch = 19,
        cex = 0.5,
        col = "grey"
    )
    best <- prarpx_best_rows(df)
    graphics::points(best$best_mse$prarpx, best$best_mse$mse,
                     pch = 19, col = "#1f77b4", cex = 1.2)
    graphics::points(best$best_prarpx$prarpx, best$best_prarpx$mse,
                     pch = 19, col = "#ff7f0e", cex = 1.2)
    graphics::legend(
        "topright",
        legend = c("Best MSE", "Best PRARPX"),
        col = c("#1f77b4", "#ff7f0e"),
        pch = 19,
        bty = "n"
    )

    mse_scaled <- df$mse
    rng <- range(mse_scaled, na.rm = TRUE)
    if (rng[1] != rng[2]) {
        mse_scaled <- (mse_scaled - rng[1]) / (rng[2] - rng[1])
    } else {
        mse_scaled <- rep(0, nrow(df))
    }
    graphics::plot(
        df$prarpx, mse_scaled,
        xlab = "PRARPX",
        ylab = "MSE (scaled)",
        pch = 19,
        cex = 0.5,
        col = "grey"
    )
    graphics::points(best$best_mse$prarpx, mse_scaled[which.min(df$mse)],
                     pch = 19, col = "#1f77b4", cex = 1.2)
    graphics::points(best$best_prarpx$prarpx, mse_scaled[which.max(df$prarpx)],
                     pch = 19, col = "#ff7f0e", cex = 1.2)
    grDevices::dev.off()
}
