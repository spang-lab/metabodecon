# Work in Progress #####

#' @title Plot Aligned Spectra
#' @description Plots the aligned and unaligned spectra for comparison.
#' @param YA Matrix. `YA[i,j]` == Signal Intensity of spectrum i at index j AFTER alignmnent.
#' @param YB Matrix. `YB[i,j]` == Signal Intensity of spectrum i at index j BEFORE alignment.
#' @param PA List of vectors. `PA[[i]][j]` == index of peak j of spectrum i AFTER alignment.
#' @param PB List of vectors. `PB[[i]][j]` == index of peak j of spectrum i BEFORE alignment.
#' @param mfcol Vector of two integers specifying the number of rows and columns of the plot grid
#' @examples
#' x <- seq(1.5 * pi, 9.5 * pi, length.out = 90)
#' y <- 10 * sin(x) # y without noise
#' p <- sort(order(y, decreasing = TRUE)[1:4]) # peaks without noise
#' Y <- lapply(1:4, function(i) smooth(smooth(y + rnorm(90)))) # add noise
#' Y <- do.call(rbind, Y)
#' Y <- Y - min(Y)
#' YB <- rbind(
#'     c(rep(0, 1), Y[1, ], rep(0, 9)),
#'     c(rep(0, 4), Y[1, ], rep(0, 6)),
#'     c(rep(0, 8), Y[1, ], rep(0, 2)),
#'     c(rep(0, 3), Y[1, ], rep(0, 7))
#' )
#' YA <- rbind(
#'     c(rep(0, 5), Y[1, ], rep(0, 5)),
#'     c(rep(0, 5), Y[1, ], rep(0, 5)),
#'     c(rep(0, 5), Y[1, ], rep(0, 5)),
#'     c(rep(0, 5), Y[1, ], rep(0, 5))
#' )
#' PA <- list(p + 5, p + 5, p + 5, p + 5)
#' PB <- list(p + 1, p + 4, p + 8, p + 3)
#' plot_aligned_spectra(YA, YB, PA, PB)
plot_aligned_spectra <- function(YA, YB, PA, PB, mfcol = c(nrow(YA), 1)) {
    s <- nrow(YA)
    if (!is.null(mfcol)) {
        opar <- par(mfcol = mfcol, mar = c(0, 2, 0, 0), oma = c(4.1, 2.1, 0, 0))
        on.exit(par(opar), add = TRUE)
    }
    for (i in seq_len(s)) {
        plot(x = seq_len(ncol(YB)), y = YB[i, ],
            type = "l", lty = 1, col = "darkgrey",
            xlim = c(1, ncol(YB)), ylim = c(0, max(YB)),
            xaxt = if (i == s) "s" else "n",
            xlab = "Datapoint Number", ylab = "Signal Intensity"
        )
        lines(x = seq_len(ncol(YA)), y = YA[i, ], col = "blue", lty = 1)
        abline(v = PB[[i]], col = transp("darkgrey", 0.5), lty = 2)
        abline(v = PA[[i]], col = transp("blue", 0.5), lty = 2)
    }
    mtext("Datapoint Number", side = 1, outer = TRUE, line = 3)
    mtext("Signal Intensity", side = 2, outer = TRUE, line = 1)
}
