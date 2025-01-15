testthat::skip_on_cran()
testthat::skip_on_ci()

#' @noRd
#' @title Setup a development environment for `plot_spectrum`
mkenv_plot_spectrum <- function() {
    target <- c("sim1", "sap")[1]
    decon <- switch(
        target,
        "sim1" = deconvolute(sim[[1]], sfr = c(3.55, 3.35)),
        "sap" = deconvolute(sap[[1]], sfr = c(3.2, -3.2), smopts = c(1, 3), delta = 3)
    )
    args <- stub(func = plot_spectrum, x = decon, ... = NULL, envir = .GlobalEnv)
    deferred_run()
    if (FALSE) plot_spectrum(decon, frame = TRUE)
}

#' @noRd
#' @examples
#' test_plot_spectrum(1, 2) # first two plots
#' test_plot_spectrum(2:4) # second to fourth plot
test_plot_spectrum <- function(figs = 1:6) {
    if (environment() %==% .GlobalEnv) {
        figs <- 1:6
        deferred_run()
    }
    n <- length(figs)
    nr <- ceiling(sqrt(n))
    nc <- if ((nr - 1) * nr >= n) nr - 1 else nr
    spec <- read_spectrum(metabodecon_file("sim/sim_01"))
    decon <- generate_lorentz_curves_sim(spec)
    local_par(mfrow = c(nr, nc))

    # Plot the full (non-deconvoluted) spectrum
    if (1 %in% figs) plot_spectrum(
        spec,
        sub1 = FALSE
    )

    # Remove connecting lines, and focus on a specific region specified in ppm
    if (2 %in% figs) plot_spectrum(
        decon,
        foc_rgn = c(3.49, 3.45),
        con_lines = FALSE
    )

    # Show second derivative and focus on a specific region specified as fraction
    if (3 %in% figs) plot_spectrum(
        decon,
        sub2 = TRUE,
        foc_frac = c(0.40, 0.30)
    )

    # Change color of focus rectangle and margins of sub figure 1
    if (4 %in% figs) plot_spectrum(
        decon,
        sub1 = list(mar = c(3, 6, 3, 6), lt_axis = list(col = "violet")),
        foc_rect = list(border = "violet", col = transp("violet")),
        con_lines = list(col = "violet")
    )

    # Hide xlab and show second derivative
    if (5 %in% figs) plot_spectrum(
        decon,
        sub2 = TRUE,
        sub3 = list(bt_axis = list(text = "")),
        frame = TRUE,
        con_lines = FALSE
    )

    # Change the figure region for sub figure 1
    if (6 %in% figs) plot_spectrum(
        decon,
        sub1 = list(
            fig_rgn_npc = c(0, 1.0, 0.3, 1.0),
            mar = c(0, 5, 0, 0)
        )
    )

    # Test all combinations of (de-)activated sub figures
    if (7 %in% figs)  plot_spectrum(decon, sub1 = FALSE, sub2 = FALSE, sub3 = FALSE, frame = TRUE)
    if (8 %in% figs)  plot_spectrum(decon, sub1 = TRUE,  sub2 = FALSE, sub3 = FALSE)
    if (9 %in% figs)  plot_spectrum(decon, sub1 = FALSE, sub2 = TRUE,  sub3 = FALSE)
    if (10 %in% figs) plot_spectrum(decon, sub1 = FALSE, sub2 = FALSE, sub3 = TRUE)
    if (11 %in% figs) plot_spectrum(decon, sub1 = TRUE,  sub2 = TRUE,  sub3 = FALSE)
    if (12 %in% figs) plot_spectrum(decon, sub1 = TRUE,  sub2 = FALSE, sub3 = TRUE)
    if (13 %in% figs) plot_spectrum(decon, sub1 = FALSE, sub2 = TRUE,  sub3 = TRUE)
    if (14 %in% figs) plot_spectrum(decon, sub1 = TRUE,  sub2 = TRUE,  sub3 = TRUE)
}

test_grconvert <- function() {
    par(mfrow = c(1, 2), xpd = TRUE)
    plot_dummy()
    x <- c(0.25, 0.75)
    xds <- list()
    units <- c("in", "dev", "ndc", "nfc", "npc", "nic", "lines", "chars")
    for (i in seq_along(units)) {
        y <- 0.1 * i
        unit <- units[i]
        xu <- grconvertX(x, from = unit, to = "user")
        xds[[i]] <- grconvertX(x, from = unit, to = "ndc")
        lines(x = xu, y = rep(y, 2))
        xur <- collapse(round(xu, 2), " - ")
        label <- sprintf("0.25 - 0.75 %s == %s %s", unit, xur, "user")
        text(x = 0.5, y = y, labels = label, pos = 3)
    }
    with_fig(
        fig = c(0, 1, 0, 1),
        expr = {
            plot_empty() # important to setup user coordinates
            for (i in seq_along(units)) {
                y <- 0.1 * i + 0.05
                xu <- grconvertX(xds[[i]], from = "ndc", to = "user")
                lines(x = xu, y = rep(y, 2), col = "red")
            }
        }
    )
    par(mfrow = c(1, 1), xpd = FALSE)
}

test_result <- test_that("plot_spectrum works", {
    tmp <- vdiffr::expect_doppelganger(
        title = "plot_spectrum",
        fig = test_plot_spectrum,
        writer = function(plot, file, title = "") {
            with_svg(file, plot(), width = 12, height = 16)
        }
    )
})
