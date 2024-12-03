#' @noRd
#' @title Setup a development environment for `plot_spectrum`
mkenv_plot_spectrum <- function() {
    target <- c("sim1", "sap2")[2]
    decon <- switch(target,
        "sim1" = get_sim1_decon1(),
        "sap2" = get_sap2_idecon()
    )
    args <- stub(
        func = plot_spectrum,
        obj = decon,
        foc_rgn = "auto",
        envir = .GlobalEnv
    )
    width <- options("width")
    line <- paste0(collapse(rep("-", width), ""), "\n")
    cat(line)
    cat("Assigned to .GlobalEnv:\n")
    cat(line)
    str(args, 1, give.attr = FALSE)
    cat(line)
    invisible(args)
}

test_interactive <- function() {
    target <- c("sim1", "sap2")[2]
    decon <- switch(target,
        "sim1" = get_sim1_decon1(),
        "sap2" = get_sap2_idecon()
    )
    plot_spectrum(decon)
}

#' @noRd
#' @examples
#' test_plot_spectrum(1, 2) # first two plots
#' test_plot_spectrum(2:4) # second to fourth plot
test_plot_spectrum <- function(figs = 1:5) {
    n <- length(figs)
    nr <- ceiling(sqrt(n))
    nc <- if (nr^2 > n) nr - 1 else nr
    spec <- read_spectrum(metabodecon_file("sim/sim_01"))
    decon <- generate_lorentz_curves_sim(spec)
    opar <- par(mfrow = c(nr, nc))
    on.exit(par(opar))

    # Plot the full (non-deconvoluted) spectrum
    if (1 %in% figs) plot_spectrum(spec)

    # Focus on a specific region, specified in ppm or as fraction
    if (2 %in% figs) plot_spectrum(decon, foc_rgn = c(3.49, 3.45), foc_unit = "ppm")
    if (3 %in% figs) plot_spectrum(decon, foc_rgn = c(0.40, 0.30))

    # Change margin, color and position of the focus region
    if (4 %in% figs) {
        plot_spectrum(decon,
            sub_mar = c(4, 4, 0, 0),
            rct_fill = rgb(0.9, 0.5, 0.9, alpha = 0.1),
            rct_col = "violet",
            sub_rgn = c(x1 = 0.1, x2 = 0.9, y1 = 0.4, y2 = 0.9)
        )
    }

    # Remove connecting lines and fill colors
    if (5 %in% figs) {
        plot_spectrum(decon,
            rct_fill = NULL,
            sub_fill_col = NULL,
            cnct_col = NULL
        )
    }

    # Hide xlab and ylab
    if (6 %in% figs) {
        plot_spectrum(decon,
            xlab = "",
            ylab = "",
            mar = c(2, 2, 0, 1)
        )
    }
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
