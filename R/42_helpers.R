# Get #####

get_foc_frac <- function(obj, foc_rgn = NULL) {
    stopifnot(is_num(obj$cs))
    stopifnot(is.null(foc_rgn) || (is_num(foc_rgn, 2)))
    if (is.null(foc_rgn)) {
        n <- length(obj$cs)
        width <- min(256 / n, 0.5)
        center <- if (n <= 2048) 0.5 else 0.75
        c(center - width, center + width)
    } else {
        convert_pos(foc_rgn, range(obj$cs), c(0, 1))
    }
}

get_foc_rgn <- function(obj, foc_frac = NULL) {
    stopifnot(is_num(obj$cs))
    stopifnot(is.null(foc_frac) || (is_num(foc_frac, 2)))
    if (is.null(foc_frac)) foc_frac <- get_foc_frac(obj)
    quantile(obj$cs, foc_frac)
}

get_sub_fig_args <- function(obj, foc_frac, foc_rgn, sub1, sub2, sub3, dot_args) {
    if ((n <- length(dot_args)) != (k <- length(names(dot_args)))) {
        if (n == 1) dot_args <- dot_args[[1]]
        else stop("... must contain tag=value pairs or a single list of such pairs")
    }
    usr_args <- lapply(list(sub1, sub2, sub3), function(x)
        modifyList(dot_args, x)
    )
    layout <- paste0(
        if (isFALSE(sub1$show)) "_" else "1",
        if (isFALSE(sub2$show)) "_" else "2",
        if (isFALSE(sub3$show)) "_" else "3"
    )
    fig_rgns <- switch(layout,           # x1 x2 y1    y2
        "___" = list(c(0.00, 0.00, 0.00, 0.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 0.00, 0.00, 0.00)),
        "1__" = list(c(0.00, 1.00, 0.00, 1.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 0.00, 0.00, 0.00)),
        "_2_" = list(c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.00, 1.00), c(0.00, 0.00, 0.00, 0.00)),
        "__3" = list(c(0.00, 1.00, 0.00, 1.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.00, 1.00)),
        "12_" = list(c(0.00, 1.00, 0.20, 1.00), c(0.00, 1.00, 0.00, 0.20), c(0.00, 0.00, 0.00, 0.00)),
        "1_3" = list(c(0.00, 1.00, 0.40, 1.00), c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.00, 0.20)),
        "_23" = list(c(0.00, 0.00, 0.00, 0.00), c(0.00, 1.00, 0.50, 1.00), c(0.00, 1.00, 0.00, 0.30)),
        "123" = list(c(0.00, 1.00, 0.50, 1.00), c(0.00, 1.00, 0.30, 0.50), c(0.00, 1.00, 0.00, 0.20)),
        stop("Invalid layout")
    )
    mars <- switch(layout,          # b  l  t  r
        "1_3" = list(c(1, 6, 2, 2), c(0, 0, 0, 0), c(0, 0, 1, 0)),
        "_23" = list(c(0, 0, 0, 0), c(1, 2, 2, 2), c(0, 0, 1, 0)),
        "123" = list(c(0, 6, 2, 2), c(1, 6, 0, 2), c(0, 0, 1, 0)),
        list(c(0, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0))
    )
    for (i in 1:3) {
        fig_rgns[[i]] <- usr_args[[i]]$fig_rgn_npc %||% fig_rgns[[i]]
        usr_args[[i]]$fig_rgn_npc <- NULL
    }
    def_args <- lapply(1:3, function(i) list(
        obj = obj,
        foc_rgn = foc_rgn,
        foc_frac = foc_frac,
        add = TRUE,
        fig_rgn = npc_to_ndc(fig_rgns[[i]]),
        mar = mars[[i]]
    ))
    def_args[[1]]$bt_axis <- if (sub2$show) FALSE else if (sub3$show) list(text = "")
    def_args[[2]]$lt_axis <- if (layout == "_2_") list() else FALSE
    def_args[[2]]$bt_axis <-if (sub3$show) list(text = "")
    def_args[[2]]$show_d2 <- TRUE
    def_args[[3]]$bg_rect <- if (layout != "__3") list(border = NA)
    def_args[[3]]$foc_only <- FALSE
    def_args[[3]]$foc_rect <- if (layout == "__3") FALSE
    def_args[[3]]$lgd <- FALSE
    def_args[[3]]$lt_axis <- if (layout != "__3") list(text = ""  )
    structure(
        mapply(modifyList, def_args, usr_args, SIMPLIFY = FALSE),
        names = c("sub1", "sub2", "sub3")
    )
}

get_ds_def_args <- function(show_d2 = FALSE, foc_only = TRUE) {
    stopifnot(is_bool(show_d2), is_bool(foc_only))
    show_si  <- !show_d2
    ylab <- if (show_si) "Signal Intensity [au] / 1e6" else "Second Derivative / 1e6"
    lgdx <- if (show_si) "topright" else "bottomright"
    list(
        lgd       = list(show = TRUE, x = lgdx, bg = "white"),
        d2_line   = list(show = show_d2),
        si_line   = list(show = show_si, col = "black"),
        sm_line   = list(show = show_si, col = "blue"),
        sp_line   = list(show = show_si, col = "red"),
        lc_lines  = list(show = show_si && foc_only, col = "darkgrey"),
        tp_lines  = list(show = FALSE, col = "darkgreen"),
        tp_verts  = list(show = FALSE, col = "darkgreen"),
        lc_verts  = list(show = FALSE, col = "darkgrey"),
        cent_pts  = list(show = show_si && foc_only, col = "blue", bg = "blue", pch = 24),
        bord_pts  = list(show = FALSE, col = "blue", pch = 124),
        norm_pts  = list(show = FALSE, col = "black", pch = 124),
        bg_rect   = list(show = TRUE, col = "white"),
        foc_rect  = list(show = TRUE, col = transp("yellow")),
        lc_rects  = list(show = FALSE, col = transp("black"), border = NA),
        tp_rects  = list(show = FALSE, col = transp("darkgreen", 0.12), border = NA),
        bt_axis   = list(show = TRUE, text = "Chemical Shift [ppm]"),
        lt_axis   = list(show = TRUE, text = ylab, sf = 1e6, las  = 1),
        tp_axis   = list(show = FALSE),
        rt_axis   = list(show = FALSE)
    )
}

combine <- function(defaults, x) {
    name <- deparse(substitute(x))
    x <- if (is.null(x)) list()
        else if (isFALSE(x)) list(show = FALSE)
        else if (isTRUE(x)) list(show = TRUE)
        else if (is.list(x)) x
        else stopf("%s must be a bool or list", name)
    modifyList(defaults, x)
}

# Convert #####

usr_to_ndc <- function(usr = par("usr")) {
    x_ndc <- grconvertX(usr[1:2], from = "user", to = "ndc")
    y_ndc <- grconvertY(usr[3:4], from = "user", to = "ndc")
    c(x_ndc, y_ndc)
}

ndc_to_usr <- function(ndc) {
    x_usr <- grconvertX(ndc[1:2], from = "ndc", to = "user")
    y_usr <- grconvertY(ndc[3:4], from = "ndc", to = "user")
    c(x_usr, y_usr)
}

npc_to_ndc <- function(npc = c(0, 1, 0, 1)) {
    if (is.null(npc)) return(NULL)
    x_ndc <- grconvertX(npc[1:2], from = "npc", to = "ndc")
    y_ndc <- grconvertY(npc[3:4], from = "npc", to = "ndc")
    c(x_ndc, y_ndc)
}

# Multifigure #####

#' @noRd
#'
#' @title Set Figure Region
#'
#' @description
#' Calling `par(fig=xxyy)` resets the current multi-figure configuration (MFC)
#' to one row and one column. `set_fig()` handles this scenario, by first
#' storing the current MFC, then calling `par(fig=xxyy)` and finally returning a
#' function that can be used to restore the MFC. See 'Details' for further
#' information.
#'
#' @param fig Region to draw into, given as normalized device coordinates.
#'
#' @param add If TRUE, the new plot is added to the existing plot.
#'
#' @return Function to reset the MFC.
#'
#' @details
#' Note 1: Setting `par(fig=xxyy)` resets the current MFC to one row and one
#' column. I.e., array layouts defined by setting `mfrow` or `mfcol`, must be
#' saved before changing `fig` and restored after the plot has been drawn. The
#' same applies to the current figure number `mfg` (see Note 2).
#'
#' Note 2: When restoring `mfg` it's important to additionally advance one
#' frame, because when querying `mfg`, the "current figure number" is returned,
#' but when setting `mfg`, the value is interpreted as "next figure number".
#' I.e. without advancing one frame, the next figure would be drawn into the
#' spot of the current figure.
#'
#' Note 3: If `fig=xxyy` and `add=FALSE`, it is still necessary to use
#' `par(fig=xxyy, new=TRUE)` to prevent the device from being cleared (which
#' would be bad in a multi-figure environment). I.e., the figure number must be
#' advanced manually, as would have been done by a normal plot. This manual
#' increment must be done before the MFC is reset by `par(fig=xxyy, new=TRUE)`.
#' @examples
#' plot_dummy <- function() {
#'     plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
#'     text(0.5, 0.5, "dummy")
#' }
#' p <- local({
#'     opar <- par(mfrow = c(2, 2), mar = c(2, 2, 0.5, 0.5))
#'     on.exit(par(opar))
#'
#'     topleft <- plot_dummy()
#'     reset_mfc <- set_fig(fig = c(0.25, 0.50, 0.50, 0.75), add = TRUE)
#'     topleft2 <- plot_dummy()
#'     reset_mfc()
#'
#'     topright <- plot_dummy()
#'
#'     reset_mfc <- set_fig(fig = c(0.1, 0.4, 0.1, 0.4), add = FALSE)
#'     bottom_left <- plot_dummy()
#'     reset_mfc()
#'
#'     bottom_right <- plot_dummy()
#' })
set_fig <- function(fig = NULL, add = TRUE) {
    if (is.null(fig)) return(function() {}) # Nothing to do if figure region is NULL
    if (isFALSE(add)) plot_empty() # Advance one frame if `add=FALSE` (Note 3)
    op <- par(c("mar", "mfrow", "mfcol", "mfg", "fig")) # Store MF conf (Note 1)
    byrow <- mf_filled_by_row() # Store MF conf (Note 1)
    par(fig = fig, new = TRUE) # Set new figure region (Note 3)
    reset_mfc <- function() {
        if (byrow) par(mfrow = op$mfrow) else par(mfcol = op$mfcol) # Restore MF Layout (Note 1)
        par(mfg = op$mfg) # Restore current figure number (Note 1)
        plot_empty() # Advance one frame (Note 2)
    }
    reset_mfc
}

#' @noRd
#' @title Plot into specific figure region
#' @description For Details see [set_fig()].
#' @examples
#' plot_dummy <- function() {
#'     plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xaxs = "i", yaxs = "i")
#'     text(0.5, 0.5, "dummy")
#' }
#' p <- local({
#'     opar <- par(mfrow = c(2, 2), mar = c(2, 2, 0.5, 0.5))
#'     on.exit(par(opar))
#'     topleft <- plot_dummy()
#'     topleft2 <- with_fig(fig = c(0.25, 0.50, 0.50, 0.75), add = TRUE, plot_dummy())
#'     topright <- plot_dummy()
#'     bottom_left <- with_fig(fig = c(0.1, 0.4, 0.1, 0.4), add = FALSE, plot_dummy())
#'     bottom_right <- plot_dummy()
#' })
with_fig <- function(expr, fig = NULL, pos = NULL, add = TRUE) {
    reset_mfc <- set_fig(fig = fig, add = add)
    on.exit(reset_mfc())
    expr
}

local_fig <- function(fig = NULL, add = TRUE, envir = parent.frame()) {
  reset_mfc <- set_fig(fig = fig, add = add)
  defer(reset_mfc(), envir = envir)
}

#' @noRd
#' @description
#' Returns TRUE if the current multi-figure gets filled by row, else FALSE.
#' @examples
#' callwith <- function(by_row = TRUE) {
#'     grid <- c(2, 2)
#'     opar <- if (by_row) par(mfrow = grid) else par(mfcol = grid)
#'     on.exit(opar)
#'     plot(1:10)
#'     by_row <- mf_filled_by_row()
#'     plot(1:5)
#'     by_row
#' }
#' x <- callwith(by_row = TRUE)
#' y <- callwith(by_row = FALSE)
#' stopifnot(isTRUE(x) && isFALSE(y))
mf_filled_by_row <- function() {
    mfg <- par("mfg")
    row <- mfg[1]
    col <- mfg[2]
    nrows <- mfg[3]
    ncols <- mfg[4]
    if (nrows == 1 || ncols == 1) {
        # In this case it doesn't matter, as the figure spots in the grid will
        # be filled top-to-bottom / left-to-right in both cases.
        return(TRUE)
    } else {
        # We have at least two rows AND two cols. So what we can do is to set
        # c(1, 1) as next figure, and then advance one frame. If we end up in
        # c(1, 2) we are row-oriented, if we end up in c(2, 1) we are
        # column-oriented. After doing this we can reset the current figure
        # number to the original value.
        par(mfg = c(1, 1, nrows, ncols))
        on.exit({
            par(mfg = mfg) # Reset current figure number
            plot_empty() # (1)
            # (1) When querying `mfg` we get the "current figure number". But
            # when setting it, we set the "next figure number". I.e. we need to
            # advance one frame, or we would set the "current figure number" as
            # "next figure number".
        })
        plot_empty() # Draw into c(1, 1)
        plot_empty() # Draw into c(1, 2) or c(2, 1)
        mfg2 <- par("mfg") # Query current position
        return(if (mfg2[1] == 1) TRUE else FALSE)
    }
}
