#' @export
#' @title Convert from unit A to unit B
#'
#' @description
#' Converts positions/widths from unit A to unit B. If the direction of units A
#' and B is reversed, the width's sign will be reversed as well. To keep widths
#' strictly positive, wrap the result with `abs()`.
#'
#' @param xa A numeric vector specifying widths/postions in unit A.
#' @param ya,yb A numeric vector specifying the positions of at least two points
#' in unit A / unit B.
#'
#' @return A numeric vector of values converted from unit A to unit B.
#'
#' @examples
#' ya <- c(244, 246, 248, 250, 252)
#' yb <- c(15, 10, 5, 0, -5)
#' convert_width(c(2, 4, 8), ya, yb)
#' convert_pos(c(247, 249), ya, yb)
convert_pos <- function(xa, ya, yb) {
    # Example:
    #     |-------------w----------|
    #     |----d----|              |
    # xa  |        247             |
    # ya  244   246 | 248   250   252
    # yb  15    10  | 5     0     -5
    # ==>
    # wa  =  ya[n] - ya[1] = 252 - 244   = 8
    # wb  =  yb[n] - yb[1] = (-5) - 15   = 20
    # da  =  xa    - ya[1] = 247 - 244   = 3
    # db  =  da/wa * wb    = (3/8) * -20 = -7.5 (because da/wa == db/wb)
    # xb  =  yb[1] + db    = 15 + (-7.5) = 7.5
    if (length(ya) != length(yb)) stop("ya and yb must have the same length.")
    n <- length(ya)
    wa <- ya[n] - ya[1]
    wb <- yb[n] - yb[1]
    da <- xa - ya[1]
    db <- (da / wa) * wb
    xb <- yb[1] + db
    xb
}

#' @export
#' @rdname convert_pos
convert_width <- function(xa, ya, yb) {
    n <- length(ya)
    wa <- ya[n] - ya[1]
    wb <- yb[n] - yb[1]
    xa * (wb / wa)
}

#' @export
#' @title Calculate the Width of a Numeric Vector
#' @description
#' Calculates the width of a numeric vector by computing the difference between
#' the maximum and minimum values in the vector.
#' @param x A numeric vector.
#' @return
#' The width of the vector, calculated as the difference between its maximum and
#' minimum values.
#' @examples
#' vec <- c(1, 3, 5, 7, 9)
#' width(vec)
width <- function(x) {
    diff(range(x))
}

#' @noRd
#' @description
#' Converts a vector of chemical shifts given in ppm to a vector of frequencies
#' in Hz.
#' @param cs Vector of chemical shifts in ppm.
#' @param fqref Frequency of the reference molecule in Hz.
in_hz <- function(cs, fqref) {
    fqmax <- fqref - (min(cs) * 1e-6 * fqref) # Highest frequency in Hz
    fqmin <- fqref - (max(cs) * 1e-6 * fqref) # Lowest frequency in Hz
    fq <- seq(fqmin, fqmax, length.out = length(cs))
    fq <- fq[rev(order(cs))] # (1)
    # (1) Sort frequencies in descending order of chemical shifts, as the
    #     highest chemical shift corresponds to the lowest frequency.
    fq
}

#' @noRd
#' @description Functions to convert between SDP and PPM in a backwards
#' compatible way, i.e. this function does the same mistakes and numerical
#' errors as the original conversion in MetaboDecon1D.
#' @param sfr_sdp Signal free region borders (SFR) in scaled data point numbers
#' (SDP).
#' @param sdp Scaled data point numbers for all datapoints.
#' @param ppm Chemical shifts (CS) in parts per million (PPM) for all datapoints.
#' datapoints.
#' @details See 'CHECK-2: Signal free region (SFR) calculation' in `TODOS.md`.
sfr_in_ppm_bwc <- function(sfr_sdp, sdp, ppm) {
    stopifnot(
        is.numeric(sfr_sdp), length(sfr_sdp) == 2,
        is.numeric(sdp), length(sdp) >= 5,
        is.numeric(ppm), length(ppm) == length(sdp)
    )
    n <- length(sdp)
    sdp_step <- diff(range(sdp)) / (n - 1)
    sdp_nstep <- diff(range(sdp)) / n
    ppm_step <- diff(range(ppm)) / (n - 1)
    ppm_nstep <- diff(range(ppm)) / n
    sfr_dp <- sfr_sdp / sdp_step
    sfr_ppm <- max(ppm) - (n + 1 - sfr_dp) * ppm_nstep
    sfr_ppm
}

#' @noRd
#' @rdname sfr_in_ppm_bwc
#' @param sfr_ppm Signal free region borders (SFR) in parts per million (PPM).
sfr_in_sdp_bwc <- function(sfr_ppm, ppm, sf) {
    sfr_left  <- max(sfr_ppm)
    sfr_right <- min(sfr_ppm)
    ppm_left  <- max(ppm)
    ppm_right <- min(ppm)
    n <- length(ppm)
    ppm_nstep <- (ppm_left - ppm_right) / n
    dp_left   <- (n + 1) - (ppm_left - sfr_left)  / ppm_nstep
    dp_right  <- (n + 1) - (ppm_left - sfr_right) / ppm_nstep
    sdp_left  <- dp_left  / sf[1]
    sdp_right <- dp_right / sf[1]
    c(sdp_left, sdp_right)
}
