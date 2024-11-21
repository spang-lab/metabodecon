# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Helpers for `get_sfr()`, `get_wshw()` and `filter_peaks()` #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @description Repeatedly ask the user to confirm/refine SFR borders.
#' @noRd
confirm_sfr <- function(x, sfr = c(11.44494, -1.8828)) {
    plot_sfr(x, sfr[1], sfr[2])
    sfr_ok <- get_yn_input("Signal free region correctly selected?")
    while (!sfr_ok) {
        get_border <- function(msg) get_num_input(msg, x$ppm_min, x$ppm_max)
        sfr[1] <- get_border("Choose another left border: [e.g. 12]")
        sfr[2] <- get_border("Choose another right border: [e.g. -2]")
        plot_sfr(x, sfr[1], sfr[2])
        sfr_ok <- get_yn_input("Signal free region correctly selected?")
    }
    sfr
}

#' @description Repeatedly ask the user to confirm/refine the WSHW.
#' @noRd
confirm_wshw <- function(x, wshw) {
    plot_ws(x, wshw)
    ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    while (!ws_ok) {
        wshw <- get_num_input("Choose another half width range (in ppm) for the water artefact: [e.g. 0.1222154]")
        plot_ws(x, wshw)
        ws_ok <- get_yn_input("Water artefact fully inside red vertical lines?")
    }
    wshw
}

#' @noRd
#' @description
#' Takes the SFR in PPM and returns the SFR in PPM, DP and SDP.
#' @note
#' Because the conversion from PPM to DP/SDP is slightly off (by 1-2 data
#' points), the SFR borders in DP/SDP returned by this function are also
#' incorrect. However, to maintain backwards compatibility with the old
#' MetaboDecon1D function, the behaviour is not changed in this function.
#' Instead, to only work with the correct ppm values, set `bwc = 2` in
#' [filter_peaks()]. For details see `CHECK-2: signal free region calculation`
#' in `TODOS.md`.
enrich_sfr <- function(sfr, x) {
    stopifnot(is_ispec(x) || is_idecon(x))
    left_ppm <- sfr[1]
    right_ppm <- sfr[2]
    left_dp <- (x$n + 1) - (x$ppm_max - left_ppm) / x$ppm_nstep
    left_sdp <- left_dp / x$sf[1]
    right_dp <- (x$n + 1) - (x$ppm_max - right_ppm) / x$ppm_nstep
    right_sdp <- right_dp / x$sf[1]
    named(left_ppm, right_ppm, left_dp, right_dp, left_sdp, right_sdp)
}

#' @noRd
#' @description
#' Calculates the WSR in dp and ppm from the WSHW in ppm.
#' @note
#' Because the conversion from PPM to DP/SDP is slightly off (by 1-2 data
#' points), the SFR borders in DP/SDP returned by this function are also
#' incorrect. However, to maintain backwards compatibility with the old
#' MetaboDecon1D function, the behaviour is not changed in this function.
#' Instead, to only work with the correct ppm values, set `bwc = 2` in
#' [rm_water_signal()]. For details see `CHECK-3: water signal calculation` in
#' `TODOS.md`.
enrich_wshw <- function(wshw, x) {
    stopifnot(is_ispec(x) || is_idecon(x))
    x <- as_ispec(x)
    hwidth_ppm <- wshw
    hwidth_dp <- hwidth_ppm / x$ppm_nstep
    center_dp <- x$n / 2
    right_dp <- center_dp + hwidth_dp
    left_dp <- center_dp - hwidth_dp
    center_ppm <- x$ppm[center_dp]
    right_ppm <- x$ppm[right_dp]
    left_ppm <- x$ppm[left_dp]
    if (left_dp <= 1 || right_dp >= x$n) stop("WSR is out of range")
    named(
        left_ppm, right_ppm, center_ppm, hwidth_ppm,
        left_dp, right_dp, center_dp, hwidth_dp
    )
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Helpers for `find_peak()` #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

calc_second_derivative <- function(y) {
    n <- length(y)
    x <- c(NA, y[-n]) # x[i] == y[i-1]
    z <- c(y[-1], NA) # z[i] == y[i+1]
    d <- x + z - 2 * y
    d
}

get_right_border <- function(j, d, m) {
    r <- j + 1
    while (r < m) { # use r<m instead of r<=m because c4 requires d[r+1]
        c1 <- d[r] > d[r - 1]
        c2 <- d[r] >= d[r + 1]
        c3 <- d[r] < 0
        c4 <- d[r + 1] >= 0
        is_right_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_right_border)) {
            return(r)
        }
        r <- r + 1
    }
    return(NA)
}

get_left_border <- function(j, d) {
    l <- j - 1
    while (l > 1) { # use l>1 instead of l>=1 because c4 requires d[l-1]
        c1 <- d[l] > d[l + 1]
        c2 <- d[l] >= d[l - 1]
        c3 <- d[l] < 0
        c4 <- d[l - 1] >= 0
        is_left_border <- (c1 && c2) || (c1 && c3 && c4)
        if (isTRUE(is_left_border)) {
            return(l)
        }
        l <- l - 1
    }
    return(NA)
}

#' @noRd
#'
#' @title Get Peak Score
#'
#' @description
#' Calculate the score of a peak based on the sum of absolute second derivative
#' values of its datapoints.
#'
#' @param j <- Index of the peak center
#' @param l <- Index of the left border
#' @param r <- Index of the right border
#' @param a <- Absolute values of the second derivative for all data points
#'
#' @return The score of the peak.
#'
#' @examples
#'
#' #      ____________________________________________________
#' #     |____2________5___________9_______12____14____16_____|
#' #     |             x                                      |
#' #     |          x  x  x  x                    x           |
#' #     |       x  x  x  x  x  x              x  x           |
#' #     |_.__x__x__x__x__x__x__x__x__.__.__x__x__x__x__x__.__|
#' #          |----2---|-----4-----|        |--3--|--6--|
#'
#' y <- c( 0, 1, 2, 3, 4, 3, 3, 2, 1, 0, 0, 1, 2, 3, 1, 1, 0  )
#' a <- c(NA, 0, 0, 0, 2, 1, 1, 0, 0, 1, 1, 0, 0, 3, 2, 1, NA )
#' all.equal(a, abs(calc_second_derivative(y)))
#'
#' s1 <- get_peak_score( 5, 2,   9, a)
#' s2 <- get_peak_score(14, 12, 16, a)
#' stopifnot(s1 == min(sum(a[2:5]), sum(a[5:9])))
#' stopifnot(s2 == min(sum(a[12:14]), sum(a[14:16])))
get_peak_score <- function(j, l, r, a) {
    if (any(is.na(a[c(l, j, r)]))) {
        0
    } else {
        min(sum(a[l:j]), sum(a[j:r]))
    }
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# Helpers for `fit_lorentz_curves()` #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @title Initialize Lorentz Curve Parameters
#' @param spec List with elements: `x`, `y`, `peak` where `peak` is a list with
#' elements `center`, `left`, `right` and `high`.
#' @noRd
init_lc <- function(spec, verbose = TRUE) {

    # Init values
    p <- spec$peak
    ir <- p$right[p$high]  #
    ic <- p$center[p$high] # Index of each peak triplet position (PTP)
    il <- p$left[p$high]   #
    lmr <- sort(unique(c(il, ic, ir))) # Combined PTP indices
    rr <- match(ir, lmr) #
    rc <- match(ic, lmr) # Rank of each PTP
    rl <- match(il, lmr) #
    x <- spec$sdp; y <- spec$y_smooth # X and Y value for each data point
    xlmr <- x[lmr]; # X value for each PTP
    yr <- y[ir]; yc <- y[ic]; yl <- y[il]; # Intensity of each PTP
    xr <- x[ir]; xc <- x[ic]; xl <- x[il]; # Position of each PTP

    # Replace shoulders
    as <- (yr > yc) & (yc > yl) # Ascending shoulders (AS)
    ds <- (yr < yc) & (yc < yl) # Descending shoulders (DS)
    xl[ds] <- 2 * xc[ds] - xr[ds] # Replace DS with mirrored points (MP)
    xr[as] <- 2 * xc[as] - xl[as] # Replace AS with MP
    yl[ds] <- yr[ds] # Replace DS with MP
    yr[as] <- yl[as] # Replace AS with MP

    # Calculate distances
    wr  <- xr - xr #
    wc  <- xc - xr # Express positions wr/wc/wl as "distance to right border"
    wl  <- xl - xr #
    wrc <- wr - wc; wrl <- wr - wl; wcl <- wc - wl # x - distance between PTPs
    yrc <- yr - yc; yrl <- yr - yl; ycl <- yc - yl # y - distance between PTPs

    # Estimate parameters
    w <- calc_w(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr)
    lambda <- calc_lambda(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr)
    A <- calc_A(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, lambda, w, xr)

    # Calculate contribution of each lorentz curve to each PTP data point
    Z <- matrix(0, nrow = length(lmr), ncol = length(ic)) # 3614 x 1227 urine_1
    idx_A_non_zero <- which(A != 0)
    for (j in idx_A_non_zero) {
        Z[, j] <- abs(A[j] * (lambda[j] / (lambda[j]^2 + (xlmr - w[j])^2)))
    }

    # Print MSE
    mse <- mean((y[lmr] - rowSums(Z))^2)
    if (verbose) logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, as, ds)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

refine_lc_v14 <- function(spec, Z) {

    # Init x and y values
    x <- spec$sdp; y <- spec$y_smooth # x and y value for each data point

    # Init peak related variables
    p <- spec$peak
    ir <- p$right[p$high]; ic <- p$center[p$high]; il <- p$left[p$high] # index of each peak triplet position (PTP)
    lmr <- sort(unique(c(il, ic, ir))) # combined PTP indices
    rr <- match(ir, lmr);  rc <- match(ic, lmr);   rl <- match(il, lmr) # rank  of each PTP
    xlmr <- x[lmr]; # x value for each PTP
    yr <- y[ir]; yc <- y[ic]; yl <- y[il]; # intensity of each PTP
    xr <- x[ir]; xc <- x[ic]; xl <- x[il]; # position of each PTP
    sl <- sc <- sr <- numeric(length(ic)); # sum of lorentz curves (SLC) at each PTP
    ql <- qc <- qr <- numeric(length(ic)); # ratio (SLC / original spectrum) at each PTP

    # Init distance related variables
    wr  <- wc  <- wl  <- numeric(length(ic));
    wrc <- wrl <- wcl <- numeric(length(ic));
    yrc <- yrl <- ycl <- numeric(length(ic));

    # Init lorentz curves parameters and matrices
    A <- lambda <- w <- numeric(length(ic))

    for (i in seq_along(il)) {

        # Calculate the sum of all lorentz curves for each data point
        sl[i] <- sum(Z[rl[i], ])
        sc[i] <- sum(Z[rc[i], ])
        sr[i] <- sum(Z[rr[i], ])

        # Calculate the proportion between original spectrum an the sum of the
        # lorentz curves for each peak triplets position
        ql[i] <- yl[i] / sl[i]
        qc[i] <- yc[i] / sc[i]
        qr[i] <- yr[i] / sr[i]

        # Calculate the new heights of the peak triplets
        yl[i] <- Z[rl[i], i] * ql[i]
        yc[i] <- Z[rc[i], i] * qc[i]
        yr[i] <- Z[rr[i], i] * qr[i]

        # Calculate mirrored points for ascending and descending shoulders
        if ((yl[i] < yc[i]) && (yc[i] < yr[i])) { # Ascending shoulder
            xr[i] <- 2 * xc[i] - xl[i]
            yr[i] <- yl[i]
        }
        if ((yl[i] > yc[i]) && (yc[i] > yr[i])) { # Descending shoulder
            xl[i] <- 2 * xc[i] - xr[i]
            yl[i] <- yr[i]
        }

        # Calculate distances between peak triplet positions and intensities
        wr[i] <- xr[i] - xr[i]
        wc[i] <- xc[i] - xr[i]
        wl[i] <- xl[i] - xr[i]
        wrc[i] <- wr[i] - wc[i]
        wrl[i] <- wr[i] - wl[i]
        wcl[i] <- wc[i] - wl[i]
        yrc[i] <- yr[i] - yc[i]
        yrl[i] <- yr[i] - yl[i]
        ycl[i] <- yc[i] - yl[i]

        # Estimate parameters
        w[i] <- calc_w(
            wr[i], wc[i], wl[i],
            yr[i], yc[i], yl[i],
            wrc[i], wrl[i], wcl[i],
            yrc[i], yrl[i], ycl[i],
            xr[i]
        )
        lambda[i] <- calc_lambda(
            wr[i], wc[i], wl[i],
            yr[i], yc[i], yl[i],
            wrc[i], wrl[i], wcl[i],
            yrc[i], yrl[i], ycl[i],
            xr[i]
        )
        A[i] <- calc_A(
            wr[i], wc[i], wl[i],
            yr[i], yc[i], yl[i],
            wrc[i], wrl[i], wcl[i],
            yrc[i], yrl[i], ycl[i],
            lambda[i], w[i], xr[i]
        )

        # Calculate contribution of each lorentz curve to each data point
        cond <- (w[i] == 0) || (lambda[i] == 0) || (A[i] == 0)
        Z[, i] <- if (cond) 0 else abs(A[i] * (lambda[i] / (lambda[i]^2 + (xlmr - w[i])^2)))

        if (w[i] == 0 || lambda[i] == 0 || A[i] == 0) {
            Z[, i] <- 0
        } else {
            Z[, i] <- abs(A[i] * (lambda[i] / (lambda[i]^2 + (xlmr - w[i])^2)))
        }
    }

    # Print MSE
    mse <- mean((y[lmr] - rowSums(Z))^2)
    logf("MSE at peak tiplet positions: %.22f", mse)

    # Create return list
    P <- data.frame(il, ic, ir, rl, rc, rr, xl, xc, xr, yl, yc, yr, sl, sc, sr, ql, qc, qr)
    D <- data.frame(wl, wc, wr, wrc, wrl, wcl, yrc, yrl, ycl)
    named(A, lambda, w, Z, D, P) # nolint: object_usage_linter
}

# Taken from Appendix E of Koh et. al. 2009
calc_w <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    t1 <- wr^2 * yr * ycl
    t2 <- wl^2 * yl * yrc
    t3 <- wc^2 * yc * yrl
    t4 <- 2 * wrc * yr * yc
    t5 <- 2 * wcl * yc * yl
    t6 <- 2 * wrl * yr * yl
    w <- (t1 + t2 - t3) / (t4 + t5 - t6) + xr
    w[is.nan(w)] <- 0 # If (t4 + t5 - t6) is 0, then w is NaN. In this case we set w to 0.
    w
}

# Taken from Appendix E of Koh et. al. 2009
calc_lambda <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, xr) {
    lambda <- -((sqrt(abs((-wc^4 * yc^2 * yrl^2 - wr^4 * yr^2 * ycl^2 - wl^4 * yrc^2 * yl^2 + 4 * wc * wl^3 * yc * ((-yr) + yc) * yl^2 + 4 * wc^3 * wl * yc^2 * yl * ((-yr) + yl) + 4 * wr^3 * yr^2 * ycl * (wc * yc - wl * yl) + 4 * wr * yr * (wc^3 * yc^2 * yrl - wc * wl^2 * yc * (yr + yc - 2 * yl) * yl + wl^3 * yrc * yl^2 - wc^2 * wl * yc * yl * (yr - 2 * yc + yl)) + 2 * wc^2 * wl^2 * yc * yl * (yr^2 - 3 * yc * yl + yr * (yc + yl)) + 2 * wr^2 * yr * (-2 * wc * wl * yc * yl * (-2 * yr + yc + yl) + wl^2 * yl * (yr * (yc - 3 * yl) + yc * (yc + yl)) + wc^2 * yc * (yr * (-3 * yc + yl) + yl * (yc + yl)))))))) / (2 * sqrt((wr * yr * ycl + wl * yrc * yl + wc * yc * ((-yr) + yl))^2))
    lambda[is.nan(lambda)] <- 0
    lambda
}

# Taken from Appendix E of Koh et. al. 2009
calc_A <- function(wr, wc, wl, yr, yc, yl, wrc, wrl, wcl, yrc, yrl, ycl, lambda, w, xr) {
    A <- (-4 * wrc * wrl * wcl * yr * yc * yl * (wr * yr * ycl + wl * yl * yrc + wc * yc * (-yrl)) * lambda) / (wrc^4 * yr^2 * yc^2 - 2 * wrc^2 * yr * yc * (wrl^2 * yr + wcl^2 * yc) * yl + (wrl^2 * yr - wcl^2 * yc)^2 * yl^2)
    A[is.nan(A)] <- 0
    A
}

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# General Helpers #####
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~

#' @noRd
#' @title Calculate Lorentz Curve values
#'
#' @description
#' Calculates the values of a Lorentz Curve for a vector of input values `x`.
#' The Lorentz Curve is defined as \mjeqn{A \cdot \frac{\lambda}{\lambda^2 +
#' (x_i - x_0)^2}}.
#'
#' @param x Numeric vector of x values.
#' @param x0 Center of the Lorentz curve.
#' @param A Amplitude parameter of the Lorentz curve.
#' @param lambda Half width at half height of the Lorentz curve.
#'
#' @return Numeric vector of y values.
#'
#' @examples
#' x <- 1:10
#' x0 <- 5
#' A <- 10
#' lambda <- 2
#' y1 <- lorentz(x, x0, A, lambda)
#' y2 <- A * pi * dcauchy(x, location = x0, scale = lambda)
#' stopifnot(all.equal(y1, y2))
lorentz <- function(x, x0, A, lambda) {
    # For details see [Wikipedia > Cauchy_distribution > Properties_of_PDF]
    # (https://en.wikipedia.org/wiki/Cauchy_distribution#Properties_of_PDF)
    # in particular the formula below sentence "In physics, a three-parameter
    # Lorentzian function is often used".
    A * (lambda / (lambda^2 + (x - x0)^2))
}

lorentz_sup <- function(x, x0, A, lambda, lcpar = NULL) {
    if (is.list(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar$A
        if ("lambda" %in% nams) lambda <- lcpar$lambda
        if ("x_0" %in% nams) x0 <- lcpar$x_0
        if ("x0" %in% nams) x0 <- lcpar$x0
        if ("w" %in% nams) x0 <- lcpar$w
    }
    sapply(x, function(xi) {
        sum(abs(A * (lambda / (lambda^2 + (xi - x0)^2))))
    })
}

#' @noRd
#' @title Calculate Lorentz Curve Integrals
#' @description
#' Calculates the integral of a Lorentz curve for a vector of input values `x`.
lorentz_int <- function(x0, A, lambda, lcpar = NULL, limits = NULL) {
    if (is.list(lcpar)) {
        nams <- names(lcpar)
        if ("A" %in% nams) A <- lcpar$A
        if ("lambda" %in% nams) lambda <- lcpar$lambda
        if ("x_0" %in% nams) x0 <- lcpar$x_0
        if ("x0" %in% nams) x0 <- lcpar$x0
        if ("w" %in% nams) x0 <- lcpar$w
    }
    if (is.null(limits)) {
        A * pi
    } else {
        a <- min(limits)
        b <- max(limits)
        A * (atan((b - x0) / lambda) - atan((a - x0) / lambda))
    }
}

#' @noRd
#' @description
#' Before version 1.2 of 'metabodecon', the deconvolution functions
#' `generate_lorentz_curves()` and `MetaboDecon1D()` wrote their output
#' partially as txt files to their input folder. The txt files were named
#' "SPEC_NAME parameter.txt" and "SPEC_NAME approximated_spectrum.txt". Since
#' version 1.2 these txt files are no longer created by default, to prevent
#' accidental modifications of the input folders. However, to stay backwards
#' compatible, functions that used to read "SPEC_NAME parameter.txt" and
#' "SPEC_NAME approximated_spectrum.txt" still accept them as input (e.g.
#' `gen_feat_mat()`). I.e., in order to test this functionality, we still need a
#' way to create the corresponding txt files (which is no longer done by
#' `generate_lorentz_curves()`). That's the purpose of this function: it takes
#' the output of `generate_lorentz_curves()` as input and creates the (now
#' deprecated) "SPEC_NAME parameter.txt" and "SPEC_NAME
#' approximated_spectrum.txt" in folder `outdir`.
write_parameters_txt <- function(decon, outdir, verbose = FALSE) {
    if (is_decon_list(decon)) {
        for (obj in decon) write_parameters_txt(obj, outdir)
        return(invisible(NULL))
    }
    name <- decon$filename
    w_new <- decon$x_0
    lambda_new <- decon$lambda
    A_new <- decon$A
    noise_threshold <- rep(NA, length(A_new))
    pardf <- data.frame(rbind(w_new, lambda_new, A_new, noise_threshold))
    supdf <- data.frame(t(decon$spectrum_superposition))
    parfile <- file.path(outdir, paste(name, "parameters.txt"))
    supfile <- file.path(outdir, paste(name, "approximated_spectrum.txt"))
    if (verbose) cat(sprintf("Creating: %s\n", parfile))
    utils::write.table(pardf, parfile, sep = ",", col.names = FALSE, append = FALSE)
    if (verbose) cat(sprintf("Creating: %s\n", parfile))
    utils::write.table(supdf, supfile, sep = ",", col.names = FALSE, append = FALSE)
}

store_as_rds <- function(decons, make_rds, data_path) {
    if (is.character(make_rds)) {
        cat("Saving results as", make_rds, "\n")
        saveRDS(decons, make_rds)
    } else if (isTRUE(make_rds)) {
        rdspath <- file.path(data_path, "spectrum_data.rds")
        if (interactive()) {
            yes <- get_yn_input(sprintf("Save results as '%s'?", rdspath))
            if (yes) saveRDS(decons, rdspath)
        } else {
            logf("Skipping RDS save: confirmation required but not in interactive mode. For details see `help('generate_lorentz_curves')`.")
        }
    }
}
