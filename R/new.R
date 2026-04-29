# Experimental alignment helpers.
# These functions implement an alternative peak-shift based alignment algorithm
# (method = 3 / align_fast). They are kept here for research purposes and are
# not part of the public API.

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
align_fast <- function(peakData,
                       refInd = 1,
                       maxShift = 50,
                       verbose = TRUE,
                       nworkers = 1) {
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    startTime <- proc.time()
    nspec <- length(peakData)
    logf("Running align_fast with maxShift = %d on %d spectra", maxShift, nspec)
    new_peakList <- lapply(peakData, function(x) x$x0)
    tarInd <- setdiff(seq_len(nspec), refInd)
    if (length(tarInd) > 0) {
        nw <- min(nworkers, length(tarInd))
        est_shift <- estimate_peak_shift
        match_pairs <- match_peak_pairs
        w_mean <- weighted_mean
        w_med <- weighted_median
        interp_shift <- function(x,
                                 shift,
                                 weight,
                                 xout,
                                 maxShift = 50,
                                 fallback = 0) {
            n <- length(x)
            if (n == 0) return(rep(fallback, length(xout)))
            if (n == 1) return(rep(shift[1], length(xout)))
            ord <- order(x)
            x <- x[ord]
            shift <- shift[ord]
            weight <- weight[ord]
            ng <- min(12L, max(1L, floor(n / 4L)))
            if (ng == 1L) return(rep(w_med(shift, weight), length(xout)))
            grp <- cut(seq_len(n), breaks = ng, labels = FALSE)
            knot_x <- vapply(split(seq_len(n), grp), function(i) {w_mean(x[i], weight[i])}, numeric(1))
            knot_s <- vapply(split(seq_len(n), grp), function(i) {w_med(shift[i], weight[i])}, numeric(1))
            if (length(knot_x) == 1) return(rep(knot_s[1], length(xout)))
            y <- stats::approx(knot_x, knot_s, xout = xout, rule = 2)$y
            y <- pmin(maxShift, pmax(-maxShift, y))
            y[!is.finite(y)] <- fallback
            y
        }
        worker_fun <- function(tarData, refData, maxShift) {
            shift0 <- est_shift(refData, tarData, maxShift)
            match <- match_pairs(refData, tarData, shift = shift0)
            if (length(match$idx_tar) < 3) {
                out <- tarData$x0 - shift0
                out <- pmin(tarData$n, pmax(1, out))
                return(out)
            }
            shift <- interp_shift(
                x = tarData$x0[match$idx_tar],
                shift = match$delta,
                weight = match$weight,
                xout = tarData$x0,
                maxShift = maxShift,
                fallback = shift0
            )
            out <- tarData$x0 - shift
            out <- pmin(tarData$n, pmax(1, out))
            out
        }
        res <- mcmapply(
            nw,
            worker_fun,
            peakData[tarInd],
            MoreArgs = list(refData = peakData[[refInd]], maxShift = maxShift),
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE,
            loadpkg = FALSE
        )
        new_peakList[tarInd] <- res
    }
    elapsed <- (proc.time() - startTime)[3]
    logf("Finished align_fast in %.1f s", elapsed)
    list(Y = NULL, new_peakList = new_peakList)
}

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
weighted_mean <- function(x, w) {
    sum(x * w) / sum(w)
}

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
weighted_median <- function(x, w) {
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    i <- which(cumsum(w) >= sum(w) / 2)[1]
    x[i]
}

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
get_peak_align_data <- function(decon2) {
    cs <- decon2$cs
    dp <- seq_along(cs)
    list(
        x0 = convert_pos(decon2$lcpar$x0, cs, dp),
        lambda = abs(convert_width(decon2$lcpar$lambda, cs, dp)),
        A = abs(decon2$lcpar$A),
        n = length(cs)
    )
}

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
match_peak_pairs <- function(refData, tarData, shift = 0, tol_mult = 2) {
    rx <- refData$x0
    tx <- tarData$x0
    xadj <- tx - shift
    idx <- findInterval(xadj, rx)
    lo <- pmax(idx, 1L)
    hi <- pmin(idx + 1L, length(rx))
    dlo <- abs(xadj - rx[lo])
    dhi <- abs(xadj - rx[hi])
    idx_ref <- lo
    use_hi <- dhi < dlo
    idx_ref[use_hi] <- hi[use_hi]
    resid <- xadj - rx[idx_ref]
    tl <- pmax(tarData$lambda, 1)
    rl <- pmax(refData$lambda[idx_ref], 1)
    tol <- pmax(1, tol_mult * (tl + rl))
    keep <- abs(resid) <= tol
    if (!any(keep)) {
        return(list(idx_tar = integer(0), delta = numeric(0), weight = numeric(0)))
    }
    idx_tar <- which(keep)
    idx_ref <- idx_ref[keep]
    rl <- rl[keep]
    delta <- tx[idx_tar] - rx[idx_ref]
    sim <- pmin(tl[keep], rl) / pmax(tl[keep], rl)
    w <- pmin(abs(tarData$A[keep]), abs(refData$A[idx_ref])) * sim
    ord <- order(abs(resid[keep]), -w)
    idx_tar <- idx_tar[ord]
    idx_ref <- idx_ref[ord]
    delta <- delta[ord]
    w <- w[ord]
    ok <- !duplicated(idx_ref)
    list(
        idx_tar = idx_tar[ok],
        delta = delta[ok],
        weight = pmax(w[ok], .Machine$double.eps)
    )
}

#' @noRd
#' @author 2026 Tobias Schmidt: initial version.
estimate_peak_shift <- function(refData, tarData, maxShift = 50) {
    rx <- refData$x0
    tx <- tarData$x0
    rA <- pmax(abs(refData$A), .Machine$double.eps)
    tA <- pmax(abs(tarData$A), .Machine$double.eps)
    rl <- pmax(refData$lambda, 1)
    tl <- pmax(tarData$lambda, 1)
    scores <- numeric(2 * maxShift + 1)
    bins <- seq.int(-maxShift, maxShift)
    for (i in seq_along(tx)) {
        j1 <- findInterval(tx[i] - maxShift, rx) + 1L
        j2 <- findInterval(tx[i] + maxShift, rx)
        j1 <- max(1L, j1)
        j2 <- min(length(rx), j2)
        if (j1 > j2) next
        idx <- j1:j2
        delta <- tx[i] - rx[idx]
        sim <- pmin(tl[i], rl[idx]) / pmax(tl[i], rl[idx])
        w <- pmin(tA[i], rA[idx]) * sim
        pos <- round(delta) + maxShift + 1L
        for (k in seq_along(pos)) {
            scores[pos[k]] <- scores[pos[k]] + w[k]
        }
    }
    if (all(scores == 0)) return(0)
    bins[which.max(scores)]
}
