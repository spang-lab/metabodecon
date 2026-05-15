
# API #####

#' @export
#'
#' @title Align Spectra
#'
#' @description
#' Align signals across  a  list  of  deconvoluted  spectra  using  the  'CluPA'
#' algorithm from the 'speaq' package, described  in  Beirnaert  et  al.  (2018)
#' <doi:10.1371/journal.pcbi.1006018>     and     Vu     et      al.      (2011)
#' <doi:10.1186/1471-2105-12-405> plus the additional peak combination described
#' in [metabodecon::combine_peaks()].
#'
#' @param x
#' An object of type `decons2` or `aligns`, as described in [metabodecon::metabodecon-classes].
#'
#' @param ref
#' Optional reference spectrum of type `align` or `decon2`. When supplied,
#' all spectra in `x` are aligned towards this reference. The reference is
#' prepended to `x` internally and removed from the result. If `NULL`
#' (default), the reference is chosen automatically.
#'
#' @param maxShift
#' Maximum number of datapoints a peak center may be shifted during CluPA
#' alignment. 50 is a suitable starting value for plasma spectra with a digital
#' resolution of 128K. Increase for urine or other sample types with larger
#' chemical-shift variation. Use `maxShift = 0L` to skip alignment entirely:
#' each peak's aligned center `x0al` is set equal to its fitted center `x0`.
#'
#' @param verbose
#' Whether to print progress messages during alignment.
#'
#' @param nworkers
#' Number of parallel workers. Default is 1 (no parallelism).
#'
#' @return
#' An object of type `aligns` as described in [metabodecon::metabodecon-classes].
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#' aligned <- align(decons)
align <- function(x, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1) {
    stopifnot(
        inherits(x, "decons2"),
        is_int(maxShift,1), is_bool(verbose,1), is_int(nworkers,1),
        is.null(ref) || inherits(ref, "decon2")
    )
    clupa(x, ref, maxShift, verbose, nworkers)
}

#' @export
#' @name alignment_funs
#' @rdname alignment_funs
#'
#' @title Alignment functions for fit_mdm
#'
#' @description
#' Pluggable alignment backends accepted by the `align_fun` argument of
#' [metabodecon::fit_mdm()] / [metabodecon::benchmark()].
#'
#' - [metabodecon::clupa()]: hierarchical-clustering peak alignment
#'   (recursive FFT shifts on spectrum sub-segments,
#'   Beirnaert et al. 2018, Vu et al. 2011). Default of
#'   [metabodecon::align()].
#' - [metabodecon::vopa()]: vote-based peak alignment. Estimates one
#'   global integer DP shift per spectrum from a weighted vote over
#'   pairwise peak matches, then refines per-peak shifts via
#'   interpolation.
#' - [metabodecon::glopa()]: global peak alignment. Per spectrum, find
#'   the single integer DP shift in `[-maxShift, maxShift]` that
#'   maximises FFT cross-correlation of the raw signal intensities
#'   against the reference, then apply it to all peak centers.
#' - [metabodecon::identity_align()]: no-op. Returns its first argument
#'   unchanged.
#'
#' @param x A `decons2` or `aligns` object.
#' @param ref Optional reference spectrum (`align` or `decon2`). When
#'   `NULL`, chosen automatically.
#' @param maxShift Maximum number of datapoints a peak center may be
#'   shifted.
#' @param verbose Print progress messages?
#' @param nworkers Number of parallel workers.
#' @param full If `TRUE` also recompute the aligned superposition.
#' @param use_speaq Use `speaq::dohCluster` instead of the bundled
#'   implementation. Only used by `clupa`.
#' @param ... Ignored (`identity_align` only).
#' @return An object of class `aligns`.
clupa <- function(
    x, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1,
    full=TRUE, use_speaq=FALSE
) {
    ndps <- vapply(x, function(s) length(s$cs), integer(1))
    if (length(unique(ndps)) > 1) stop("All spectra must have the same number of data points.")
    ref <- ref %||% find_ref(x)
    aligns <- mcmapply(
        nworkers, align_decon, x,
        MoreArgs=list(ref, maxShift, full=full, use_speaq=use_speaq, method="clupa")
    )
    class(aligns) <- c("aligns", "decons2", "spectra")
    aligns
}

#' @export
#' @rdname alignment_funs
vopa <- function(
    x, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1, full=TRUE
) {
    ndps <- vapply(x, function(s) length(s$cs), integer(1))
    if (length(unique(ndps)) > 1) stop("All spectra must have the same number of data points.")
    ref <- ref %||% find_ref(x)
    aligns <- mcmapply(nworkers, align_decon, x,
        MoreArgs=list(ref, maxShift, full=full, use_speaq=FALSE, method="shift"))
    class(aligns) <- c("aligns", "decons2", "spectra")
    aligns
}

#' @export
#' @rdname alignment_funs
glopa <- function(
    x, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1, full=TRUE
) {
    ndps <- vapply(x, function(s) length(s$cs), integer(1))
    if (length(unique(ndps)) > 1) stop("All spectra must have the same number of data points.")
    ref <- ref %||% find_ref(x)
    aligns <- mcmapply(nworkers, glopa_one, x,
        MoreArgs=list(ref=ref, maxShift=maxShift, full=full))
    class(aligns) <- c("aligns", "decons2", "spectra")
    aligns
}

glopa_one <- function(x, ref, maxShift, full=TRUE) {
    if (maxShift == 0L) {
        x$lcpar$x0al <- x$lcpar$x0
        x$lcpar$pcial <- round(convert_pos(x$lcpar$x0, x$cs, seq_along(x$cs)))
        if (full) x$sit$supal <- lorentz_sup(x$cs, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
        class(x) <- c("align", "decon2", "spectrum")
        return(x)
    }
    sig_ref <- ref$si %||% ref$sit$sup
    sig_tar <- x$si %||% x$sit$sup
    adj <- fft_shift(sig_ref, sig_tar, maxShift=maxShift)
    pci_x <- round(convert_pos(x$lcpar$x0, x$cs, seq_along(x$cs)))
    pcial <- pmin(length(x$cs), pmax(1L, pci_x - adj$stepAdj))
    x$lcpar$x0al <- x$cs[pcial]
    x$lcpar$pcial <- pcial
    if (full) x$sit$supal <- lorentz_sup(x$cs, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

#' @export
#' @rdname alignment_funs
identity_align <- function(x, ...) x

# Internal #####

align_decon <- function(x, ref, maxShift, full=TRUE, use_speaq=FALSE, method="clupa") {
    if (maxShift == 0L) {
        x$lcpar$x0al <- x$lcpar$x0
        x$lcpar$pcial <- round(convert_pos(x$lcpar$x0, x$cs, seq_along(x$cs)))
        if (full) x$sit$supal <- lorentz_sup(x$cs, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
        class(x) <- c("align", "decon2", "spectrum")
        return(x)
    }
    pci_x <- round(convert_pos(x$lcpar$x0, x$cs, seq_along(x$cs)))
    pci_ref <- round(convert_pos(ref$lcpar$x0, ref$cs, seq_along(ref$cs)))

    if (method == "clupa") {
        np_x <- length(pci_x); np_ref <- length(pci_ref)
        obj <- hclust_align(
            refSpec=ref$sit$sup, tarSpec=x$sit$sup,
            peakList=c(pci_ref, pci_x),
            peakLabel=c(rep(1, np_ref), rep(0, np_x)),
            startP=1, endP=length(x$sit$sup),
            maxShift=maxShift, use_speaq=use_speaq
        )
        if (length(obj$peakList) != np_ref + np_x) stop("Lost peaks during alignment")
        pcial <- obj$peakList[(np_ref + 1):(np_ref + np_x)]
    } else {
        lam_x <- abs(round(convert_width(x$lcpar$lambda, x$cs, seq_along(x$cs))))
        lam_ref <- abs(round(convert_width(ref$lcpar$lambda, ref$cs, seq_along(ref$cs))))
        ref_pd <- list(x0=pci_ref, lambda=pmax(lam_ref, 1),
                       A=abs(ref$lcpar$A), n=length(ref$cs))
        tar_pd <- list(x0=pci_x, lambda=pmax(lam_x, 1),
                       A=abs(x$lcpar$A), n=length(x$cs))
        shift0 <- estimate_peak_shift(ref_pd, tar_pd, maxShift)
        m <- match_peak_pairs(ref_pd, tar_pd, shift=shift0)
        if (length(m$idx_tar) < 3) {
            pcial <- pmin(length(x$cs), pmax(1L, round(pci_x - shift0)))
        } else {
            shift <- interp_shift(x=pci_x[m$idx_tar], shift=m$delta,
                weight=m$weight, xout=pci_x, maxShift=maxShift, fallback=shift0)
            pcial <- pmin(length(x$cs), pmax(1L, round(pci_x - shift)))
        }
    }

    x$lcpar$x0al <- x$cs[pcial]
    x$lcpar$pcial <- pcial
    if (full) x$sit$supal <- lorentz_sup(x$cs, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

find_ref <- function(x) {
    pci <- lapply(x, function(s) round(convert_pos(s$lcpar$x0, s$cs, seq_along(s$cs))))
    x[[find_ref_ind(pci)$refInd]]
}

# Shift alignment helpers #####

weighted_mean <- function(x, w) sum(x * w) / sum(w)

weighted_median <- function(x, w) {
    ord <- order(x); x <- x[ord]; w <- w[ord]
    x[which(cumsum(w) >= sum(w) / 2)[1]]
}

# Estimates a global integer DP shift from ref to tar by a weighted vote
# over all pairwise peak matches within maxShift.
estimate_peak_shift <- function(ref_pd, tar_pd, maxShift=50) {
    rx <- ref_pd$x0; tx <- tar_pd$x0
    rA <- pmax(abs(ref_pd$A), .Machine$double.eps)
    tA <- pmax(abs(tar_pd$A), .Machine$double.eps)
    rl <- pmax(ref_pd$lambda, 1); tl <- pmax(tar_pd$lambda, 1)
    scores <- numeric(2 * maxShift + 1)
    bins <- seq.int(-maxShift, maxShift)
    for (i in seq_along(tx)) {
        j1 <- max(1L, findInterval(tx[i] - maxShift, rx) + 1L)
        j2 <- min(length(rx), findInterval(tx[i] + maxShift, rx))
        if (j1 > j2) next
        idx <- j1:j2
        delta <- tx[i] - rx[idx]
        sim <- pmin(tl[i], rl[idx]) / pmax(tl[i], rl[idx])
        w <- pmin(tA[i], rA[idx]) * sim
        pos <- round(delta) + maxShift + 1L
        for (k in seq_along(pos)) scores[pos[k]] <- scores[pos[k]] + w[k]
    }
    if (all(scores == 0)) return(0)
    bins[which.max(scores)]
}

# Matches each tar peak to its nearest ref peak (after applying shift0).
# Returns idx_tar, delta (tar - ref) and weight for each matched pair.
match_peak_pairs <- function(ref_pd, tar_pd, shift=0, tol_mult=2) {
    rx <- ref_pd$x0; tx <- tar_pd$x0; xadj <- tx - shift
    idx <- findInterval(xadj, rx)
    lo <- pmax(idx, 1L); hi <- pmin(idx + 1L, length(rx))
    idx_ref <- lo; idx_ref[abs(xadj - rx[hi]) < abs(xadj - rx[lo])] <- hi[abs(xadj - rx[hi]) < abs(xadj - rx[lo])]
    tl <- pmax(tar_pd$lambda, 1); rl <- pmax(ref_pd$lambda[idx_ref], 1)
    keep <- abs(xadj - rx[idx_ref]) <= pmax(1, tol_mult * (tl + rl))
    if (!any(keep)) return(list(idx_tar=integer(0), delta=numeric(0), weight=numeric(0)))
    idx_tar <- which(keep); idx_ref <- idx_ref[keep]; rl <- rl[keep]
    delta <- tx[idx_tar] - rx[idx_ref]
    sim <- pmin(tl[keep], rl) / pmax(tl[keep], rl)
    w <- pmin(abs(tar_pd$A[idx_tar]), abs(ref_pd$A[idx_ref])) * sim
    ord <- order(abs(delta), -w)
    idx_tar <- idx_tar[ord]; idx_ref <- idx_ref[ord]
    delta <- delta[ord]; w <- w[ord]
    ok <- !duplicated(idx_ref)
    list(idx_tar=idx_tar[ok], delta=delta[ok],
         weight=pmax(w[ok], .Machine$double.eps))
}

# Interpolates per-peak shifts at xout positions from sparse matched-peak
# (x, shift, weight) triples. Groups peaks into up to 12 knots, fits a
# piecewise-linear interpolant, and clamps to [-maxShift, maxShift].
interp_shift <- function(x, shift, weight, xout, maxShift=50, fallback=0) {
    n <- length(x)
    if (n == 0) return(rep(fallback, length(xout)))
    if (n == 1) return(rep(shift[1], length(xout)))
    ord <- order(x); x <- x[ord]; shift <- shift[ord]; weight <- weight[ord]
    ng <- min(12L, max(1L, floor(n / 4L)))
    if (ng == 1L) return(rep(weighted_median(shift, weight), length(xout)))
    grp <- cut(seq_len(n), breaks=ng, labels=FALSE)
    knot_x <- vapply(split(seq_len(n), grp),
                     function(i) weighted_mean(x[i], weight[i]), numeric(1))
    knot_s <- vapply(split(seq_len(n), grp),
                     function(i) weighted_median(shift[i], weight[i]), numeric(1))
    if (length(knot_x) == 1) return(rep(knot_s[1], length(xout)))
    y <- stats::approx(knot_x, knot_s, xout=xout, rule=2)$y
    y <- pmin(maxShift, pmax(-maxShift, y))
    y[!is.finite(y)] <- fallback
    y
}

# Speaq #####

# Lightweight replacements for the speaq functions used by metabodecon. Only the
# subset of functionality actually needed is implemented. Original speaq
# package: Beirnaert et al. (2018) <doi:10.1371/journal.pcbi.1006018> and Vu et
# al. (2011) <doi:10.1186/1471-2105-12-405>.

#' @noRd
#'
#' @description
#' Find the reference spectrum from a list of peak indices. The
#' reference is the spectrum whose peaks have the smallest total
#' distance to all other spectra's peaks.
#'
#' Replacement for `speaq::findRef()`.
#'
#' @param peakList
#' A list of integer vectors, each containing peak indices for
#' one spectrum.
#'
#' @return
#' A list with elements `refInd` (index of the best reference)
#' and `orderSpec` (all indices ordered by suitability).
#'
#' @author 2025 Tobias Schmidt: initial version.
find_ref_ind <- function(peakList) {
    n <- length(peakList)
    sumDis <- double(n)
    for (r in seq_len(n)) {
        rp <- sort(peakList[[r]])
        for (t in seq_len(n)) {
            if (r == t) next
            tp <- peakList[[t]]
            # For each target peak, find nearest ref peak
            # using binary search on the sorted ref peaks.
            idx <- findInterval(tp, rp)
            lo <- pmax(idx, 1L)
            hi <- pmin(idx + 1L, length(rp))
            d <- pmin(abs(tp - rp[lo]), abs(tp - rp[hi]))
            sumDis[r] <- sumDis[r] + sum(d)
        }
    }
    ord <- order(sumDis)
    list(refInd = ord[1], orderSpec = ord)
}

#' @noRd
#'
#' @description
#' Compute the optimal integer shift between a reference and
#' target spectrum segment using FFT cross-correlation.
#'
#' Replacement for `speaq::findShiftStepFFT()`.
#'
#' @param refSpec Numeric vector (reference segment).
#' @param tarSpec Numeric vector (target segment, same length).
#' @param maxShift Maximum allowed shift in either direction.
#'
#' @return
#' A list with `stepAdj` (integer shift) and `corValue`
#' (cross-correlation value at that shift).
#'
#' @author 2025 Tobias Schmidt: initial version.
fft_shift <- function(refSpec, tarSpec, maxShift) {
    M <- length(refSpec)
    pad <- 2^ceiling(log2(M)) - M
    r <- c(refSpec * 1e6, double(pad))
    s <- c(tarSpec * 1e6, double(pad))
    N <- M + pad
    R <- stats::fft(r) * Conj(stats::fft(s)) / N
    vals <- Re(stats::fft(R, inverse = TRUE)) / N
    if (maxShift == 0 || maxShift > M) maxShift <- M
    if (anyNA(vals)) return(list(corValue = -1, stepAdj = 0L))
    nv <- length(vals)
    # Interleave forward/backward indices to preserve the same
    # tie-breaking order as the original speaq loop:
    # lag 0, -1, 1, -2, 2, -3, ...
    fwd <- seq_len(maxShift)
    bwd <- seq.int(nv, nv - maxShift + 1L)
    idx <- as.vector(rbind(fwd, bwd))
    best <- which.max(vals[idx])
    maxpos <- idx[best]
    maxi <- vals[maxpos]
    if (maxi < 0.1) return(list(corValue = maxi, stepAdj = 0L))
    lag <- if (maxpos > nv / 2) maxpos - nv - 1L else maxpos - 1L
    list(corValue = maxi, stepAdj = lag)
}

#' @noRd
#'
#' @description
#' Shift a spectrum segment by `step` positions, padding the
#' vacated side with the nearest edge value.
#'
#' Replacement for `speaq::doShift()`.
#'
#' @param seg Numeric vector.
#' @param step Integer shift (positive = shift right, negative = shift left).
#'
#' @return Shifted numeric vector of the same length.
#'
#' @author 2025 Tobias Schmidt: initial version.
do_shift <- function(seg, step) {
    n <- length(seg)
    out <- double(n)
    # Copy shifted values
    src <- seq_len(n)
    dst <- src + step
    valid <- dst >= 1L & dst <= n
    out[dst[valid]] <- seg[src[valid]]
    # Pad edges with nearest boundary value
    if (step > 0) {
        out[seq_len(step)] <- out[step + 1L]
    } else if (step < 0) {
        start <- n + step
        out[start:n] <- out[start - 1L]
    } else {
        # step == 0: replicate speaq quirk where last element
        # gets overwritten by second-to-last
        out[n] <- out[n - 1L]
    }
    out
}

#' @noRd
#'
#' @description
#' Align a target spectrum to a reference spectrum using
#' recursive hierarchical-clustering-based segmentation with
#' FFT cross-correlation shifts.
#'
#' Replacement for `speaq::hClustAlign()` with
#' `acceptLostPeak = FALSE` and `distanceMethod = "average"`.
#'
#' @param refSpec Numeric vector (full reference spectrum).
#' @param tarSpec Numeric vector (full target spectrum).
#' @param peakList Integer vector of peak positions (ref then target, interleaved via labels).
#' @param peakLabel Integer vector, 1 for ref peaks, 0 for target peaks.
#' @param startP Start index of the segment to align.
#' @param endP End index of the segment to align.
#' @param maxShift Maximum shift per recursion level.
#'
#' @return
#' A list with `tarSpec` (aligned target spectrum) and
#' `peakList` (updated peak positions).
#'
#' @author 2025 Tobias Schmidt: initial version.
hclust_align <- function(
    refSpec, tarSpec, peakList, peakLabel, startP, endP, maxShift, use_speaq = FALSE
) {

    if (use_speaq) return(
        speaq::hClustAlign(
            refSpec=refSpec, tarSpec=tarSpec, peakList=peakList,
            peakLabel=peakLabel, startP=startP, endP=endP,
            distanceMethod="average", maxShift=maxShift, acceptLostPeak=FALSE
        )
    )

    minPk <- min(peakList)
    maxPk <- max(peakList)

    # Narrow the active region to signal boundaries
    startCheckP <- startP + which.min(tarSpec[startP:(minPk - 1L)]) - 1L
    if (is.na(startCheckP) || startCheckP < 1L) startCheckP <- startP
    endCheckP <- maxPk + which.min(tarSpec[(maxPk + 1L):endP])
    if (is.na(endCheckP) || endCheckP > length(tarSpec)) endCheckP <- endP

    if ((endCheckP - startCheckP) < 2L) return(list(tarSpec = tarSpec, peakList = peakList))

    # FFT cross-correlation to find the best shift
    adj <- fft_shift(
        refSpec[startCheckP:endCheckP],
        tarSpec[startCheckP:endCheckP],
        maxShift = maxShift
    )

    if (adj$stepAdj != 0L) {
        # acceptLostPeak = FALSE: only shift if no peaks are
        # pushed outside the region
        ok <- (adj$stepAdj < 0 && adj$stepAdj + minPk >= startCheckP) ||
              (adj$stepAdj > 0 && adj$stepAdj + maxPk <= endCheckP)
        if (ok) {
            seg <- tarSpec[startCheckP:endCheckP]
            tarSpec[startCheckP:endCheckP] <- do_shift(seg, adj$stepAdj)
            tar_idx <- which(peakLabel == 0L)
            peakList[tar_idx] <- peakList[tar_idx] + adj$stepAdj
            lost <- which(peakList <= 0L | peakList > length(tarSpec))
            if (length(lost) > 0L) {
                peakList <- peakList[-lost]
                peakLabel <- peakLabel[-lost]
            }
        }
    }

    # Need >= 3 peaks total (ref + target) for splitting
    if (length(peakList) < 3L) {
        return(list(tarSpec = tarSpec, peakList = peakList))
    }

    # Split peaks into 2 groups using average-linkage
    # hierarchical clustering, matching speaq::hClustAlign.
    hc <- stats::hclust(stats::dist(peakList), method = "average")
    cl <- stats::cutree(hc, h = hc$height[length(hc$height) - 1])
    if (length(unique(cl)) < 2L) {
        return(list(tarSpec = tarSpec, peakList = peakList))
    }
    left_set <- which(cl == 1)
    right_set <- which(cl == 2)

    sub1 <- peakList[left_set]
    lab1 <- peakLabel[left_set]
    id1 <- left_set
    sub2 <- peakList[right_set]
    lab2 <- peakLabel[right_set]
    id2 <- right_set

    max1 <- max(sub1)
    min2 <- min(sub2)

    # speaq handles both orderings (cluster 1 left or right)
    if (max1 < min2) {
        endP1 <- max1 + which.min(tarSpec[(max1 + 1L):(min2 - 1L)])
        if (is.na(endP1) || endP1 > length(tarSpec)) endP1 <- max1
        startP2 <- endP1 + 1L
    } else {
        # Cluster 1 is right, cluster 2 is left — swap
        tmp_set <- left_set; left_set <- right_set
        right_set <- tmp_set
        sub1 <- peakList[left_set]
        lab1 <- peakLabel[left_set]
        id1 <- left_set
        sub2 <- peakList[right_set]
        lab2 <- peakLabel[right_set]
        id2 <- right_set
        max1 <- max(sub1); min2 <- min(sub2)
        endP1 <- max1 + which.min(tarSpec[(max1 + 1L):(min2 - 1L)])
        if (is.na(endP1) || endP1 > length(tarSpec)) endP1 <- max1
        startP2 <- endP1 + 1L
    }
    if (length(unique(lab1)) > 1L) {
        res <- hclust_align(refSpec, tarSpec, sub1, lab1, startP, endP1, maxShift)
        tarSpec <- res$tarSpec
        peakList[id1] <- pad_peaks(res$peakList, length(id1))
    }
    if (length(unique(lab2)) > 1L) {
        res <- hclust_align(refSpec, tarSpec, sub2, lab2, startP2, endP, maxShift)
        tarSpec <- res$tarSpec
        peakList[id2] <- pad_peaks(res$peakList, length(id2))
    }
    list(tarSpec = tarSpec, peakList = peakList)
}

#' @noRd
#' @description
#' Pad or truncate a peak vector to length `n`, replicating
#' the first element when peaks were lost during alignment.
#' Matches the speaq convention.
pad_peaks <- function(peaks, n) {
    if (length(peaks) >= n) {
        peaks[seq_len(n)]
    } else {
        c(peaks, rep(peaks[1L], n - length(peaks)))
    }
}
