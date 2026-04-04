# Lightweight replacements for the speaq functions used by metabodecon.
# Only the subset of functionality actually needed is implemented.
# Original speaq package: Beirnaert et al. (2018) <doi:10.1371/journal.pcbi.1006018>
# and Vu et al. (2011) <doi:10.1186/1471-2105-12-405>.

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
find_ref <- function(peakList) {
    n <- length(peakList)
    sumDis <- double(n)
    for (r in seq_len(n)) {
        for (t in seq_len(n)) {
            if (r == t) next
            for (pk in peakList[[t]]) {
                sumDis[r] <- sumDis[r] + min(abs(pk - peakList[[r]]))
            }
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
    maxi <- -1
    maxpos <- 1L
    nv <- length(vals)
    for (i in seq_len(maxShift)) {
        if (vals[i] > maxi) {
            maxi <- vals[i]
            maxpos <- i
        }
        j <- nv - i + 1L
        if (vals[j] > maxi) {
            maxi <- vals[j]
            maxpos <- j
        }
    }
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
#' @param step Integer shift (positive = shift right, negative
#'   = shift left).
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
#' @param peakList Integer vector of peak positions (ref then
#'   target, interleaved via labels).
#' @param peakLabel Integer vector, 1 for ref peaks, 0 for
#'   target peaks.
#' @param startP Start index of the segment to align.
#' @param endP End index of the segment to align.
#' @param maxShift Maximum shift per recursion level.
#'
#' @return
#' A list with `tarSpec` (aligned target spectrum) and
#' `peakList` (updated peak positions).
#'
#' @author 2025 Tobias Schmidt: initial version.
hclust_align <- function(refSpec, tarSpec, peakList, peakLabel,
                         startP, endP, maxShift) {

    minPk <- min(peakList)
    maxPk <- max(peakList)

    # Narrow the active region to signal boundaries
    startCheckP <- startP +
        which.min(tarSpec[startP:(minPk - 1L)]) - 1L
    if (is.na(startCheckP) || startCheckP < 1L) {
        startCheckP <- startP
    }
    endCheckP <- maxPk +
        which.min(tarSpec[(maxPk + 1L):endP])
    if (is.na(endCheckP) || endCheckP > length(tarSpec)) {
        endCheckP <- endP
    }

    if ((endCheckP - startCheckP) < 2L) {
        return(list(tarSpec = tarSpec, peakList = peakList))
    }

    # FFT cross-correlation to find the best shift
    adj <- fft_shift(
        refSpec[startCheckP:endCheckP],
        tarSpec[startCheckP:endCheckP],
        maxShift = maxShift
    )

    if (adj$stepAdj != 0L) {
        # acceptLostPeak = FALSE: only shift if no peaks are
        # pushed outside the region
        ok <- (adj$stepAdj < 0 &&
                   adj$stepAdj + minPk >= startCheckP) ||
              (adj$stepAdj > 0 &&
                   adj$stepAdj + maxPk <= endCheckP)
        if (ok) {
            seg <- tarSpec[startCheckP:endCheckP]
            tarSpec[startCheckP:endCheckP] <- do_shift(
                seg, adj$stepAdj
            )
            tar_idx <- which(peakLabel == 0L)
            peakList[tar_idx] <- peakList[tar_idx] + adj$stepAdj
            lost <- which(peakList <= 0L |
                              peakList > length(tarSpec))
            if (length(lost) > 0L) {
                peakList <- peakList[-lost]
                peakLabel <- peakLabel[-lost]
            }
        }
    }

    # Need >= 3 peaks total (ref + target) for clustering
    if (length(peakList) < 3L) {
        return(list(tarSpec = tarSpec, peakList = peakList))
    }

    # Hierarchical clustering to split into 2 segments
    hc <- stats::hclust(stats::dist(peakList), method = "average")
    cl <- stats::cutree(hc, h = hc$height[length(hc$height) - 1L])
    if (length(unique(cl)) < 2L) {
        return(list(tarSpec = tarSpec, peakList = peakList))
    }

    id1 <- which(cl == 1L)
    sub1 <- peakList[id1]
    lab1 <- peakLabel[id1]
    id2 <- which(cl == 2L)
    sub2 <- peakList[id2]
    lab2 <- peakLabel[id2]

    max1 <- max(sub1)
    min2 <- min(sub2)

    if (max1 < min2) {
        # Cluster 1 is left, cluster 2 is right
        endP1 <- max1 + which.min(tarSpec[(max1 + 1L):(min2 - 1L)])
        if (is.na(endP1) || endP1 > length(tarSpec)) endP1 <- max1
        startP2 <- endP1 + 1L
        if (length(unique(lab1)) > 1L) {
            res <- hclust_align(refSpec, tarSpec, sub1, lab1,
                startP, endP1, maxShift)
            tarSpec <- res$tarSpec
            peakList[id1] <- pad_peaks(res$peakList, length(id1))
        }
        if (length(unique(lab2)) > 1L) {
            res <- hclust_align(refSpec, tarSpec, sub2, lab2,
                startP2, endP, maxShift)
            tarSpec <- res$tarSpec
            peakList[id2] <- pad_peaks(res$peakList, length(id2))
        }
    } else {
        # Cluster 2 is left, cluster 1 is right
        max2 <- max(sub2)
        min1 <- min(sub1)
        endP2 <- max2 + which.min(tarSpec[(max2 + 1L):(min1 - 1L)])
        if (is.na(endP2) || endP2 > length(tarSpec)) endP2 <- max2
        startP1 <- endP2 + 1L
        if (length(unique(lab2)) > 1L) {
            res <- hclust_align(refSpec, tarSpec, sub2, lab2,
                startP, endP2, maxShift)
            tarSpec <- res$tarSpec
            peakList[id2] <- pad_peaks(res$peakList, length(id2))
        }
        if (length(unique(lab1)) > 1L) {
            res <- hclust_align(refSpec, tarSpec, sub1, lab1,
                startP1, endP, maxShift)
            tarSpec <- res$tarSpec
            peakList[id1] <- pad_peaks(res$peakList, length(id1))
        }
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
