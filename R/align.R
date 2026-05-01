
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
#' An object of type `decons2` or `aligns`, as described in [metabodecon-classes].
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
#' chemical-shift variation.
#'
#' @param verbose
#' Whether to print progress messages during alignment.
#'
#' @param nworkers
#' Number of parallel workers. Default is 1 (no parallelism).
#'
#' @return
#' An object of type `aligns` as described in [metabodecon-classes].
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
    align_decons(x, ref, maxShift, verbose, nworkers)
}

# Internal #####

align_decons <- function(
    x, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1,
    full=TRUE, use_speaq=FALSE
) {

    # Validate consistent number of data points across spectra
    ndps <- vapply(x, function(s) length(s$cs), integer(1))
    if (length(unique(ndps)) > 1) {
        stop("All spectra must have the same number of data points.")
    }

    # Do alignments
    ref <- ref %||% find_ref(x)
    aligns <- mcmapply(
        nworkers, align_decon, x,
        MoreArgs = list(ref, maxShift, full=full, use_speaq=use_speaq)
    )
    class(aligns) <- c("aligns", "decons2", "spectra")
    aligns
}

align_decon <- function(x, ref, maxShift, full=TRUE, use_speaq=FALSE) {

    # Init Helpers
    pci_x <- round(convert_pos(x$lcpar$x0, x$cs, seq_along(x$cs)))
    pci_ref <- round(convert_pos(ref$lcpar$x0, ref$cs, seq_along(ref$cs)))
    np_x <- length(pci_x)
    np_ref <- length(pci_ref)
    np_tot <- np_ref + np_x

    # Do 'Cluster Based Peak Alignment' (CluPA)
    obj <- hclust_align(
        refSpec = ref$sit$sup,
        tarSpec = x$sit$sup,
        peakList = c(pci_ref, pci_x),
        peakLabel = c(rep(1, np_ref), rep(0, np_x)),
        startP = 1,
        endP = length(x$sit$sup),
        maxShift = maxShift,
        use_speaq = use_speaq
    )
    if(length(obj$peakList) != np_tot) stop("Lost peaks during alignment")

    # Prepare return object: pcial = integer index, x0al = aligned ppm
    pcial <- obj$peakList[(np_ref+1):np_tot]
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
