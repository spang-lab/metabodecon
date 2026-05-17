
# API #####

#' @export
#'
#' @title Align deconvoluted spectra
#'
#' @description
#' Aligns peaks across a set of deconvoluted spectra by chaining two
#' stages:
#'
#' 1. **CluPA** ([metabodecon::clupa()]) shifts peak centers
#'    continuously toward a reference using hierarchical-clustering FFT
#'    segment shifts (Beirnaert et al. 2018, Vu et al. 2011). Adds
#'    `x0al` and `pcial` (post-CluPA center and cssh column index) to
#'    each peak; original `x0`, `A`, `lambda`, `pcide` are preserved.
#' 2. **RefPA** ([metabodecon::snap_to_ref()]) records, for each peak,
#'    the nearest reference column on `cssh` within `maxCombine` as
#'    `pcisn` / `x0sn`. Peaks farther than `maxCombine` from every
#'    reference column get `pcisn = NA` / `x0sn = NA`. No peaks are
#'    dropped and amplitudes are not summed here — collisions on the
#'    same `pcisn` are aggregated downstream by [metabodecon::si_mat()].
#'
#' RefPA collapses the continuous peak grid produced by CluPA into the
#' reference's discrete peak grid, which is what downstream feature
#' extraction (e.g. [metabodecon::peak_mat()]) expects.
#'
#' @param x A `decons2` (or `aligns`) object.
#' @param ref
#' Optional reference spectrum (`align` or `decon2`). When
#' `NULL` (default) the reference is chosen by [metabodecon::find_ref()].
#' @param maxShift
#' Maximum number of datapoints a peak center may be
#' shifted by CluPA. `maxShift = 0L` skips CluPA (sets `x0al = x0`).
#' @param maxCombine
#' Maximum snap distance for RefPA in chemical-shift
#' columns. `maxCombine = 0L` skips RefPA (no snapping). A negative
#' value is treated as `maxShift`.
#' @param verbose
#' Print progress messages?
#' @param nworkers
#' Number of parallel workers.
#' @param full
#' If `TRUE` also recompute the aligned superposition during
#' CluPA. RefPA always drops `sit$supal` (the post-snap peak list is
#' no longer Lorentz-compatible).
#' @param use_speaq
#' Use `speaq::dohCluster` instead of the bundled CluPA
#' implementation.
#'
#' @return An object of class `aligns`.
#'
#' @author 2024-2026 Tobias Schmidt: initial version.
#'
#' @examples
#' \dontrun{
#'   decons <- deconvolute(sim[1:5], sfr=c(3.55, 3.35))
#'   aligned <- align(decons, maxShift=50, maxCombine=20)
#' }
align <- function(x, ref=NULL, maxShift=50, maxCombine=0,
                  verbose=TRUE, nworkers=1, full=TRUE, use_speaq=FALSE,
                  supsh="lorentz") {
    stopifnot(
        inherits(x, "decons2"),
        is_int(maxShift, 1), maxShift >= 0,
        is_int(maxCombine, 1),
        is_bool(verbose, 1), is_int(nworkers, 1),
        is.null(ref) || inherits(ref, "decon2")
    )
    if (maxCombine < 0L) maxCombine <- as.integer(maxShift)
    a <- clupa(x, ref=ref, maxShift=maxShift, verbose=verbose,
               nworkers=nworkers, full=full, use_speaq=use_speaq,
               supsh=supsh)
    if (maxCombine > 0L) a <- snap_to_ref(a, ref=ref, maxCombine=maxCombine)
    a
}

#' @export
#' @name alignment_funs
#' @rdname alignment_funs
#'
#' @title Alignment building blocks
#'
#' @description
#' Pluggable alignment stages used by [metabodecon::align()] and the
#' `align_fun` argument of [metabodecon::fit_mdm()] /
#' [metabodecon::benchmark()].
#'
#' - [metabodecon::clupa()]: **CluPA** — hierarchical-clustering peak
#'   alignment (recursive FFT segment shifts, Beirnaert et al. 2018, Vu
#'   et al. 2011).
#' - [metabodecon::snap_to_ref()]: **RefPA** — reference-based peak
#'   alignment: snap each peak to the nearest reference column within
#'   `maxCombine`.
#' - [metabodecon::identity_align()]: no-op. Returns its argument unchanged.
#'
#' @param x A `decons2` or `aligns` object.
#' @param ref Optional reference spectrum (`align` or `decon2`). When
#'   `NULL`, chosen by [metabodecon::find_ref()].
#' @param maxShift Maximum CluPA shift in datapoints.
#' @param maxCombine Maximum RefPA snap distance in datapoints.
#' @param verbose Print progress messages?
#' @param nworkers Number of parallel workers.
#' @param full If `TRUE` also recompute the aligned superposition.
#' @param use_speaq Use `speaq::dohCluster` (CluPA only).
#' @param ... Ignored.
#' @return An object of class `aligns`.
clupa <- function(
    x, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1,
    full=TRUE, use_speaq=FALSE, supsh="lorentz"
) {
    supsh <- match.arg(supsh, c("lorentz", "sparse"))
    x <- ensure_cssh(x)
    x <- ensure_supsh(x, supsh=supsh)
    if (!is.null(ref) && is.null(ref$cssh)) ref$cssh <- x[[1]]$cssh
    if (!is.null(ref) && is.null(ref$sit$supsh)) {
        ref$sit$supsh <- make_supsh(ref$cssh, ref$lcpar, supsh)
    }
    if (maxShift == 0L) return(noshift_align(x, full=full))
    ref <- ref %||% find_ref(x)
    aligns <- mcmapply(
        nworkers, align_decon, x,
        MoreArgs=list(ref, maxShift, full=full, use_speaq=use_speaq)
    )
    class(aligns) <- c("aligns", "decons2", "spectra")
    aligns
}

# Make sure every spectrum in `x` carries a shared cssh grid. If none
# is present, derive one from the input `cs` ranges. If some carry
# cssh, require all to agree (otherwise alignment indices wouldn't be
# comparable). Does NOT compute sit$supsh — see ensure_supsh().
ensure_cssh <- function(x) {
    have_cssh <- vapply(x, function(s) !is.null(s$cssh), logical(1))
    if (!all(have_cssh)) {
        cssh <- make_cssh(x)
        for (i in seq_along(x)) x[[i]]$cssh <- cssh
        return(x)
    }
    cssh <- x[[1]]$cssh
    for (i in seq_along(x)) {
        if (!isTRUE(all.equal(x[[i]]$cssh, cssh))) {
            stop("Spectra carry different cssh; alignment requires a shared grid.")
        }
    }
    x
}

# Compute the shared-grid input vector sit$supsh for CluPA. Two modes:
#
#   "lorentz" (default): full Lorentz superposition evaluated at cssh.
#     Smooth, expensive: O(length(cssh) * npeaks).
#
#   "sparse": delta comb — zeros everywhere except at the cssh column
#     nearest each peak center, which carries the peak area A. Cheap:
#     O(npeaks). Lets speaq's FFT cross-correlator align directly on
#     amplitude-weighted peak positions.
#
# Only fills in sit$supsh when it is missing. Callers that want to
# *switch* the cached mode must clear sit$supsh first (or re-run
# deconvolute_spectra with the desired `supsh=`). Assumes
# ensure_cssh() has already run.
ensure_supsh <- function(x, supsh="lorentz") {
    supsh <- match.arg(supsh, c("lorentz", "sparse"))
    for (i in seq_along(x)) {
        if (is.null(x[[i]]$sit$supsh)) {
            x[[i]]$sit$supsh <- make_supsh(x[[i]]$cssh, x[[i]]$lcpar, supsh)
        }
    }
    x
}

# Build the sit$supsh vector from a peak list. Branches on `supsh`
# mode; see ensure_supsh() for semantics. Shared between
# deconvolute_spectrum_{r,rust}() and ensure_supsh().
make_supsh <- function(cssh, lcpar, supsh="lorentz") {
    supsh <- match.arg(supsh, c("lorentz", "sparse"))
    if (supsh == "lorentz") return(lorentz_sup(cssh, lcpar=lcpar))
    out <- numeric(length(cssh))
    if (nrow(lcpar) == 0L) return(out)
    idx <- round(convert_pos(lcpar$x0, cssh, seq_along(cssh)))
    idx <- pmin(pmax(as.integer(idx), 1L), length(cssh))
    s <- tapply(as.numeric(lcpar$A), idx, sum)
    out[as.integer(names(s))] <- s
    out
}

# Integer cssh-column index for each value in `vals`. Clamped to
# [1, length(cssh)]. Used to build pcide / pcial / pcisn from the
# corresponding x0 / x0al / x0sn vectors.
pci_on_cssh <- function(vals, cssh) {
    idx <- round(convert_pos(vals, cssh, seq_along(cssh)))
    pmin(pmax(as.integer(idx), 1L), length(cssh))
}

#' @export
#' @rdname alignment_funs
identity_align <- function(x, ...) x

#' @export
#' @rdname alignment_funs
#'
#' @description
#' [metabodecon::snap_to_ref()] applies the RefPA step on its own: for
#' each peak in each spectrum, finds the nearest reference column on
#' `cssh` and records that column as `pcisn` (and its ppm value as
#' `x0sn`). Peaks farther than `maxCombine` columns from every
#' reference column get `pcisn = NA` / `x0sn = NA`. Original `x0`,
#' `x0al`, `A`, `lambda`, `pcide` and `pcial` are preserved — RefPA
#' only *adds* the snapped fields. Collisions on the same `pcisn`
#' column are not merged here; [metabodecon::si_mat()] sums their
#' areas when rasterising the feature matrix. `sit$supal` is cleared
#' because the post-snap superposition would need recomputing.
#'
#' [metabodecon::combine_peaks()] is an alternative post-CluPA fine
#' tuner that, unlike `snap_to_ref`, does not require the target
#' columns to come from a reference spectrum: it greedily merges
#' neighbouring `cssh` columns whose non-zero rows do not collide,
#' within a window of `maxCombine`. Operates on the cross-spectrum
#' peak-area matrix built from `lcpar$pcial` / `lcpar$A`; sets
#' `lcpar$pcisn` / `lcpar$x0sn` per peak from the discovered column
#' mapping. No peaks are dropped. `sit$supal` is cleared.
snap_to_ref <- function(x, ref=NULL, maxCombine=20, ...) {
    stopifnot(inherits(x, "decons2"), is_int(maxCombine, 1), maxCombine >= 0)
    if (maxCombine == 0L) return(x)
    x <- ensure_cssh(x)
    if (!is.null(ref) && is.null(ref$cssh)) ref$cssh <- x[[1]]$cssh
    ref <- ref %||% find_ref(x)
    if (is.null(ref$lcpar$pcial)) {
        ref$lcpar$pcial <- pci_on_cssh(ref$lcpar$x0, ref$cssh)
    }
    cssh <- x[[1]]$cssh
    nc <- length(cssh)
    pp <- sort(unique(as.integer(ref$lcpar$pcial)))
    pp <- pp[pp >= 1L & pp <= nc]
    for (s in seq_along(x)) {
        x[[s]]$lcpar <- snap_lcpar(x[[s]]$lcpar, pp, maxCombine, cssh)
        x[[s]]$sit$supal <- NULL
        class(x[[s]]) <- c("align", "decon2", "spectrum")
    }
    class(x) <- c("aligns", "decons2", "spectra")
    x
}

# Per-spectrum peak-list snap: add `pcisn` (nearest reference column
# index on cssh) and `x0sn` (= cssh[pcisn]) to each row of `lcpar`,
# keeping all original columns. Peaks farther than `maxCombine` from
# every reference column get pcisn = NA / x0sn = NA. Amplitudes are
# NOT summed here; collisions on the same pcisn are aggregated by
# si_mat() at rasterisation time.
snap_lcpar <- function(lcpar, pp, maxCombine, cssh) {
    n <- nrow(lcpar)
    lcpar$pcisn <- rep(NA_integer_, n)
    lcpar$x0sn  <- rep(NA_real_,    n)
    if (n == 0L || length(pp) == 0L) return(lcpar)
    pcial <- as.integer(lcpar$pcial)
    # For each peak, find nearest reference column (within maxCombine).
    idx <- findInterval(pcial, pp)
    lo <- pmax(idx, 1L); hi <- pmin(idx + 1L, length(pp))
    dlo <- abs(pcial - pp[lo]); dhi <- abs(pcial - pp[hi])
    nearest <- ifelse(dlo <= dhi, pp[lo], pp[hi])
    dist <- pmin(dlo, dhi)
    keep <- dist <= maxCombine
    lcpar$pcisn[keep] <- as.integer(nearest[keep])
    lcpar$x0sn[keep]  <- cssh[nearest[keep]]
    lcpar
}

#' @export
#' @rdname alignment_funs
combine_peaks <- function(x, ref=NULL, maxCombine=20, ...) {
    stopifnot(inherits(x, "decons2"), is_int(maxCombine, 1), maxCombine >= 0)
    if (maxCombine == 0L) return(x)
    x <- ensure_cssh(x)
    cssh <- x[[1]]$cssh
    nc <- length(cssh)
    ns <- length(x)
    # Build cross-spectrum peak-area matrix (rows = spectra, cols = cssh)
    # from each spectrum's pcial / A. Collisions on the same column are
    # summed, mirroring how si_mat rasterizes downstream.
    M <- matrix(0, nrow=ns, ncol=nc)
    for (s in seq_len(ns)) {
        lcpar <- x[[s]]$lcpar
        n <- nrow(lcpar)
        if (n == 0L) next
        if (is.null(lcpar$pcial)) {
            lcpar$pcial <- pci_on_cssh(lcpar$x0, cssh)
            x[[s]]$lcpar <- lcpar
        }
        pcial <- as.integer(lcpar$pcial)
        A <- as.numeric(lcpar$A)
        # tapply -> faster than a per-peak loop when many peaks land
        # on the same column.
        keep <- pcial >= 1L & pcial <= nc
        if (!any(keep)) next
        s_by_col <- tapply(A[keep], pcial[keep], sum)
        M[s, as.integer(names(s_by_col))] <- s_by_col
    }
    map <- combine_peaks_mat(M, maxCombine)$map
    for (s in seq_len(ns)) {
        lcpar <- x[[s]]$lcpar
        n <- nrow(lcpar)
        if (n == 0L) {
            lcpar$pcisn <- integer(0)
            lcpar$x0sn  <- numeric(0)
        } else {
            pcial <- as.integer(lcpar$pcial)
            ok <- pcial >= 1L & pcial <= nc
            pcisn <- rep(NA_integer_, n)
            x0sn  <- rep(NA_real_, n)
            pcisn[ok] <- map[pcial[ok]]
            x0sn[ok]  <- cssh[pcisn[ok]]
            lcpar$pcisn <- pcisn
            lcpar$x0sn  <- x0sn
        }
        x[[s]]$lcpar <- lcpar
        x[[s]]$sit$supal <- NULL
        class(x[[s]]) <- c("align", "decon2", "spectrum")
    }
    class(x) <- c("aligns", "decons2", "spectra")
    x
}

# Greedy column-merge on the cross-spectrum peak-area matrix `M`.
# Two columns may merge only if no row has non-zero entries in both
# (no collision). Within each pass we pick the most "beneficial"
# neighbour (column with the most non-zero entries) within
# `maxCombine` columns. Returns the merged matrix together with a
# per-column destination map: `map[c]` is the column that the original
# column `c` ended up in.
#
# 2021-2024 Wolfram Gronwald: initial version.
# 2024-2025 Tobias Schmidt: refactored initial version.
combine_peaks_mat <- function(M, maxCombine=5, lower_bound=1) {
    U <- M != 0
    uu <- colSums(U)
    nc <- ncol(M)
    map <- seq_len(nc)
    if (nrow(M) <= lower_bound) return(list(M=M, map=map))
    for (i in (nrow(M) - 1):lower_bound) {
        for (j in which(uu == i)) {
            if (uu[j] == 0) next
            nn <- seq(max(1, j - maxCombine), min(nc, j + maxCombine))
            nn <- nn[nn != j]
            if (length(nn) == 0) next
            mj <- M[, j]; uj <- U[, j]
            repeat {
                nn <- nn[uu[nn] > 0]
                if (length(nn) == 0) break
                cc <- combine_scores(U, uu, j, nn, uj=uj)
                if (max(cc) == 0) break
                n <- nn[which.max(cc)]
                mj <- mj + M[, n]
                uj <- uj | U[, n]
                uu[j] <- uu[j] + uu[n]
                M[, n] <- 0; U[, n] <- FALSE
                uu[n] <- 0
                map[map == n] <- j
                nn <- nn[nn != n]
                if (length(nn) == 0) break
            }
            M[, j] <- mj; U[, j] <- uj
        }
    }
    list(M=M, map=map)
}

# Combine score for each neighbour `nn` of column `j`: how many
# non-zero rows the neighbour contributes, or 0 if it collides with
# `j` (i.e. some row is non-zero in both). `U`, `uu`, `uj` are the
# precomputed non-zero mask / per-column counts / column-`j` mask
# (passed in so the caller can update them in place between calls).
combine_scores <- function(U, uu, j, nn, uj=NULL) {
    nn <- nn[nn >= 1 & nn <= ncol(U)]
    if (length(nn) == 0) return(numeric(0))
    if (is.null(uj)) uj <- U[, j]
    overlaps <- .colSums(U[, nn, drop=FALSE] & uj, nrow(U), length(nn))
    cc <- uu[nn]
    cc[overlaps > 0] <- 0
    unname(cc)
}

# Internal #####

# No-op CluPA: set x0al = x0 (no shift) for every spectrum and return
# an aligns object. Used by clupa() when maxShift = 0 so the value is a
# valid grid-search point alongside positive shifts.
noshift_align <- function(x, full=TRUE) {
    aligns <- lapply(x, noshift_one, full=full)
    class(aligns) <- c("aligns", "decons2", "spectra")
    aligns
}

noshift_one <- function(x, full=TRUE) {
    cssh <- x$cssh
    if (is.null(x$lcpar$pcide)) x$lcpar$pcide <- pci_on_cssh(x$lcpar$x0, cssh)
    x$lcpar$x0al <- x$lcpar$x0
    x$lcpar$pcial <- x$lcpar$pcide
    if (full) x$sit$supal <- lorentz_sup(cssh, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

align_decon <- function(x, ref, maxShift, full=TRUE, use_speaq=FALSE) {
    cssh <- x$cssh
    pci_x <- lcpar_pci(x$lcpar, cssh)
    pci_ref <- lcpar_pci(ref$lcpar, cssh)
    np_x <- length(pci_x); np_ref <- length(pci_ref)
    obj <- hclust_align(
        refSpec=ref$sit$supsh, tarSpec=x$sit$supsh,
        peakList=c(pci_ref, pci_x),
        peakLabel=c(rep(1, np_ref), rep(0, np_x)),
        startP=1, endP=length(x$sit$supsh),
        maxShift=maxShift, use_speaq=use_speaq
    )
    if (length(obj$peakList) != np_ref + np_x) stop("Lost peaks during alignment")
    pcial <- obj$peakList[(np_ref + 1):(np_ref + np_x)]
    x$lcpar$x0al <- cssh[pcial]
    x$lcpar$pcial <- pcial
    if (full) x$sit$supal <- lorentz_sup(cssh, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

find_ref <- function(x) {
    x <- ensure_cssh(x)
    cssh <- x[[1]]$cssh
    pci <- lapply(x, function(s) lcpar_pci(s$lcpar, cssh))
    x[[find_ref_ind(pci)$refInd]]
}

# Datapoint indices on `cssh` for the peaks in `lcpar`. Prefers the
# cached `pcide` (set at deconvolution time); otherwise computes from
# `x0`; otherwise falls back to `pcial`. The latter two paths exist
# only for backwards compatibility with objects saved before `pcide`
# was added.
lcpar_pci <- function(lcpar, cssh) {
    pcide <- lcpar[["pcide"]]
    if (!is.null(pcide)) return(as.integer(pcide))
    x0 <- lcpar[["x0"]]
    if (!is.null(x0)) return(pci_on_cssh(x0, cssh))
    as.integer(lcpar[["pcial"]])
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
