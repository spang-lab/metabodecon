
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
#'    `x0al` and `pcial` (post-CluPA center and column index) to each
#'    peak; original `x0`, `A`, `lambda`, `pcide` are preserved.
#' 2. **RefPA** ([metabodecon::snap_to_ref()]) records, for each peak,
#'    the nearest reference column within `maxCombine` as `pcisn` /
#'    `x0sn`. Peaks farther than `maxCombine` from every reference
#'    column get `pcisn = NA` / `x0sn = NA`. No peaks are dropped and
#'    amplitudes are not summed here — collisions on the same `pcisn`
#'    are aggregated downstream by [metabodecon::si_mat()].
#'
#' All spectra in `x` must already live on the same chemical-shift
#' grid (identical `$cs` vector across spectra). Call
#' [metabodecon::harmonize_grid()] upstream if your inputs come from
#' different acquisitions with slight calibration offsets.
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
#' Use `speaq::hClustAlign` instead of the bundled CluPA
#' implementation. Defaults to `FALSE`; the bundled implementation is
#' byte-equivalent to the speaq one (see `tests/testthat/test-speaq.R`).
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
align <- function(x, y=NULL, ref=NULL, maxShift=50, maxCombine=0,
                  verbose=TRUE, nworkers=1, full=TRUE, use_speaq=FALSE,
                  gap_tol=NULL) {
    stopifnot(
        inherits(x, "decons2"),
        is_int(maxShift, 1), maxShift >= 0,
        is_int(maxCombine, 1),
        is_bool(verbose, 1), is_int(nworkers, 1),
        is.null(ref) || inherits(ref, "decon2"),
        is.null(y) || (is.factor(y) && length(y) == length(x))
    )
    if (maxCombine < 0L) maxCombine <- as.integer(maxShift)
    a <- clupa(
        x, y=y, ref=ref, maxShift=maxShift, verbose=verbose,
        nworkers=nworkers, full=full, use_speaq=use_speaq,
        gap_tol=gap_tol
    )
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
#'   et al. 2011). Operates on the Lorentz reconstruction `sit$sup`
#'   already attached to each spectrum by deconvolution.
#' - [metabodecon::snap_to_ref()]: **RefPA** — reference-based peak
#'   alignment: snap each peak to the nearest reference column within
#'   `maxCombine`.
#' - [metabodecon::identity_align()]: no-op. Returns its argument unchanged.
#'
#' All these functions require every input spectrum to share the same
#' `$cs` grid; an explicit `stop()` is raised otherwise. Call
#' [metabodecon::harmonize_grid()] upstream to enforce that invariant.
#'
#' @param x A `decons2` or `aligns` object.
#' @param ref Optional reference spectrum (`align` or `decon2`). When
#'   `NULL`, chosen by [metabodecon::find_ref()].
#' @param maxShift Maximum CluPA shift in datapoints.
#' @param maxCombine Maximum RefPA snap distance in datapoints.
#' @param verbose Print progress messages?
#' @param nworkers Number of parallel workers.
#' @param full If `TRUE` also recompute the aligned superposition.
#' @param use_speaq Use `speaq::hClustAlign` (CluPA only).
#' @param ... Ignored.
#' @return An object of class `aligns`.
clupa <- function(
    x, y=NULL, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1,
    full=TRUE, use_speaq=FALSE, gap_tol=NULL
) {
    # 1) Assert shared grid.
    cs <- ensure_shared_cs(x)

    # 2) Resolve reference. If none supplied, pick via find_ref (or
    #    build_clupa_consensus when class labels are provided).
    if (is.null(ref)) {
        ref <- if (is.null(y)) find_ref(x) else build_clupa_consensus(
            x, y, maxShift=maxShift, use_speaq=use_speaq, gap_tol=gap_tol
        )
    }
    # The reference must live on the same grid as the rest of x. A
    # consensus carries no $cs of its own (it's built on the shared
    # grid); use cs in that case.
    if (is.null(ref$cs)) {
        ref$cs <- cs
    } else if (!isTRUE(all.equal(ref$cs, cs))) {
        stop(
            "clupa: reference cs does not match the shared cs of x. ",
            "Call harmonize_grid(x, target=ref$cs) first.", call.=FALSE
        )
    }

    # 3) Ensure every spectrum carries pcide and sit$sup.
    x <- mcmapply(nworkers, ensure_align_aux, x, MoreArgs=list(cs=cs))
    if (is.null(ref$lcpar$pcide)) ref$lcpar$pcide <- pci_on_cs(ref$lcpar$x0, cs)
    if (is.null(ref$sit$sup)) ref$sit$sup <- lorentz_sup(cs, lcpar=ref$lcpar)

    # 4) Align (or short-circuit when no shift is requested).
    aligns <- if (maxShift == 0L) {
        noshift_align(x, full=full)
    } else {
        mcmapply(
            nworkers, align_decon, x,
            MoreArgs=list(ref, maxShift, full=full, use_speaq=use_speaq)
        )
    }
    class(aligns) <- c("aligns", "decons2", "spectra")
    attr(aligns, "ref") <- ref
    aligns
}

# Build a CluPA-aligned class consensus reference.
#
# 1) Pick one representative per class via find_ref().
# 2) Run a 2-pass clupa() over the K representatives so their peak
#    lists are on a common scale before union.
# 3) Union all aligned peak lists, then collapse near-coincident peaks
#    within `gap_tol` ppm via dedupe_peaks() (amplitude-weighted x0,
#    arithmetic-mean A / lambda).
# 4) Return a `consensus` object shaped like an aligned spectrum so
#    align_decon can use it as a reference (carries lcpar with pcide,
#    sit$sup).
build_clupa_consensus <- function(x, y, maxShift, use_speaq, gap_tol=NULL) {
    stopifnot(is.factor(y), length(y) == length(x))
    cs <- ensure_shared_cs(x)
    if (is.null(gap_tol)) gap_tol <- 2 * abs(cs[2] - cs[1])

    lvs <- levels(y)
    reps <- lapply(lvs, function(lv) {
        ix <- which(y == lv); if (length(ix) == 0L) NULL else find_ref(x[ix])
    })
    reps <- reps[!vapply(reps, is.null, logical(1))]
    if (length(reps) <= 1L) {
        return(if (length(reps) == 1L) reps[[1]] else find_ref(x))
    }

    class(reps) <- c("decons2", "spectra")
    reps_al <- clupa(reps, maxShift=maxShift, verbose=FALSE, nworkers=1,
                     full=FALSE, use_speaq=use_speaq)

    pos_of <- function(s) {
        p <- s$lcpar$x0al
        if (is.null(p)) p <- s$lcpar$x0
        as.numeric(p)
    }
    union_lcpar <- data.frame(
        x0     = unlist(lapply(reps_al, pos_of)),
        A      = unlist(lapply(reps_al, function(s) as.numeric(s$lcpar$A))),
        lambda = unlist(lapply(reps_al, function(s) as.numeric(s$lcpar$lambda)))
    )
    union_lcpar <- dedupe_peaks(union_lcpar, gap_tol)
    union_lcpar$pcide <- pci_on_cs(union_lcpar$x0, cs)

    sit <- list(sup=lorentz_sup(cs, lcpar=union_lcpar))
    structure(
        list(cs=cs, lcpar=union_lcpar, sit=sit),
        class=c("consensus", "align", "decon2", "spectrum")
    )
}

# Assert that every spectrum in `x` shares the same $cs grid. Returns
# the shared cs vector. Stops with an actionable error message
# otherwise; the recommended remedy is to call harmonize_grid(x)
# upstream.
ensure_shared_cs <- function(x) {
    stopifnot(is.list(x), length(x) > 0L)
    cs <- x[[1]]$cs
    if (is.null(cs)) stop(
        "Spectra have no $cs vector; cannot align.", call.=FALSE
    )
    for (i in seq_along(x)) {
        if (!isTRUE(all.equal(x[[i]]$cs, cs))) stop(
            "Spectra do not share a common chemical-shift grid. ",
            "Call harmonize_grid(x) to bring them onto a single grid ",
            "before alignment.", call.=FALSE
        )
    }
    cs
}

# Per-spectrum auxiliary fields needed for CluPA: pcide (datapoint
# index of each peak on the shared grid) and sit$sup (Lorentz
# reconstruction). Both are usually already attached at deconvolution
# time; this helper fills them in for callers that built a decons2 by
# hand.
ensure_align_aux <- function(s, cs) {
    if (is.null(s$lcpar$pcide)) s$lcpar$pcide <- pci_on_cs(s$lcpar$x0, cs)
    if (is.null(s$sit$sup)) s$sit$sup <- lorentz_sup(cs, lcpar=s$lcpar)
    s
}

# Integer column index for each value in `vals`, computed against `cs`.
# Clamped to [1, length(cs)]. Used to build pcide / pcial / pcisn from
# the corresponding x0 / x0al / x0sn vectors.
pci_on_cs <- function(vals, cs) {
    idx <- round(convert_pos(vals, cs, seq_along(cs)))
    pmin(pmax(as.integer(idx), 1L), length(cs))
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
#' the shared `cs` grid and records that column as `pcisn` (and its
#' ppm value as `x0sn`). Peaks farther than `maxCombine` columns from
#' every reference column get `pcisn = NA` / `x0sn = NA`. Original
#' `x0`, `x0al`, `A`, `lambda`, `pcide` and `pcial` are preserved —
#' RefPA only *adds* the snapped fields. Collisions on the same
#' `pcisn` column are not merged here; [metabodecon::si_mat()] sums
#' their areas when rasterising the feature matrix. `sit$supal` is
#' cleared because the post-snap superposition would need recomputing.
#'
#' [metabodecon::combine_peaks()] is an alternative post-CluPA fine
#' tuner that, unlike `snap_to_ref`, does not require the target
#' columns to come from a reference spectrum: it greedily merges
#' neighbouring `cs` columns whose non-zero rows do not collide,
#' within a window of `maxCombine`. Operates on the cross-spectrum
#' peak-area matrix built from `lcpar$pcial` / `lcpar$A`; sets
#' `lcpar$pcisn` / `lcpar$x0sn` per peak from the discovered column
#' mapping. No peaks are dropped. `sit$supal` is cleared.
snap_to_ref <- function(x, ref=NULL, maxCombine=20, ...) {
    stopifnot(inherits(x, "decons2"), is_int(maxCombine, 1), maxCombine >= 0)
    if (maxCombine == 0L) return(x)
    cs <- ensure_shared_cs(x)
    ref <- ref %||% find_ref(x)
    if (is.null(ref$lcpar$pcial)) {
        ref$lcpar$pcial <- pci_on_cs(ref$lcpar$x0, ref$cs %||% cs)
    }
    nc <- length(cs)
    pp <- sort(unique(as.integer(ref$lcpar$pcial)))
    pp <- pp[pp >= 1L & pp <= nc]
    for (s in seq_along(x)) {
        x[[s]]$lcpar <- snap_lcpar(x[[s]]$lcpar, pp, maxCombine, cs)
        x[[s]]$sit$supal <- NULL
        class(x[[s]]) <- c("align", "decon2", "spectrum")
    }
    class(x) <- c("aligns", "decons2", "spectra")
    attr(x, "ref") <- ref
    x
}

# Per-spectrum peak-list snap: add `pcisn` (nearest reference column
# index) and `x0sn` (= cs[pcisn]) to each row of `lcpar`, keeping all
# original columns. Peaks farther than `maxCombine` from every
# reference column get pcisn = NA / x0sn = NA. Amplitudes are NOT
# summed here; collisions on the same pcisn are aggregated by
# si_mat() at rasterisation time.
snap_lcpar <- function(lcpar, pp, maxCombine, cs) {
    n <- nrow(lcpar)
    lcpar$pcisn <- rep(NA_integer_, n)
    lcpar$x0sn  <- rep(NA_real_,    n)
    if (n == 0L || length(pp) == 0L) return(lcpar)
    pcial <- as.integer(lcpar$pcial)
    idx <- findInterval(pcial, pp)
    lo <- pmax(idx, 1L); hi <- pmin(idx + 1L, length(pp))
    dlo <- abs(pcial - pp[lo]); dhi <- abs(pcial - pp[hi])
    nearest <- ifelse(dlo <= dhi, pp[lo], pp[hi])
    dist <- pmin(dlo, dhi)
    keep <- dist <= maxCombine
    lcpar$pcisn[keep] <- as.integer(nearest[keep])
    lcpar$x0sn[keep]  <- cs[nearest[keep]]
    lcpar
}

#' @export
#' @rdname alignment_funs
combine_peaks <- function(x, ref=NULL, maxCombine=20, ...) {
    stopifnot(inherits(x, "decons2"), is_int(maxCombine, 1), maxCombine >= 0)
    if (maxCombine == 0L) return(x)
    cs <- ensure_shared_cs(x)
    nc <- length(cs)
    ns <- length(x)
    M <- matrix(0, nrow=ns, ncol=nc)
    for (s in seq_len(ns)) {
        lcpar <- x[[s]]$lcpar
        n <- nrow(lcpar)
        if (n == 0L) next
        if (is.null(lcpar$pcial)) {
            lcpar$pcial <- pci_on_cs(lcpar$x0, cs)
            x[[s]]$lcpar <- lcpar
        }
        pcial <- as.integer(lcpar$pcial)
        A <- as.numeric(lcpar$A)
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
            x0sn[ok]  <- cs[pcisn[ok]]
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
    cs <- x$cs
    if (is.null(x$lcpar$pcide)) x$lcpar$pcide <- pci_on_cs(x$lcpar$x0, cs)
    x$lcpar$x0al <- x$lcpar$x0
    x$lcpar$pcial <- x$lcpar$pcide
    if (full) x$sit$supal <- lorentz_sup(cs, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

# Per-spectrum CluPA kernel. Operates on the shared `cs` grid via
# `x$cs` (which equals `ref$cs` by the time we get here). The FFT
# input is `x$sit$sup`, the Lorentz reconstruction already attached at
# deconvolution time — this is the speaq-equivalent shape (matches the
# `get_sup_mat(decons2)` input that v1.7.0 fed to `dohCluster`).
align_decon <- function(x, ref, maxShift, full=TRUE, use_speaq=FALSE) {
    cs <- x$cs
    pci_x <- lcpar_pci(x$lcpar, cs)
    pci_ref <- lcpar_pci(ref$lcpar, cs)
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
    x$lcpar$x0al <- cs[pcial]
    x$lcpar$pcial <- pcial
    if (full) x$sit$supal <- lorentz_sup(cs, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

find_ref <- function(x) {
    # Compare candidate references in ppm space — no shared grid is
    # needed because `find_ref_ind` only uses pairwise distances, and
    # those are directly comparable across spectra with different `cs`
    # grids when expressed in ppm.
    x0 <- lapply(x, function(s) as.numeric(s$lcpar$x0))
    x[[find_ref_ind(x0)$refInd]]
}

# Datapoint indices on `cs` for the peaks in `lcpar`. Prefers the
# cached `pcide` (set at deconvolution time); otherwise computes from
# `x0`; otherwise falls back to `pcial`. The latter two paths exist
# only for backwards compatibility with objects saved before `pcide`
# was added.
lcpar_pci <- function(lcpar, cs) {
    pcide <- lcpar[["pcide"]]
    if (!is.null(pcide)) return(as.integer(pcide))
    x0 <- lcpar[["x0"]]
    if (!is.null(x0)) return(pci_on_cs(x0, cs))
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
    refSpec, tarSpec, peakList, peakLabel, startP, endP, maxShift,
    use_speaq = FALSE
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

    startCheckP <- startP + which.min(tarSpec[startP:(minPk - 1L)]) - 1L
    if (is.na(startCheckP) || startCheckP < 1L) startCheckP <- startP
    endCheckP <- maxPk + which.min(tarSpec[(maxPk + 1L):endP])
    if (is.na(endCheckP) || endCheckP > length(tarSpec)) endCheckP <- endP

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
        ok <- (adj$stepAdj < 0 && adj$stepAdj + minPk >= startCheckP) ||
              (adj$stepAdj > 0 && adj$stepAdj + maxPk <= endCheckP)
        if (ok) {
            tar_idx <- which(peakLabel == 0L)
            peakList[tar_idx] <- peakList[tar_idx] + adj$stepAdj
            lost <- which(peakList <= 0L | peakList > length(tarSpec))
            if (length(lost) > 0L) {
                peakList <- peakList[-lost]
                peakLabel <- peakLabel[-lost]
            }
            seg <- tarSpec[startCheckP:endCheckP]
            tarSpec[startCheckP:endCheckP] <- do_shift(seg, adj$stepAdj)
        }
    }

    if (length(peakList) < 3L) {
        return(list(tarSpec = tarSpec, peakList = peakList))
    }

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

    if (max1 < min2) {
        endP1 <- max1 + which.min(tarSpec[(max1 + 1L):(min2 - 1L)])
        if (is.na(endP1) || endP1 > length(tarSpec)) endP1 <- max1
        startP2 <- endP1 + 1L
    } else {
        tmp_set <- left_set; left_set <- right_set
        right_set <- tmp_set
        sub1 <- peakList[left_set]; lab1 <- peakLabel[left_set]; id1 <- left_set
        sub2 <- peakList[right_set]; lab2 <- peakLabel[right_set]; id2 <- right_set
        max1 <- max(sub1); min2 <- min(sub2)
        endP1 <- max1 + which.min(tarSpec[(max1 + 1L):(min2 - 1L)])
        if (is.na(endP1) || endP1 > length(tarSpec)) endP1 <- max1
        startP2 <- endP1 + 1L
    }
    if (length(unique(lab1)) > 1L) {
        res <- hclust_align(refSpec, tarSpec, sub1, lab1, startP, endP1,
                            maxShift, use_speaq=use_speaq)
        tarSpec <- res$tarSpec
        peakList[id1] <- pad_peaks(res$peakList, length(id1))
    }
    if (length(unique(lab2)) > 1L) {
        res <- hclust_align(refSpec, tarSpec, sub2, lab2, startP2, endP,
                            maxShift, use_speaq=use_speaq)
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

# Needleman-Wunsch alignment #####

#' @export
#' @rdname alignment_funs
#'
#' @title Pairwise Needleman-Wunsch snap of peak lists onto a reference
#'
#' @description
#' [metabodecon::snap_nw()] is a drop-in alternative to
#' [metabodecon::snap_to_ref()] that aligns each spectrum's peak list to
#' the reference's peak list by global pairwise Needleman-Wunsch on `x0`
#' (chemical shift, ppm). For each spectrum, matched peaks inherit the
#' reference's `cs` column index as their `pcisn`; unmatched peaks get
#' `pcisn = NA`. The match cost is `|x0_spec - x0_ref|` and the gap cost is
#' a constant `gap_tol` (ppm), so a match is rejected in favour of two gaps
#' whenever its position difference exceeds `2 * gap_tol`.
#'
#' Unlike [metabodecon::snap_to_ref()], which is a per-peak nearest-neighbour
#' lookup (greedy, many-to-one), `snap_nw` enforces 1-to-1 pairing via the
#' Needleman-Wunsch DP. This avoids the failure mode where two adjacent
#' spectrum peaks both collapse onto a single reference column.
#'
#' The DP runs in C (`align_dp_c` in `src/align_dp.c`); the R wrapper
#' builds the cost matrix via [base::outer] and dispatches one `.Call`
#' per spectrum.
#'
#' @param x A `decons2` / `aligns` object.
#' @param ref Reference. Either an `aligns`/`decons2` *spectrum* with an
#'   `lcpar` element (in which case `lcpar$x0` is the alignment target), or
#'   the consensus list returned by [metabodecon::build_consensus()] (which
#'   carries the same `lcpar` shape). If `NULL`, picked via
#'   [metabodecon::find_ref()] on `x`.
#' @param gap_tol Numeric scalar. Matching tolerance in ppm. Default 0.02.
#' @param pos_field Which `lcpar` column on `x` to use as the peak position
#'   for alignment. `"x0"` (default) is the raw deconvolution position;
#'   `"x0al"` is the post-CluPA aligned position. The reference is always
#'   matched against on its `lcpar$x0`, so callers building a reference
#'   from aligned data should store the aligned positions in its `x0`
#'   column (which is what [metabodecon::build_consensus()] does when
#'   given `pos_field="x0al"`).
#' @param ... Ignored (signature compatibility with `snap_to_ref`).
#'
#' @return An `aligns` object with `pcisn` and `x0sn` set on each spectrum's
#'   `lcpar`. `sit$supal` is cleared (the snap output is no longer
#'   Lorentz-compatible until rebuilt downstream).
#'
#' @author 2026 Tobias Schmidt: initial version.
snap_nw <- function(x, ref=NULL, gap_tol=0.02, pos_field="x0", w_A=0, ...) {
    stopifnot(inherits(x, "decons2"), is_num(gap_tol, 1), gap_tol >= 0,
              is_num(w_A, 1), w_A >= 0)
    if (gap_tol == 0) return(x)
    cs <- ensure_shared_cs(x)
    ref <- ref %||% find_ref(x)
    if (is.null(ref$lcpar$pcide)) {
        ref$lcpar$pcide <- pci_on_cs(ref$lcpar$x0, ref$cs %||% cs)
    }
    ro   <- order(ref$lcpar$x0)
    rx0  <- as.numeric(ref$lcpar$x0[ro])
    rcol <- as.integer(ref$lcpar$pcide[ro])
    rA   <- if (w_A > 0) normalize_A(as.numeric(ref$lcpar$A[ro])) else NULL
    for (s in seq_along(x)) {
        x[[s]]$lcpar <- snap_nw_lcpar(x[[s]]$lcpar, rx0, rcol, cs, gap_tol,
                                       pos_field=pos_field, w_A=w_A, rA=rA)
        x[[s]]$sit$supal <- NULL
        class(x[[s]]) <- c("align", "decon2", "spectrum")
    }
    class(x) <- c("aligns", "decons2", "spectra")
    x
}

normalize_A <- function(A) {
    if (length(A) == 0L) return(A)
    pos <- A[A > 0]
    if (length(pos) == 0L) return(A)
    A / stats::median(pos)
}

snap_nw_lcpar <- function(lcpar, rx0, rcol, cs, gap_tol, pos_field="x0",
                           w_A=0, rA=NULL) {
    n <- nrow(lcpar)
    lcpar$pcisn <- rep(NA_integer_, n)
    lcpar$x0sn  <- rep(NA_real_,    n)
    if (n == 0L || length(rcol) == 0L) return(lcpar)
    pos <- lcpar[[pos_field]] %||% lcpar$x0
    o   <- order(pos)
    sx0 <- as.numeric(pos[o])
    M   <- abs(outer(sx0, rx0, "-"))
    if (w_A > 0 && !is.null(rA) && !is.null(lcpar$A)) {
        sA <- normalize_A(as.numeric(lcpar$A[o]))
        eps <- 1e-12
        ratio <- outer(pmax(sA, eps), pmax(rA, eps), "/")
        M_amp <- abs(log(ratio))
        M <- M + w_A * gap_tol * M_amp
    }
    storage.mode(M) <- "double"
    gp  <- rep_len(as.double(gap_tol), length(sx0))
    gq  <- rep_len(as.double(gap_tol), length(rx0))
    ans <- .Call(align_dp_c, M, gp, gq)
    al  <- ans$alignment
    mt  <- !is.na(al[, 1]) & !is.na(al[, 2])
    if (any(mt)) {
        si <- o[al[mt, 1]]
        ci <- rcol[al[mt, 2]]
        lcpar$pcisn[si] <- ci
        lcpar$x0sn[si]  <- cs[ci]
    }
    lcpar
}

#' @export
#'
#' @title Build a consensus peak-list reference for NW snapping
#'
#' @description
#' Constructs a single consensus peak list to be used as the alignment
#' target by [metabodecon::snap_nw()]. The consensus represents the union
#' of training-set peak positions, deduplicated within `gap_tol` ppm so
#' near-coincident peaks collapse to a single column.
#'
#' When `y` is supplied, the consensus is class-aware: one
#' [metabodecon::find_ref()] is picked per level of `y`, the per-class
#' references are merged into a seed consensus (so class-specific peaks
#' are present from the start), and every training spectrum is then
#' NW-snapped to this seed; their `x0` values contribute additional
#' consensus positions wherever they did not match the seed. The final
#' consensus is again deduplicated within `gap_tol`.
#'
#' The returned object is shaped like an aligned spectrum
#' (`list(cs, lcpar, ...)`) so it can be passed directly as `ref` to
#' [metabodecon::snap_nw()] at prediction time.
#'
#' @param x A `decons2` (or `aligns`) object.
#' @param y Optional factor of class labels (length == `length(x)`). If
#'   supplied, the seed consensus is built per-class.
#' @param gap_tol Numeric scalar. Tolerance in ppm for both the per-spectrum
#'   NW snap and for deduplicating consensus positions. Default 0.02.
#' @param pos_field Which `lcpar` column to use as the peak position.
#'   `"x0"` (default) for raw deconvolution positions; `"x0al"` to build
#'   the consensus from CluPA-aligned positions (then the consensus lives
#'   in aligned space and should be snapped against using `pos_field="x0al"`).
#'
#' @return A list with components `cs`, `lcpar` (data frame with `x0`,
#'   `A`, `lambda`, `pcide`), plus the class attributes
#'   `c("consensus", "align", "decon2", "spectrum")` so it behaves like a
#'   single-spectrum reference for downstream code.
#'
#' @author 2026 Tobias Schmidt: initial version.
build_consensus <- function(x, y=NULL, gap_tol=0.02, pos_field="x0") {
    stopifnot(inherits(x, "decons2"), is_num(gap_tol, 1), gap_tol > 0)
    cs <- ensure_shared_cs(x)

    swap_field <- function(xx, fld) {
        if (fld == "x0") return(xx)
        for (s in seq_along(xx)) {
            lc <- xx[[s]]$lcpar
            xx[[s]]$lcpar$x0 <- lc[[fld]] %||% lc$x0
        }
        xx
    }
    x_pos <- swap_field(x, pos_field)

    if (is.null(y)) {
        seed <- find_ref(x_pos)
    } else {
        stopifnot(is.factor(y), length(y) == length(x_pos))
        lvs <- levels(y)
        reps <- lapply(lvs, function(lv) {
            ix <- which(y == lv)
            if (length(ix) == 0L) return(NULL)
            find_ref(x_pos[ix])
        })
        reps <- reps[!vapply(reps, is.null, logical(1))]
        seed_x0  <- unlist(lapply(reps, function(r) r$lcpar$x0))
        seed_A   <- unlist(lapply(reps, function(r) r$lcpar$A))
        seed_lam <- unlist(lapply(reps, function(r) r$lcpar$lambda))
        seed_lcpar <- data.frame(x0=seed_x0, A=seed_A, lambda=seed_lam)
        seed_lcpar <- dedupe_peaks(seed_lcpar, gap_tol)
        seed <- list(cs=cs, lcpar=seed_lcpar)
    }
    if (is.null(seed$lcpar$pcide)) {
        seed$lcpar$pcide <- pci_on_cs(seed$lcpar$x0, cs)
    }

    snapped <- snap_nw(x_pos, ref=seed, gap_tol=gap_tol)
    extra_x0  <- c(); extra_A <- c(); extra_lam <- c()
    for (s in seq_along(snapped)) {
        lc <- snapped[[s]]$lcpar
        un <- is.na(lc$pcisn)
        if (any(un)) {
            extra_x0  <- c(extra_x0,  lc$x0[un])
            extra_A   <- c(extra_A,   lc$A[un])
            extra_lam <- c(extra_lam, lc$lambda[un])
        }
    }
    full_lcpar <- data.frame(
        x0=c(seed$lcpar$x0, extra_x0),
        A=c(seed$lcpar$A,   extra_A),
        lambda=c(seed$lcpar$lambda, extra_lam)
    )
    full_lcpar <- dedupe_peaks(full_lcpar, gap_tol)
    full_lcpar$pcide <- pci_on_cs(full_lcpar$x0, cs)
    structure(
        list(cs=cs, lcpar=full_lcpar),
        class=c("consensus", "align", "decon2", "spectrum")
    )
}

# Deduplicate a peak list: peaks within `gap_tol` ppm of each other are
# collapsed into a single peak whose x0 / A / lambda are amplitude-weighted
# means of the cluster.
dedupe_peaks <- function(lcpar, gap_tol) {
    n <- nrow(lcpar)
    if (n == 0L) return(lcpar)
    o <- order(lcpar$x0)
    x0 <- lcpar$x0[o]; A <- lcpar$A[o]; lam <- lcpar$lambda[o]
    g <- c(0, cumsum(diff(x0) > gap_tol))
    keep_x0  <- vapply(split(seq_len(n), g), function(i) {
        w <- A[i]; if (sum(w) <= 0) mean(x0[i]) else stats::weighted.mean(x0[i], w)
    }, numeric(1))
    keep_A   <- vapply(split(seq_len(n), g), function(i) mean(A[i]),   numeric(1))
    keep_lam <- vapply(split(seq_len(n), g), function(i) mean(lam[i]), numeric(1))
    data.frame(x0=as.numeric(keep_x0), A=as.numeric(keep_A),
               lambda=as.numeric(keep_lam))
}
