
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
align <- function(x, y=NULL, ref=NULL, maxShift=50, maxCombine=0,
                  verbose=TRUE, nworkers=1, full=TRUE, use_speaq=FALSE,
                  supsh="triangle", shift_method="auto", gap_tol=NULL) {
    stopifnot(
        inherits(x, "decons2"),
        is_int(maxShift, 1), maxShift >= 0,
        is_int(maxCombine, 1),
        is_bool(verbose, 1), is_int(nworkers, 1),
        is.null(ref) || inherits(ref, "decon2"),
        is.null(y) || (is.factor(y) && length(y) == length(x))
    )
    if (maxCombine < 0L) maxCombine <- as.integer(maxShift)
    a <- clupa(x, y=y, ref=ref, maxShift=maxShift, verbose=verbose,
               nworkers=nworkers, full=full, use_speaq=use_speaq,
               supsh=supsh, shift_method=shift_method, gap_tol=gap_tol)
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
    x, y=NULL, ref=NULL, maxShift=50, verbose=TRUE, nworkers=1,
    full=TRUE, use_speaq=FALSE, supsh="triangle", shift_method="auto",
    gap_tol=NULL
) {
    supsh <- match.arg(
        supsh, c("triangle", "rectangle", "lorentz", "sparse", "eiffel")
    )
    shift_method <- match.arg(shift_method, c("auto", "recompute", "slide"))
    sm <- resolve_shift_method(shift_method, supsh)

    # 1) Resolve reference. The ref's grid drives the rest of the call.
    if (is.null(ref)) {
        ref <- if (is.null(y)) find_ref(x) else build_clupa_consensus(
            x, y, maxShift=maxShift, supsh=supsh, shift_method=sm,
            use_speaq=use_speaq, gap_tol=gap_tol
        )
    }
    cssh <- ref$cssh %||% ref$cs

    # 2) Bind every spectrum to ref's cssh and compute sit$supsh in
    #    parallel. After this, x[[i]]$cssh == cssh for all i.
    x <- mcmapply(
        nworkers, bind_to_cssh, x,
        MoreArgs=list(cssh=cssh, supsh=supsh)
    )

    # 3) Backfill ref so align_decon can read cssh / pcide / supsh.
    if (is.null(ref$cssh)) ref$cssh <- cssh
    if (is.null(ref$lcpar$pcide)) {
        ref$lcpar$pcide <- pci_on_cssh(ref$lcpar$x0, cssh)
    }
    if (is.null(ref$sit$supsh)) {
        ref$sit$supsh <- make_supsh(cssh, ref$lcpar, supsh)
    }

    # 4) Align (or short-circuit when no shift is requested).
    aligns <- if (maxShift == 0L) {
        noshift_align(x, full=full)
    } else {
        mcmapply(
            nworkers, align_decon, x,
            MoreArgs=list(ref, maxShift, full=full, use_speaq=use_speaq,
                          supsh=supsh, shift_method=sm)
        )
    }
    class(aligns) <- c("aligns", "decons2", "spectra")
    attr(aligns, "ref") <- ref
    aligns
}

# Stamp `cssh` onto a spectrum, recompute the per-peak cssh column
# index `pcide`, and rebuild `sit$supsh` against the new grid.
bind_to_cssh <- function(s, cssh, supsh) {
    s$cssh <- cssh
    s$lcpar$pcide <- pci_on_cssh(s$lcpar$x0, cssh)
    s$sit$supsh <- make_supsh(cssh, s$lcpar, supsh)
    s
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
#    align_decon can use it as a reference (carries cssh, lcpar with
#    pcide, sit$supsh).
#
# Compared to picking a single class-A spectrum as the reference,
# this guarantees that peaks present only in class B still have a
# matching column on the reference grid — so class-distinguishing
# peaks are no longer silently snapped to whichever class-A column
# happens to be nearest.
build_clupa_consensus <- function(x, y, maxShift, supsh, shift_method,
                                  use_speaq, gap_tol=NULL) {
    stopifnot(is.factor(y), length(y) == length(x))
    # Anchor every per-class subset on the same cssh so the inner clupa
    # call sees a single consistent grid. clupa() already does this
    # before delegating here, but a direct caller (e.g. paper-repo
    # benchmark code) may have skipped it.
    x <- ensure_cssh(x)
    cssh <- x[[1]]$cssh
    if (is.null(gap_tol)) gap_tol <- 2 * abs(cssh[2] - cssh[1])

    # Pick one rep per class. Empty classes are skipped silently.
    lvs <- levels(y)
    reps <- lapply(lvs, function(lv) {
        ix <- which(y == lv); if (length(ix) == 0L) NULL else find_ref(x[ix])
    })
    reps <- reps[!vapply(reps, is.null, logical(1))]
    if (length(reps) <= 1L) {
        # 0 or 1 class with members: fall back to a single rep.
        return(if (length(reps) == 1L) reps[[1]] else find_ref(x))
    }

    # Align the reps to one of themselves (find_ref picks). Use the
    # same supsh / shift_method / use_speaq as the outer call so the
    # consensus is built on the same correlation geometry CluPA will
    # ultimately use.
    class(reps) <- c("decons2", "spectra")
    reps_al <- clupa(reps, maxShift=maxShift, verbose=FALSE, nworkers=1,
                     full=FALSE, use_speaq=use_speaq,
                     supsh=supsh, shift_method=shift_method)

    # Union of aligned peak lists. Use x0al when available (aligned
    # position), x0 otherwise (the chosen reference inside reps_al).
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
    union_lcpar$pcide <- pci_on_cssh(union_lcpar$x0, cssh)

    sit <- list(supsh=make_supsh(cssh, union_lcpar, supsh))
    structure(
        list(cssh=cssh, lcpar=union_lcpar, sit=sit),
        class=c("consensus", "align", "decon2", "spectrum")
    )
}

# Pick the concrete shift method when `auto` is requested. The
# peakList-driven rebuild beats slide-and-pad correctness-wise for any
# shape, but is asymptotically slower than slide for lorentz (which
# touches every column anyway). Sparse shapes (triangle / rectangle /
# sparse) are *cheaper* to rebuild than to slide, so auto picks them.
resolve_shift_method <- function(method, supsh) {
    if (method != "auto") return(method)
    if (supsh == "lorentz") "slide" else "recompute"
}

# Make sure every spectrum in `x` carries a shared cssh grid and a
# cached cssh-column index `lcpar$pcide` for each peak. If no spectrum
# carries cssh, derive one from the input `cs` ranges. If some carry
# cssh, require all to agree (otherwise alignment indices wouldn't be
# comparable). Does NOT compute sit$supsh — see ensure_supsh().
ensure_cssh <- function(x) {
    have_cssh <- vapply(x, function(s) !is.null(s$cssh), logical(1))
    if (!all(have_cssh)) {
        cssh <- make_cssh(x)
        for (i in seq_along(x)) x[[i]]$cssh <- cssh
    } else {
        cssh <- x[[1]]$cssh
        for (i in seq_along(x)) {
            if (!isTRUE(all.equal(x[[i]]$cssh, cssh))) {
                stop(
                    "Spectra carry different cssh; alignment requires a shared grid."
                )
            }
        }
    }
    for (i in seq_along(x)) {
        if (is.null(x[[i]]$lcpar$pcide)) {
            x[[i]]$lcpar$pcide <- pci_on_cssh(x[[i]]$lcpar$x0, cssh)
        }
    }
    x
}

# Build the shared chemical-shift grid used by alignment.
#
# Range: intersection of every input spectrum's `cs` range. A peak
# fitted outside the intersection is excluded from the shared
# superposition `sit$supsh` (its Lorentzian still evaluates, but
# alignment only sees the part inside the intersection).
#
# Length: max(length(cs)) across input spectra. All inputs typically
# share a length already (alignment downstream enforces it), but if
# they differ, picking the max preserves the finest available
# resolution. Spacing is uniform and decreasing to match `cs`.
make_cssh <- function(x) {
    mins <- vapply(x, function(s) min(s$cs), numeric(1))
    maxs <- vapply(x, function(s) max(s$cs), numeric(1))
    ns_pts <- vapply(x, function(s) length(s$cs), integer(1))
    lo <- max(mins); hi <- min(maxs); n <- max(ns_pts)
    if (lo >= hi) stop(
        "Spectra cs ranges have empty intersection; cannot build cssh."
    )
    seq(hi, lo, length.out=n)
}

# Compute the shared-grid input vector sit$supsh for CluPA. Four modes:
#
#   "triangle" (default): narrow isoceles triangle per peak, full-width
#     = lambda, peak height = A. Round-first semantics — supsh peak
#     position is byte-equal to the integer peakList index passed into
#     hclust_align(). Cheap: O(npeaks * w_dp).
#
#   "rectangle": constant A on [c - hw, c + hw] per peak, zero outside.
#     Fastest. Discontinuous edges; can ring under FFT.
#
#   "lorentz": full Lorentz superposition evaluated at cssh. Smooth,
#     expensive: O(length(cssh) * npeaks).
#
#   "sparse": delta comb — zeros everywhere except at the cssh column
#     nearest each peak center, which carries the peak area A. Cheap:
#     O(npeaks). Lets speaq's FFT cross-correlator align directly on
#     amplitude-weighted peak positions.
#
# Only fills in sit$supsh when it is missing. Callers that want to
# *switch* the cached mode must clear sit$supsh first. Assumes
# ensure_cssh() has already run.
ensure_supsh <- function(x, supsh="triangle") {
    supsh <- match.arg(
        supsh, c("triangle", "rectangle", "lorentz", "sparse", "eiffel")
    )
    for (i in seq_along(x)) {
        if (is.null(x[[i]]$sit$supsh)) {
            x[[i]]$sit$supsh <- make_supsh(x[[i]]$cssh, x[[i]]$lcpar, supsh)
        }
    }
    x
}

# Build the sit$supsh vector from a peak list. Branches on `supsh`
# mode; see ensure_supsh() for semantics. Shared between clupa()'s
# ensure_supsh() entry point and build_clupa_consensus().
make_supsh <- function(cssh, lcpar, supsh="triangle") {
    supsh <- match.arg(
        supsh, c("triangle", "rectangle", "lorentz", "sparse", "eiffel")
    )
    if (supsh == "triangle")  return(triangle_sup(cssh, lcpar))
    if (supsh == "rectangle") return(rect_sup(cssh, lcpar))
    if (supsh == "lorentz")   return(lorentz_sup(cssh, lcpar=lcpar))
    if (supsh == "eiffel")    return(eiffel_sup(cssh, lcpar))
    # sparse
    out <- numeric(length(cssh))
    if (nrow(lcpar) == 0L) return(out)
    idx <- pci_on_cssh(lcpar$x0, cssh)
    s <- tapply(as.numeric(lcpar$A), idx, sum)
    out[as.integer(names(s))] <- s
    out
}

# Narrow-triangle superposition on `cssh`. Each peak contributes a
# symmetric triangle centered on the nearest cssh column to `x0`, with
# full base width = `lambda` (in datapoints) and peak height = `A`.
# Returns a numeric vector of length `length(cssh)`.
triangle_sup <- function(cssh, lcpar) {
    n <- length(cssh)
    if (nrow(lcpar) == 0L) return(numeric(n))
    pc <- pci_on_cssh(lcpar$x0, cssh)
    hw <- lambda_to_hw_dp(lcpar$lambda, cssh)
    A  <- as.numeric(lcpar$A)
    .Call(triangle_sup_c, as.integer(pc), as.integer(hw), A, n)
}

# Eiffel-tower superposition on `cssh`. Each peak contributes a wide
# flat triangle (full base = 4 * lambda, apex height = 0.25 * A) PLUS
# a single-column spike at the exact peak center of additional height
# `A` (so the total amplitude at the center column is 1.25 * A). The
# wide base extends the FFT cross-correlator's basin of attraction
# (peaks now reach each other from up to 2 * lambda away vs. 0.5 *
# lambda for plain triangle), while the spike snaps in once the
# shift is within one column of the true center.
eiffel_sup <- function(cssh, lcpar) {
    n <- length(cssh)
    if (nrow(lcpar) == 0L) return(numeric(n))
    pc <- pci_on_cssh(lcpar$x0, cssh)
    A <- as.numeric(lcpar$A)
    hw_wide <- lambda_to_hw_dp(lcpar$lambda * 4, cssh)
    out <- .Call(triangle_sup_c, as.integer(pc), as.integer(hw_wide),
                 A * 0.25, n)
    ok <- pc >= 1L & pc <= n
    if (any(ok)) {
        spike <- tapply(A[ok], pc[ok], sum)
        out[as.integer(names(spike))] <-
            out[as.integer(names(spike))] + spike
    }
    out
}

# Rectangle superposition on `cssh`. Each peak contributes a constant
# `A` over the cssh columns within `lambda/2` datapoints of the nearest
# cssh column to `x0`.
rect_sup <- function(cssh, lcpar) {
    n <- length(cssh)
    if (nrow(lcpar) == 0L) return(numeric(n))
    pc <- pci_on_cssh(lcpar$x0, cssh)
    hw <- lambda_to_hw_dp(lcpar$lambda, cssh)
    A  <- as.numeric(lcpar$A)
    .Call(rect_sup_c, as.integer(pc), as.integer(hw), A, n)
}

# Convert a vector of per-peak lambdas (ppm) into integer half-widths
# in cssh datapoints. Floor 1 — even a sub-datapoint lambda gets a
# single column on each side so the supsh peak is always >=3 columns
# wide (which keeps the FFT cross-correlator's job well-conditioned).
lambda_to_hw_dp <- function(lambda, cssh) {
    n <- length(cssh)
    w_dp <- abs(convert_width(abs(as.numeric(lambda)), cssh, seq_len(n)))
    pmax(1L, as.integer(round(w_dp / 2)))
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
    attr(x, "ref") <- ref
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

align_decon <- function(x, ref, maxShift, full=TRUE, use_speaq=FALSE,
                        supsh="triangle", shift_method="slide") {
    cssh <- x$cssh
    pci_x <- lcpar_pci(x$lcpar, cssh)
    pci_ref <- lcpar_pci(ref$lcpar, cssh)
    np_x <- length(pci_x); np_ref <- length(pci_ref)
    # `orig_idx` mirrors `peakList`: 0 for ref entries, 1..np_x for the
    # target peak's row in x$lcpar. Threaded through hclust_align so a
    # rebuild_tar call can look up the original A / lambda for each
    # surviving target peak after any shift.
    orig_idx <- c(integer(np_ref), seq_len(np_x))
    rebuild_tar <- if (use_speaq || shift_method == "slide") NULL else
        make_rebuild_tar(x$lcpar, cssh, supsh)
    obj <- hclust_align(
        refSpec=ref$sit$supsh, tarSpec=x$sit$supsh,
        peakList=c(pci_ref, pci_x),
        peakLabel=c(rep(1, np_ref), rep(0, np_x)),
        orig_idx=orig_idx,
        startP=1, endP=length(x$sit$supsh),
        maxShift=maxShift, use_speaq=use_speaq,
        shift_method=shift_method, rebuild_tar=rebuild_tar
    )
    if (length(obj$peakList) != np_ref + np_x) stop("Lost peaks during alignment")
    pcial <- obj$peakList[(np_ref + 1):(np_ref + np_x)]
    x$lcpar$x0al <- cssh[pcial]
    x$lcpar$pcial <- pcial
    if (full) x$sit$supal <- lorentz_sup(cssh, x$lcpar$x0al, x$lcpar$A, x$lcpar$lambda)
    class(x) <- c("align", "decon2", "spectrum")
    x
}

# Closure factory used by align_decon under shift_method="recompute".
# Returns a function rebuild_tar(orig_tar_idx, new_pcis) that yields
# the full-length supsh vector reconstructed from the original A /
# lambda (looked up by orig index) placed at the *current* integer
# positions on `cssh`. The shape mode (`supsh`) is captured.
#
# Caller writes the relevant slice of the returned vector back into
# `tarSpec`; the unused parts cost nothing for sparse shapes (their
# support is zero outside each peak's small window).
make_rebuild_tar <- function(lcpar_full, cssh, supsh) {
    A_all   <- as.numeric(lcpar_full$A)
    lam_all <- as.numeric(lcpar_full$lambda)
    function(tar_orig, tar_pcis) {
        if (length(tar_orig) == 0L) return(numeric(length(cssh)))
        lc <- data.frame(
            x0=cssh[tar_pcis],
            A=A_all[tar_orig],
            lambda=lam_all[tar_orig]
        )
        make_supsh(cssh, lc, supsh)
    }
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
    refSpec, tarSpec, peakList, peakLabel, startP, endP, maxShift,
    use_speaq = FALSE, orig_idx = NULL,
    shift_method = "slide", rebuild_tar = NULL
) {

    if (use_speaq) return(
        speaq::hClustAlign(
            refSpec=refSpec, tarSpec=tarSpec, peakList=peakList,
            peakLabel=peakLabel, startP=startP, endP=endP,
            distanceMethod="average", maxShift=maxShift, acceptLostPeak=FALSE
        )
    )

    if (is.null(orig_idx)) orig_idx <- integer(length(peakList))
    do_rebuild <- shift_method == "recompute" && !is.null(rebuild_tar)

    minPk <- min(peakList)
    maxPk <- max(peakList)

    # Narrow the active region to signal boundaries
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
                orig_idx <- orig_idx[-lost]
            }
            if (do_rebuild) {
                # Peak-list-driven rebuild: derive the supsh from
                # current peak positions instead of sliding the old
                # vector and edge-padding the vacated columns. Write
                # the slice into tarSpec[startP:endP] — HC partitions
                # peaks into disjoint position clusters, so this
                # recursion's segment never overlaps a sibling's.
                tar_now <- which(peakLabel == 0L)
                new_full <- rebuild_tar(orig_idx[tar_now], peakList[tar_now])
                tarSpec[startP:endP] <- new_full[startP:endP]
            } else {
                seg <- tarSpec[startCheckP:endCheckP]
                tarSpec[startCheckP:endCheckP] <- do_shift(seg, adj$stepAdj)
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
    sid1 <- orig_idx[left_set]
    id1 <- left_set
    sub2 <- peakList[right_set]
    lab2 <- peakLabel[right_set]
    sid2 <- orig_idx[right_set]
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
        sid1 <- orig_idx[left_set]
        id1 <- left_set
        sub2 <- peakList[right_set]
        lab2 <- peakLabel[right_set]
        sid2 <- orig_idx[right_set]
        id2 <- right_set
        max1 <- max(sub1); min2 <- min(sub2)
        endP1 <- max1 + which.min(tarSpec[(max1 + 1L):(min2 - 1L)])
        if (is.na(endP1) || endP1 > length(tarSpec)) endP1 <- max1
        startP2 <- endP1 + 1L
    }
    if (length(unique(lab1)) > 1L) {
        res <- hclust_align(refSpec, tarSpec, sub1, lab1, startP, endP1,
                            maxShift, orig_idx=sid1,
                            shift_method=shift_method,
                            rebuild_tar=rebuild_tar)
        tarSpec <- res$tarSpec
        peakList[id1] <- pad_peaks(res$peakList, length(id1))
    }
    if (length(unique(lab2)) > 1L) {
        res <- hclust_align(refSpec, tarSpec, sub2, lab2, startP2, endP,
                            maxShift, orig_idx=sid2,
                            shift_method=shift_method,
                            rebuild_tar=rebuild_tar)
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
#' reference's `cssh` column index as their `pcisn`; unmatched peaks get
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
    x <- ensure_cssh(x)
    if (!is.null(ref) && is.null(ref$cssh)) ref$cssh <- x[[1]]$cssh
    ref <- ref %||% find_ref(x)
    if (is.null(ref$lcpar$pcide)) {
        ref$lcpar$pcide <- pci_on_cssh(ref$lcpar$x0, ref$cssh)
    }
    cssh <- x[[1]]$cssh
    ro   <- order(ref$lcpar$x0)
    rx0  <- as.numeric(ref$lcpar$x0[ro])
    rcol <- as.integer(ref$lcpar$pcide[ro])
    rA   <- if (w_A > 0) normalize_A(as.numeric(ref$lcpar$A[ro])) else NULL
    for (s in seq_along(x)) {
        x[[s]]$lcpar <- snap_nw_lcpar(x[[s]]$lcpar, rx0, rcol, cssh, gap_tol,
                                       pos_field=pos_field, w_A=w_A, rA=rA)
        x[[s]]$sit$supal <- NULL
        class(x[[s]]) <- c("align", "decon2", "spectrum")
    }
    class(x) <- c("aligns", "decons2", "spectra")
    x
}

# Per-spectrum amplitude normalization: divide by median of the
# strictly-positive entries so the resulting vector has typical magnitude
# 1, scale-invariant across spectra of different total intensity. Used to
# put target and reference amplitudes on a comparable scale before
# computing log-ratio cost terms.
normalize_A <- function(A) {
    if (length(A) == 0L) return(A)
    pos <- A[A > 0]
    if (length(pos) == 0L) return(A)
    A / stats::median(pos)
}

# Per-spectrum NW snap: returns lcpar with pcisn / x0sn columns added.
# When `w_A > 0` the cost matrix mixes in a log-ratio amplitude term:
#   M = |Δx0| + w_A * gap_tol * |log(A_t / A_r)|
# so the gap-vs-match decision favors pairings of similarly-sized peaks.
# Both A vectors are pre-normalized to unit median so the ratio is
# scale-invariant across spectra.
snap_nw_lcpar <- function(lcpar, rx0, rcol, cssh, gap_tol, pos_field="x0",
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
        lcpar$x0sn[si]  <- cssh[ci]
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
#' (`list(cssh, lcpar, ...)`) so it can be passed directly as `ref` to
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
#' @return A list with components `cssh`, `lcpar` (data frame with `x0`,
#'   `A`, `lambda`, `pcide`), plus the class attributes
#'   `c("consensus", "align", "decon2", "spectrum")` so it behaves like a
#'   single-spectrum reference for downstream code.
#'
#' @author 2026 Tobias Schmidt: initial version.
build_consensus <- function(x, y=NULL, gap_tol=0.02, pos_field="x0") {
    stopifnot(inherits(x, "decons2"), is_num(gap_tol, 1), gap_tol > 0)
    x <- ensure_cssh(x)
    cssh <- x[[1]]$cssh

    # Pull positions / A / lambda from the chosen field. find_ref always
    # operates on raw x0 via pcide; for pos_field != "x0" we temporarily
    # swap pos_field into x0 before the find_ref pick (and snap call below).
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
        seed <- list(cssh=cssh, lcpar=seed_lcpar)
    }
    if (is.null(seed$lcpar$pcide)) {
        seed$lcpar$pcide <- pci_on_cssh(seed$lcpar$x0, cssh)
    }

    # Now extend the seed with any training peaks that didn't match it. We
    # snap on the same pos_field used to build the seed.
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
    full_lcpar$pcide <- pci_on_cssh(full_lcpar$x0, cssh)
    structure(
        list(cssh=cssh, lcpar=full_lcpar),
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
    # Cluster by consecutive-gap threshold.
    g <- c(0, cumsum(diff(x0) > gap_tol))
    keep_x0  <- vapply(split(seq_len(n), g), function(i) {
        w <- A[i]; if (sum(w) <= 0) mean(x0[i]) else stats::weighted.mean(x0[i], w)
    }, numeric(1))
    keep_A   <- vapply(split(seq_len(n), g), function(i) mean(A[i]),   numeric(1))
    keep_lam <- vapply(split(seq_len(n), g), function(i) mean(lam[i]), numeric(1))
    data.frame(x0=as.numeric(keep_x0), A=as.numeric(keep_A),
               lambda=as.numeric(keep_lam))
}
