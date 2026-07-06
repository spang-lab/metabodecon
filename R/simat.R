
# API #####

#' @export
#' @title Signal-integral matrix
#'
#' @description
#' Builds a per-spectrum peak-area matrix. Each row is a spectrum, each
#' column is a chemical-shift datapoint. For each peak, the column is
#' picked from `lcpar$pcisn` (post-snap) when available, else
#' `lcpar$pcial` (post-CluPA), else `lcpar$pcide` (post-decon). Peaks
#' with `pcisn = NA` (snapped beyond `maxCombine`) are skipped.
#' Collisions on the same column have their `A * pi` summed.
#'
#' `si_mat()` is intentionally a dumb peak-list rasterizer: all
#' alignment (continuous shift via CluPA) and snapping (RefPA) must
#' have happened upstream — typically inside [metabodecon::align()].
#' To build a feature matrix where every spectrum shares the same
#' column grid, run `align(x, maxShift, maxCombine)` first.
#'
#' @param x A `decons2` or `aligns` object.
#' @param drop_zero Drop columns whose entries are all zero?
#' @param igrs List of two-element ppm intervals to zero out before
#'   returning.
#' @param peakPos Optional integer column indices. When supplied (predict
#'   mode) the matrix is subset to those columns; when `NULL` and used as
#'   a `feat_fun` the non-zero columns are kept and attached as
#'   `attr(., "peakPos")`.
#' @param ... Ignored (protocol compatibility with other `feat_fun`s).
#'
#' @return A numeric matrix with one row per spectrum and
#'   `length(x[[1]]$cs)` columns (the full cs grid). Column names are
#'   ppm values; row names are spectrum names.
#'
#' @author 2024-2026 Tobias Schmidt: initial version.
#'
#' @examples
#' \dontrun{
#'   decons <- deconvolute(sim[1:2], sfr=c(3.55, 3.35))
#'   aligned <- align(decons, maxShift=50, maxCombine=20)
#'   X <- si_mat(aligned)
#' }
si_mat <- function(x, drop_zero=FALSE, igrs=list(), peakPos=NULL, ...) {
    stopifnot(inherits(x, "decons2"))
    feat_mode <- !missing(peakPos)
    # Every spectrum shares the same $cs (enforced upstream by
    # harmonize_grid + the grid-equality assertion inside clupa /
    # snap_to_ref), so the first spectrum's $cs is canonical.
    cs <- x[[1]]$cs
    ns <- length(x)
    nc <- length(cs)
    mat <- matrix(0, nrow=ns, ncol=nc)
    for (s in seq_len(ns)) {
        lcpar <- x[[s]]$lcpar
        if (nrow(lcpar) == 0L) next
        idx <- lcpar_idx(lcpar, cs)
        A <- lcpar$A
        # Sum collisions at the same column; skip NA (out-of-range snap).
        for (p in seq_along(idx)) {
            if (is.na(idx[p])) next
            mat[s, idx[p]] <- mat[s, idx[p]] + A[p] * base::pi
        }
    }
    if (length(igrs) > 0) {
        ig_mask <- logical(nc)
        for (r in igrs) ig_mask <- ig_mask | (cs >= min(r) & cs <= max(r))
        if (any(ig_mask)) mat[, ig_mask] <- 0
    }
    colnames(mat) <- cs
    rownames(mat) <- get_names(x)
    if (feat_mode) {
        if (is.null(peakPos)) peakPos <- which(colSums(mat != 0) > 0)
        mat <- mat[, peakPos, drop=FALSE]
        attr(mat, "peakPos") <- peakPos
        return(mat)
    }
    if (drop_zero) mat <- mat[, colSums(mat != 0) > 0, drop=FALSE]
    mat
}

# Pick the most-aligned peak-column index for each peak in `lcpar`.
# Returns integer indices into the shared `cs` grid, with NA for peaks
# that were snapped out (beyond `maxCombine`). Priority: `pcisn`
# (post-RefPA) > `pcial` (post-CluPA) > `pcide` (post-decon). Falls
# back to deriving the index from `x0al`/`x0` for legacy objects that
# pre-date the pci* fields.
lcpar_idx <- function(lcpar, cs) {
    pcisn <- lcpar[["pcisn"]]
    if (!is.null(pcisn)) return(as.integer(pcisn))
    pcial <- lcpar[["pcial"]]
    if (!is.null(pcial)) return(as.integer(pcial))
    pcide <- lcpar[["pcide"]]
    if (!is.null(pcide)) return(as.integer(pcide))
    pos <- lcpar[["x0al"]] %||% lcpar[["x0"]]
    nc <- length(cs)
    idx <- match(pos, cs)
    if (anyNA(idx)) idx <- round(convert_pos(pos, cs, seq_along(cs)))
    pmin(nc, pmax(1L, as.integer(idx)))
}

#' @export
#' @title Peak feature matrix
#'
#' @description
#' Thin wrapper around [metabodecon::si_mat()] suitable as the
#' `feat_fun` argument of [metabodecon::fit_mdm()]. Equivalent to
#' `si_mat(x, igrs=igrs)`; the snapping that used to live here has
#' moved into [metabodecon::align()] (RefPA stage).
#'
#' @param x An `aligns` object (or `decons2`).
#' @param igrs List of two-element ppm intervals to ignore.
#' @param peakPos Optional integer column indices, forwarded to
#'   [metabodecon::si_mat()] for predict mode.
#' @param ... Ignored. Accepted so `peak_mat` and [metabodecon::bin()]
#'   share a single `feat_fun(x, maxCombine, igrs)` protocol;
#'   `peak_mat` ignores `maxCombine` because snapping happens upstream
#'   inside [metabodecon::align()].
#'
#' @return A numeric matrix with spectra in rows and chemical shifts as
#'   colnames.
#'
#' @author 2024-2026 Tobias Schmidt: initial version.
peak_mat <- function(x, igrs=list(), peakPos=NULL, ...) {
    si_mat(x, igrs=igrs, peakPos=peakPos)
}

#' @export
#' @title Bin a spectra-like object into a feature matrix
#' @description
#' Bins the per-spectrum signal vector left-to-right into chunks of
#' `maxCombine` chemical-shift columns and returns the per-bin sums as
#' a feature matrix. Columns whose chemical-shift falls inside any
#' `igrs` interval are removed before binning.
#'
#' Accepts three input types:
#' \itemize{
#'   \item `spectra`: uses `x[[i]]$si` directly.
#'   \item `decons2`: uses `x[[i]]$sit$sup` (smoothed reconstruction).
#'   \item `aligns`: builds a sparse vector from `lcpar$pcial` /
#'         `lcpar$A * pi`, then bins.
#' }
#'
#' Suitable as the `feat_fun` argument of [metabodecon::fit_mdm()] for
#' binning baselines.
#'
#' @param x A `spectra`, `decons2` or `aligns` object.
#' @param maxCombine Bin width in chemical-shift columns.
#' @param igrs List of two-element ppm intervals to ignore.
#' @param peakPos Optional integer column indices for predict mode; when
#'   `NULL` the non-zero bins are kept and attached as `attr(., "peakPos")`.
#' @param ... Ignored (protocol compatibility with peak_mat).
#'
#' @return A numeric matrix with one row per spectrum and one column per
#'   bin.
bin <- function(x, maxCombine=128, igrs=list(), peakPos=NULL, ...) {
    stopifnot(
        inherits(x, "spectra") || inherits(x, "decons2") ||
            inherits(x, "aligns"),
        is_int(maxCombine, 1), maxCombine >= 1
    )
    # Every spectrum in `x` shares $cs (harmonize_grid + grid-equality
    # assertion inside the alignment stages). `aligns` carries `pcial`
    # as integer indices into that shared grid.
    cs <- x[[1]]$cs
    nc <- length(cs)
    ns <- length(x)
    keep <- rep(TRUE, nc)
    for (r in igrs) keep <- keep & !(cs >= min(r) & cs <= max(r))
    kept <- which(keep)
    if (length(kept) == 0) stop("All ppm range is ignored.", call.=FALSE)
    grp <- ceiling(seq_along(kept) / maxCombine)
    groups <- split(kept, grp)

    # Build per-spectrum signal vectors of length `nc`.
    Y <- matrix(0, nrow=ns, ncol=nc)
    if (inherits(x, "aligns")) {
        for (s in seq_len(ns)) {
            lcpar <- x[[s]]$lcpar
            if (nrow(lcpar) == 0L) next
            idx <- lcpar_idx(lcpar, cs)
            A <- lcpar$A
            ok <- !is.na(idx)
            for (p in which(ok)) {
                Y[s, idx[p]] <- Y[s, idx[p]] + A[p] * base::pi
            }
        }
    } else if (inherits(x, "decons2")) {
        for (s in seq_len(ns)) Y[s, ] <- x[[s]]$sit$sup
    } else {
        for (s in seq_len(ns)) Y[s, ] <- x[[s]]$si
    }

    nb <- length(groups)
    out <- matrix(0, nrow=ns, ncol=nb)
    centers <- numeric(nb)
    for (j in seq_len(nb)) {
        cols <- groups[[j]]
        out[, j] <- rowSums(Y[, cols, drop=FALSE])
        centers[j] <- mean(cs[cols])
    }
    colnames(out) <- sprintf("%.4f", centers)
    rownames(out) <- get_names(x)
    if (is.null(peakPos)) peakPos <- which(colSums(out != 0) > 0)
    out <- out[, peakPos, drop=FALSE]
    attr(out, "peakPos") <- peakPos
    out
}

#' @export
#' @title 700-bin Zacharias 2013 feature matrix
#'
#' @description
#' Builds a feature matrix on the fixed 700-bin grid of Zacharias
#' (2013): 300 bins covering 6.5-9.5 ppm + 400 bins covering 0.5-4.5
#' ppm, both at 0.01 ppm width. The water region (4.5-6.5 ppm) is
#' excluded. Bins are ordered high-to-low — column 1 covers
#' (9.49, 9.50) ppm, column 700 covers (0.50, 0.51) ppm.
#'
#' Suitable as the `feat_fun` argument of [metabodecon::fit_mdm()].
#' Per-spectrum dispatch:
#'
#' - If `lcpar` is empty (raw spectra): bin `$si` directly.
#' - If `lcpar` is non-empty (deconvoluted / aligned spectra):
#'   reconstruct `si_hat = lorentz_sup(cs, lcpar)` on the spectrum's
#'   own `cs` grid using `x0al` when present (else `x0`), then bin
#'   the reconstruction.
#'
#' Always returns all 700 columns; `maxCombine`, `igrs` and `peakPos`
#' are accepted for `feat_fun` protocol compatibility but ignored
#' (the bin layout is hardcoded).
#'
#' @param x A `spectra`, `decons2`, or `aligns` object.
#' @param maxCombine Ignored. Accepted for protocol compatibility.
#' @param igrs Ignored. Accepted for protocol compatibility.
#' @param peakPos Ignored. Accepted for protocol compatibility.
#' @param ... Ignored.
#'
#' @return A numeric matrix with one row per spectrum and 700 columns
#'   of bin sums. Column names are `sprintf("%.4f", bin_midpoint)`;
#'   row names are spectrum names.
#'
#' @author 2026 Tobias Schmidt: initial version.
bin700 <- function(x, maxCombine=0L, igrs=list(), peakPos=NULL, ...) {
    stopifnot(is_spectra(x))
    centers <- c(seq(9.495, length.out=300, by=-0.01), seq(4.495, length.out=400, by=-0.01))
    mat <- matrix(0, nrow=length(x), ncol=700L)
    for (i in seq_along(x)) mat[i, ] <- bin700_one(x[[i]])
    rownames(mat) <- get_names(x)
    colnames(mat) <- sprintf("%.4f", centers)
    mat
}

# Bin one spectrum onto the 700-bin grid. Uses lorentz reconstruction
# when lcpar carries fitted peaks, raw $si otherwise. Per-peak position
# priority: x0sn (snapped) > x0al (aligned) > x0 (raw). Falls back to
# the next column when x0sn is NA (snap_to_ref leaves NA for peaks
# beyond maxCombine), so partial snaps still reconstruct cleanly.
bin700_one <- function(s) {
    cs <- s$cs
    if (!is.null(s$lcpar) && nrow(s$lcpar) > 0L) {
        pos <- s$lcpar$x0
        if (!is.null(s$lcpar$x0al)) pos <- s$lcpar$x0al
        if (!is.null(s$lcpar$x0sn)) {
            ok <- !is.na(s$lcpar$x0sn)
            pos[ok] <- s$lcpar$x0sn[ok]
        }
        si <- lorentz_sup(cs, x0=pos, A=s$lcpar$A, lambda=s$lcpar$lambda)
    } else {
        si <- s$si
    }
    sums <- numeric(700L)
    sums[1:300] <- bin700_range(cs, si, 6.5, 9.5, 300L)
    sums[301:700] <- bin700_range(cs, si, 0.5, 4.5, 400L)
    sums
}

# Sum `si` into `nb` equal-width bins covering ppm interval (lo, hi),
# high-ppm-first. Mirrors the Zacharias 2013 indexing
# (floor((hi - cs) / dw) + 1). Uses a cumsum + run-boundary diff to
# avoid `tapply`/`factor` overhead — ~35x faster on a 7000-point
# range. Relies on `cs` being monotonic (which it is for any sane NMR
# spectrum); that makes `bi` non-decreasing so each bin's entries are
# contiguous and its sum is a single cumsum-difference.
bin700_range <- function(cs, si, lo, hi, nb) {
    idx <- which(cs > lo & cs < hi)
    if (length(idx) == 0L) return(numeric(nb))
    dw <- (hi - lo) / nb
    bi <- pmin(pmax(floor((hi - cs[idx]) / dw) + 1L, 1L), nb)
    csum <- c(0, cumsum(si[idx]))
    bnd <- c(0L, which(diff(bi) != 0L), length(bi))
    out <- numeric(nb)
    out[bi[bnd[-1L]]] <- diff(csum[bnd + 1L])
    out
}
