
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
si_mat <- function(x, drop_zero=FALSE, igrs=list()) {
    stopifnot(inherits(x, "decons2"))
    # Use the shared chemical-shift grid (cssh) so peaks at the same
    # ppm in different spectra land in the same column even when the
    # input spectra had different cs ranges.
    cs <- x[[1]]$cssh %||% x[[1]]$cs
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
    if (drop_zero) mat <- mat[, colSums(mat != 0) > 0, drop=FALSE]
    mat
}

# Pick the most-aligned peak-column index for each peak in `lcpar`.
# Returns integer indices into `cs` (typically a cssh grid), with NA
# for peaks that were snapped out (beyond `maxCombine`). Priority:
# `pcisn` (post-RefPA) > `pcial` (post-CluPA) > `pcide` (post-decon).
# Falls back to deriving the index from `x0al`/`x0` for legacy objects
# that pre-date the pci* fields.
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
#' @title Extract matrix of aligned signal intensities
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Deprecated in favour of [metabodecon::si_mat()], which returns the
#' same data with spectra in rows and features (chemical shifts) in
#' columns.
#'
#' @inheritParams si_mat
#'
#' @return
#' A numeric matrix with chemical shifts as rownames and spectrum names
#' as colnames (the transpose of [metabodecon::si_mat()]).
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
get_si_mat <- function(x, drop_zero=FALSE) {
    lifecycle::deprecate_warn("2.0.0", "get_si_mat()", "si_mat()")
    t(si_mat(x, drop_zero=drop_zero))
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
#' @param ... Ignored. Accepted so `peak_mat` and [metabodecon::bin()]
#'   share a single `feat_fun(x, maxCombine, igrs)` protocol;
#'   `peak_mat` ignores `maxCombine` because snapping happens upstream
#'   inside [metabodecon::align()].
#'
#' @return A numeric matrix with spectra in rows and chemical shifts as
#'   colnames.
#'
#' @author 2024-2026 Tobias Schmidt: initial version.
peak_mat <- function(x, igrs=list(), ...) {
    si_mat(x, igrs=igrs)
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
#' @param ... Ignored (protocol compatibility with peak_mat).
#'
#' @return A numeric matrix with one row per spectrum and one column per
#'   bin.
bin <- function(x, maxCombine=128, igrs=list(), ...) {
    stopifnot(
        inherits(x, "spectra") || inherits(x, "decons2") ||
            inherits(x, "aligns"),
        is_int(maxCombine, 1), maxCombine >= 1
    )
    # `aligns` carries `pcial` on the shared cssh grid; raw spectra and
    # un-aligned decons2 still live on per-spectrum `cs` (and the legacy
    # assumption that `x[[1]]$cs` is representative).
    cs <- if (inherits(x, "aligns")) x[[1]]$cssh %||% x[[1]]$cs else x[[1]]$cs
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
    out
}
