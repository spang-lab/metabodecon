
# API #####

#' @export
#' @title Signal Integral Matrix
#'
#' @description
#' Extracts a peak-area matrix from aligned spectra. Rows are spectra,
#' columns are chemical-shift positions (features).
#'
#' @details
#' Rules:
#' 1. If `maxCombine==0` the peak areas are returned as is.
#' 2. If `maxCombine > 0 && length(peakPos) == 0` non-overlapping neighboring
#'    columns with `maxCombine` distance are merged first.
#' 3. If `maxCombine > 0 && length(peakPos) > 0` each peak is moved towards its
#'    closest `peakPos`, if it is within `maxCombine` distance. If multiple peaks
#'    have the same "closest peakPos", only the closest one is shifted.
#'
#' Example: assume the following Peak Area Matrix
#'
#'         1  2  3  4  5  6  7  8  9
#' si1 = c(0, 0, 0, 2, 0, 0, 0, 4, 0) --> pcial=c(4,8),   A=c(2,4)
#' si2 = c(0, 0, 3, 0, 0, 0, 0, 4, 0) --> pcial=c(3,8),   A=c(3,4)
#' si3 = c(0, 0, 2, 4, 0, 0, 5, 0, 0) --> pcial=c(3,4,7), A=c(2,4,5)
#' si4 = c(0, 0, 0, 0, 3, 0, 0, 0, 3) --> pcial=c(5,9),   A=c(3,3)
#' si5 = c(0, 0, 0, 0, 2, 0, 0, 0, 3) --> pcial=c(5,9),   A=c(2,3)
#'
#' If we call `si_mat(x, maxCombine=0)`, the matrix is returned as is.
#'
#' If we call `si_mat(x, maxCombine=1)`, we get
#'
#'         1  2  3  4  5  6  7  8  9
#' si1 = c(0, 0, 0, 2, 0, 0, 0, 4, 0)
#' si2 = c(0, 0, 3, 0, 0, 0, 0, 4, 0)
#' si3 = c(0, 0, 2, 4, 0, 0, 0, 5, 0)
#' si4 = c(0, 0, 0, 3, 0, 0, 0, 3, 0)
#' si5 = c(0, 0, 0, 2, 0, 0, 0, 3, 0)
#'
#' If we call `si_mat(x, maxCombine=1, peakPos=c(3,4,9))`, we get
#'
#'         1  2  3  4  5  6  7  8  9
#' si1 = c(0, 0, 0, 2, 0, 0, 0, 0, 4)
#' si2 = c(0, 0, 3, 0, 0, 0, 0, 0, 4)
#' si3 = c(0, 0, 2, 4, 0, 0, 5, 0, 0)
#' si4 = c(0, 0, 0, 3, 0, 0, 0, 0, 3)
#' si5 = c(0, 0, 0, 2, 0, 0, 0, 0, 3)
#'
#' @param x
#' An object of type `decons2` or `aligns`. For `aligns`, the aligned
#' peak positions are used; for `decons2`, the deconvoluted peak positions
#' are used.
#'
#' @param maxCombine
#' How many adjacent columns to consider for merging.
#'
#' @param peakPos
#' Integer vector of column indices in the `cs` grid to snap
#' peaks to. Used to align new spectra to the same features as a reference
#' matrix (see 'Details').
#'
#' @param igrs
#' List of length-2 numeric vectors `c(left, right)` (in ppm) of
#' ignore-regions. Columns whose chemical shift falls inside any region are
#' zeroed out. Use `list()` to disable.
#'
#' @param drop_zero
#' Drop columns where all values are zero?
#'
#' @return
#' A matrix with spectra in rows and chemical shifts as colnames.
#' Always has `length(x[[1]]$cs)` columns, regardless of `peakPos`.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'     decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#'     aligns <- align(decons)
#'     X0 <- si_mat(aligns)
#'     Xc <- si_mat(aligns, maxCombine = 20)
#'     # Reuse the feature grid of Xc on new spectra:
#'     pp <- which(colSums(Xc != 0) > 0)
#'     Xn <- si_mat(aligns, maxCombine = 20, peakPos = pp)
#' }
si_mat <- function(x, drop_zero=FALSE, maxCombine=0, peakPos=NULL, igrs=list()) {
    stopifnot(inherits(x, "decons2"))
    cs <- x[[1]]$cs
    ns <- length(x)
    nc <- length(cs)
    get_idx <- function(vals) {
        idx <- match(vals, cs)
        if (anyNA(idx)) idx <- round(convert_pos(vals, cs, seq_along(cs)))
        pmin(nc, pmax(1L, as.integer(idx)))
    }
    # Use aligned peak centers when available, else raw deconvolution centers.
    pos <- lapply(x, function(xi) xi$lcpar$x0al %||% xi$lcpar$x0)
    if (maxCombine == 0) {
        mat <- t(sapply(seq_len(ns), function(s) {
            al <- numeric(nc)
            idx <- x[[s]]$lcpar$pcial %||% get_idx(pos[[s]])
            al[idx] <- x[[s]]$lcpar$A * base::pi
            al
        }))
    } else {
        mat <- matrix(0, nrow = ns, ncol = nc)
        for (s in seq_len(ns)) {
            pidx <- get_idx(pos[[s]])
            A <- x[[s]]$lcpar$A
            for (p in seq_along(pidx)) {
                mat[s, pidx[p]] <- mat[s, pidx[p]] + A[p] * base::pi
            }
        }
        if (length(peakPos) == 0) {
            mat <- combine_peaks(mat, maxCombine)
        } else {
            mat <- snap_to_peakPos(mat, peakPos, maxCombine)
        }
    }
    if (length(igrs) > 0) {
        ig_mask <- logical(nc)
        for (r in igrs) ig_mask <- ig_mask | (cs >= min(r) & cs <= max(r))
        if (any(ig_mask)) mat[, ig_mask] <- 0
    }
    colnames(mat) <- cs
    rownames(mat) <- get_names(x)
    if (drop_zero) mat <- mat[, colSums(mat != 0) > 0, drop = FALSE]
    mat
}

#' @export
#' @title Extract Matrix of aligned Signal Intensities
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Deprecated in favour of [metabodecon::si_mat()], which returns the same data with
#' spectra in rows and features (chemical shifts) in columns.
#'
#' @inheritParams si_mat
#'
#' @return
#' A numeric matrix with chemical shifts as rownames and spectrum names as
#' colnames (the transpose of [metabodecon::si_mat()]).
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
get_si_mat <- function(x, drop_zero = FALSE, maxCombine = 0, peakPos = NULL) {
    lifecycle::deprecate_warn("2.0.0", "get_si_mat()", "si_mat()")
    t(si_mat(x, drop_zero = drop_zero, maxCombine = maxCombine, peakPos = peakPos))
}

#' @export
#' @title Peak-snapped feature matrix
#'
#' @description
#' Builds a feature matrix by snapping each spectrum's aligned peak
#' centers to the peak grid of a reference spectrum. When `peakPos` is
#' `NULL`, the reference is chosen via [metabodecon::find_ref()] and its
#' `lcpar$pcial` is used. Each peak is moved to its closest reference
#' position within `maxCombine` columns; collisions are summed and peaks
#' farther than `maxCombine` from any reference position are dropped.
#'
#' Equivalent to
#' `si_mat(x, maxCombine=maxCombine, peakPos=peakPos, igrs=igrs)` with
#' `peakPos` defaulted to the reference's `pcial`. Suitable as the
#' `feat_mat` argument of [metabodecon::fit_mdm()] and the recommended
#' default for the decon -> align -> classify pipeline.
#'
#' @param x
#' An `aligns` object (or `decons2`).
#'
#' @param maxCombine
#' Maximum allowed snap distance in chemical-shift columns.
#'
#' @param peakPos
#' Optional integer vector of column indices to snap peaks to.
#' Defaults to the reference spectrum's peak grid.
#'
#' @param igrs
#' List of two-element ppm intervals to ignore.
#'
#' @return
#' A numeric matrix with spectra in rows and chemical shifts as
#' colnames. Always has `length(x[[1]]$cs)` columns.
#'
#' @author 2024-2026 Tobias Schmidt: initial version.
peak_mat <- function(x, maxCombine=20, peakPos=NULL, igrs=list()) {
    stopifnot(inherits(x, "decons2"))
    if (is.null(peakPos)) {
        ref <- find_ref(x)
        peakPos <- ref$lcpar$pcial
    }
    si_mat(x, maxCombine=maxCombine, peakPos=peakPos, igrs=igrs)
}

#' @export
#' @title Bin a spectra-like object into a feature matrix
#' @description
#' Bins the per-spectrum signal vector left-to-right into chunks of
#' `maxCombine` chemical-shift columns and returns the per-bin sums as a
#' feature matrix. Columns whose chemical-shift falls inside any `igrs`
#' interval are removed before binning.
#'
#' Accepts three input types:
#' \itemize{
#'   \item `spectra`: uses `x[[i]]$si` directly.
#'   \item `decons2`: uses `x[[i]]$sit$sup` (smoothed reconstruction).
#'   \item `aligns`: builds a sparse vector from `lcpar$pcial` /
#'         `lcpar$A * pi`, then bins.
#' }
#'
#' Suitable as the `feat_mat` argument of [metabodecon::fit_mdm()].
#'
#' @param x A `spectra`, `decons2` or `aligns` object.
#' @param maxCombine Bin width in chemical-shift columns.
#' @param peakPos Ignored. Accepted for API compatibility with
#'   [metabodecon::si_mat()].
#' @param igrs List of two-element ppm intervals to ignore.
#'
#' @return A numeric matrix with one row per spectrum and one column per
#'   bin.
bin <- function(x, maxCombine=128, peakPos=NULL, igrs=list()) {
    stopifnot(
        inherits(x, "spectra") || inherits(x, "decons2") ||
            inherits(x, "aligns"),
        is_int(maxCombine, 1), maxCombine >= 1
    )
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
            idx <- x[[s]]$lcpar$pcial
            if (length(idx) > 0) Y[s, idx] <- x[[s]]$lcpar$A * base::pi
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

# Combine Peaks #####

#' @noRd
#'
#' @title Combine Peaks
#'
#' @description
#' Even after calling [metabodecon::speaq_align()], the alignment of individual
#' signals is not always perfect, as 'speaq' performs a segment-wise alignment
#' i.e. groups of signals are aligned. For further improvements, partly filled
#' neighboring columns are merged. See 'Details' for an illustrative example.
#'
#' @param M A matrix of signal integrals (spectra × datapoints).
#' @param maxCombine Amount of adjacent columns to consider for merging.
#' @param lower_bound Required amount of non-zero elements to trigger merging.
#'
#' @return
#' A data frame containing one column for each data point in the original spectrum. The second
#' data frame contains only columns where at least one entry is non-zero.
#'
#' @details
#'
#' Example of what the function does:
#'
#' ```txt
#' |            | 1    | 2    | 3    | 4    | 5    |
#' |----------- |------|------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0    | 0.11 | 0    |
#' | Spectrum 2 | 0    | 0.88 | 0    | 0.12 | 0    |
#' | Spectrum 3 | 0.07 | 0.56 | 0.30 | 0    | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0.07 | 0    | 0.07 |
#' | Spectrum 5 | 0.04 | 0    | 0    | 0.04 | 0    |
#' ```
#'
#' becomes
#'
#' ```txt
#' |            | 1    | 2    | 3    | 4    | 5    |
#' |----------- |------|------|------|------|------|
#' | Spectrum 1 | 0.13 | 0    | 0    | 0.11 | 0    |
#' | Spectrum 2 | 0    | 0.88 | 0    | 0.12 | 0    |
#' | Spectrum 3 | 0.07 | 0.56 | 0    | 0.30 | 0    |
#' | Spectrum 4 | 0.08 | 0    | 0    | 0.07 | 0.07 |
#' | Spectrum 5 | 0.04 | 0    | 0    | 0.04 | 0    |
#' ```
#'
#' I.e.
#'
#' 1. Column 1 and 2 get NOT merged, because they have a common non-zero entry.
#' 2. Column 3 and 4 get merged, because they are in `range` of each other and
#'    have no common non-zero entries.
#' 3. Column 4 and 5 get NOT merged, because it is more beneficial to merge
#'    column 3 and 4, as they have more mergeable entries and after merging
#'    column 3 and 4, column 4 and 5 have a common non-zero entry.
#'
#' @author
#'
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' deps <- c("MassSpecWavelet", "impute")
#' deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
#' if (all(deps_installed)) {
#'     # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
#'     sim_subset <- metabodecon_file("bruker/sim_subset")
#'     spectrum_data <- generate_lorentz_curves_sim(sim_subset)
#'     shifted_mat <- speaq_align(spectrum_data = spectrum_data, verbose = FALSE)
#'     range <- 5
#'     lower_bound <- 1
#'     obj <- combine_peaks(shifted_mat, range, lower_bound)
#'     str(obj)
#' }
combine_peaks <- function(M, maxCombine=5, lower_bound=1) {
    if (anyNA(M)) stop("`M` must not contain NA values.", call.=FALSE)
    U <- M != 0
    uu <- colSums(U)
    nc <- ncol(M)
    for (i in (nrow(M) - 1):lower_bound) {
        for (j in which(uu == i)) {
            if (uu[j] == 0) next
            nn <- seq(max(1, j - maxCombine), min(nc, j + maxCombine))
            nn <- nn[nn != j]
            if (length(nn) == 0) next
            mj <- M[, j]
            uj <- U[, j]
            repeat {
                nn <- nn[uu[nn] > 0]
                if (length(nn) == 0) break
                cc <- combine_scores(U, uu, j, nn, uj = uj)
                if (max(cc) == 0) break
                n <- nn[which.max(cc)]
                mj <- mj + M[, n]
                uj <- uj | U[, n]
                uu[j] <- uu[j] + uu[n]
                M[, n] <- 0
                U[, n] <- FALSE
                uu[n] <- 0
                nn <- nn[nn != n]
                if (length(nn) == 0) break
            }
            M[, j] <- mj
            U[, j] <- uj
        }
    }
    M
}

#' @noRd
#'
#' @description
#' Calculates a "combine score" for a set of columns `nn` of a Matrix `M`. The
#' score describes how "beneficial" it is to merge column `n` (from `nn`) into
#' column `j`. A column `n` is considered "combinable" if there is no row where
#' columns `n` and `j` both have a non-zero element. If a column is not
#' combinable, its combine score is 0. If a column is combinable, its combine
#' score is the amount of non-zero elements in the column. This function
#' calculates the combine score for all neighbour column `nn` and then returns
#' the calculated combine scores as vector.
#'
#' @param U Matrix describing nonzero entries of `M.` `U[i,j]` should be `TRUE`
#' if `M[i,j]` is nonzero, else `FALSE.` I.e. for any matrix `M` you can
#' generate `U` as `U <- M != 0`.
#'
#' @param uu Vector describing the amount of nonzero elements in each column of
#' `M`. I.e. for any matrix `M` you can generate `uu` as `uu <- colSums(M !=
#' 0)`. Called `uu` because it denotes the amount of elements that are *unequal
#' zero*.
#'
#' @param j Index of the column where the neighbor columns `nn` should be
#' combined into.
#'
#' @param nn Indices of the columns to calculate the combine score for. Called
#' `nn` because this should be a set of *neighbor columns*.
#'
#' @details
#' Since we only need to know whether an element is nonzero in order to
#' calculate the combine scores, the function takes the matrix `U` and the
#' vector `uu` as input instead of their common ancestor `M`. This is much more
#' efficient if the function is called multiple times, as the conversion must
#' not be done multiple times.
#'
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
#'
#' @examples
#' M <- rbind(
#'     c(2, 0, 2, 2, 0),
#'     c(2, 1, 0, 2, 0),
#'     c(0, 1, 0, 2, 0),
#'     c(0, 0, 3, 0, 1)
#' )
#' U <- M != 0
#' uu <- colSums(U) # 2 2 2 3 1
#'
#' cc <- combine_scores(U, uu, j = 2, nn = c(1, 3))
#' cc[1] == 0  # M[,1] and M[, 2] are not combinable
#' cc[2] == 2  # M[,3] and M[, 2] are combinable and M[, 3] has two nonzero elements
#'
#' cc <- combine_scores(U, uu, j = 1, nn = 2:5)
#' cc[1] == 0  # M[, 2] and M[, 1] are not combinable
#' cc[2] == 0  # M[, 3] and M[, 1] are not combinable
#' cc[3] == 0  # M[, 4] and M[, 1] are not combinable
#' cc[4] == 1  # M[, 5] and M[, 1] are combinable and M[, 5] has one nonzero element
combine_scores <- function(U, uu, j, nn, uj = NULL) {
    nn <- nn[nn >= 1 & nn <= ncol(U)]
    if (length(nn) == 0) return(numeric(0))
    if (is.null(uj)) uj <- U[, j]

    # A neighbor is combinable if it has no shared nonzero row with column j.
    overlaps <- .colSums(U[, nn, drop = FALSE] & uj, nrow(U), length(nn))
    cc <- uu[nn]
    cc[overlaps > 0] <- 0
    unname(cc)
}

#' @noRd
#' @title Snap peaks to a fixed feature grid
#'
#' @description
#' For each spectrum row of `mat` and each `peakPos`, the single closest
#' non-zero entry within `maxCombine` columns is moved to that `peakPos`
#' column in the output. All other non-zero entries are dropped. When two
#' peaks are equidistant from a `peakPos`, the leftmost one wins.
#'
#' @param mat Numeric matrix (spectra × cs grid) of raw peak areas.
#' @param peakPos Integer vector of target column indices.
#' @param maxCombine Maximum allowed snap distance in columns.
#'
#' @return A matrix with the same dimensions as `mat`.
snap_to_peakPos <- function(mat, peakPos, maxCombine) {
    ns <- nrow(mat)
    nc <- ncol(mat)
    pp <- sort(unique(pmin(nc, pmax(1L, as.integer(peakPos)))))
    out <- matrix(0, nrow = ns, ncol = nc)
    if (length(pp) == 0) return(out)
    for (s in seq_len(ns)) {
        nz <- which(mat[s, ] != 0)
        if (length(nz) == 0) next
        for (p in pp) {
            d <- nz - p                    # signed distance
            ad <- abs(d)
            ok <- which(ad <= maxCombine)
            if (length(ok) == 0) next
            # Among candidates, pick leftmost of those with minimum distance.
            min_d <- min(ad[ok])
            candidates <- ok[ad[ok] == min_d]
            winner <- candidates[which.min(nz[candidates])]
            out[s, p] <- mat[s, nz[winner]]
        }
    }
    out
}