
# API #####

#' @export
#' @title Signal Integral Matrix
#'
#' @description
#' Extracts a peak-area matrix from aligned spectra. Rows are spectra,
#' columns are chemical-shift positions (features).
#'
#' - `maxCombine = 0`: returns the raw aligned integral vectors.
#' - `maxCombine > 0`, `intervals = NULL`: applies [combine_peaks()] to merge
#'   neighboring columns before returning.
#' - `maxCombine > 0`, `intervals` provided: returns one column per interval,
#'   mapping each spectrum's nearest in-interval peak to that column.
#'
#' @param x
#' An object of type `aligns`.
#'
#' @param maxCombine
#' Maximum snap distance in datapoints. `0` (default): raw integrals are
#' returned. `> 0` with `intervals = NULL`: applies `combine_peaks`. `> 0`
#' with `intervals`: maps peaks to interval centers.
#'
#' @param intervals
#' A data frame with integer columns `center`, `min`, and `max` (datapoint
#' indices). Each row defines one interval. A peak of spectrum `i` is assigned
#' to interval `j` if it falls within `[min[j], max[j]]` and is the closest
#' peak to `center[j]`; all others are discarded. When `NULL` (default) and
#' `maxCombine > 0`, `combine_peaks` is used instead.
#'
#' @param drop_zero
#' If `TRUE`, columns where all values are zero are removed.
#'
#' @return
#' A numeric matrix with spectrum names as rownames and chemical shifts
#' (or interval centers) as colnames.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'     decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#'     aligns <- align(decons)
#'     X0 <- si_mat(aligns)
#'     Xc <- si_mat(aligns, maxCombine = 20)
#' }
si_mat <- function(x, drop_zero = FALSE, maxCombine = 0, intervals = NULL) {
    stopifnot(is_aligns(x))
    cs <- x[[1]]$cs
    ns <- length(x)
    nc <- length(cs)
    if (maxCombine == 0) {
        mat <- t(sapply(x, function(xi) {
            al <- numeric(nc)
            al[xi$lcpar$cial] <- xi$lcpar$A * base::pi
            al
        }))
        colnames(mat) <- cs
    } else {
        get_idx <- function(vals) {
            idx <- match(vals, cs)
            if (anyNA(idx)) idx <- round(convert_pos(vals, cs, seq_along(cs)))
            pmin(nc, pmax(1L, as.integer(idx)))
        }
        mat <- matrix(0, nrow = ns, ncol = nc)
        for (s in seq_len(ns)) {
            pidx <- get_idx(x[[s]]$lcpar$x0al)
            A <- x[[s]]$lcpar$A
            for (p in seq_along(pidx)) mat[s, pidx[p]] <- mat[s, pidx[p]] + A[p] * base::pi
        }
        if (is.null(intervals)) {
            mat <- combine_peaks(mat, maxCombine)
            colnames(mat) <- cs
        } else {
            ni <- nrow(intervals)
            result <- matrix(0, nrow = ns, ncol = ni)
            for (s in seq_len(ns)) {
                row <- mat[s, ]
                nz <- which(row != 0)
                if (length(nz) == 0) next
                for (j in seq_len(ni)) {
                    lo <- intervals$min[j]; hi <- intervals$max[j]
                    cand <- nz[nz >= lo & nz <= hi]
                    if (length(cand) == 0) next
                    k <- cand[which.min(abs(cand - intervals$center[j]))]
                    result[s, j] <- row[k]
                }
            }
            mat <- result
            colnames(mat) <- cs[intervals$center]
        }
    }
    if (drop_zero) mat <- mat[, colSums(mat != 0) > 0, drop = FALSE]
    rownames(mat) <- get_names(x)
    mat
}

#' @export
#' @title Extract Matrix of aligned Signal Intensities
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Deprecated in favour of [si_mat()], which returns the same data with
#' spectra in rows and features (chemical shifts) in columns.
#'
#' @inheritParams si_mat
#'
#' @return
#' A numeric matrix with chemical shifts as rownames and spectrum names as
#' colnames (the transpose of [si_mat()]).
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
get_si_mat <- function(x, drop_zero = FALSE, maxCombine = 0, intervals = NULL) {
    lifecycle::deprecate_warn("2.0.0", "get_si_mat()", "si_mat()")
    t(si_mat(x, drop_zero = drop_zero, maxCombine = maxCombine,
             intervals = intervals))
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
# Build a data frame of non-overlapping intervals around `centers`.
# Each interval stretches `maxCombine` datapoints in each direction,
# clipped to [1, nc]. Where two intervals would overlap, both borders
# are shrunk to the midpoint (floor) of the two centers.
make_intervals <- function(centers, maxCombine, nc) {
    n <- length(centers)
    lo <- pmax(1L, centers - maxCombine)
    hi <- pmin(nc, centers + maxCombine)
    if (n >= 2) {
        mids <- floor((centers[-n] + centers[-1]) / 2)
        hi[-n] <- pmin(hi[-n], mids)
        lo[-1] <- pmax(lo[-1], mids)
    }
    data.frame(center = centers, min = lo, max = hi)
}

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


