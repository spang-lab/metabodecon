
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
#' An object of type `decons2` or `aligns`, as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @param ...
#' Reserved for internal package use. External callers should not pass
#' additional arguments.
#'
#' @param maxShift
#' Maximum number of datapoints a peak center may be shifted during CluPA
#' alignment. 50 is a suitable starting value for plasma spectra with a digital
#' resolution of 128K. Increase for urine or other sample types with larger
#' chemical-shift variation.
#'
#' @param maxCombine
#' Maximum distance in datapoints within which two adjacent peak columns may be
#' merged into one during the combine-peaks step. `0` (default) disables
#' merging. Only peaks from different spectra that have no overlap are eligible.
#'
#' @param verbose
#' Whether to print progress messages during alignment.
#'
#' @param nworkers
#' Number of parallel workers. Default is 1 (no parallelism).
#'
#' @param ref
#' Optional reference spectrum of type `align` or `decon2`. When supplied,
#' all spectra in `x` are aligned towards this reference. The reference is
#' prepended to `x` internally and removed from the result. If `NULL`
#' (default), the reference is chosen automatically.
#'
#'
#' @return
#' An object of type `aligns` as described in [Metabodecon
#' Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#' aligned <- align(decons)
align <- function(x,
                  maxShift = 50,
                  maxCombine = 0,
                  verbose = TRUE,
                  nworkers = 1,
                  ref = NULL,
                  ...) {
    opts  <- parse_align_dots(list(...), verbose)
    xx    <- chkargs_align(x, verbose, ref, opts$method)
    pciaa <- clupa_decons(xx, maxShift, ref, opts$method, verbose, nworkers)
    smat  <- discretize_pciaa(pciaa, xx)
    cmat  <- if (maxCombine >= 1) combine_peaks(smat, maxCombine) else smat
    build_aligns(xx, cmat, opts$full, !is.null(ref))
}

# Exported Helpers #####

#' @export
#' @title Extract Matrix of aligned Signal Intensities
#'
#' @description
#' Extracts a peak-area matrix from aligned spectra. Rows are
#' chemical-shift positions, columns are spectra. When
#' `maxCombine = 0` the raw aligned integral vectors are returned.
#' When `maxCombine > 0`, peaks are combined using one of two
#' methods determined by whether a reference is supplied.
#'
#' @param x
#' An object of type `aligns`.
#'
#' @param maxCombine
#' Controls peak combining in datapoints. `0` (default): off.
#' Any positive number enables peak combining. If `ref` is `NULL`, this is the
#' `range` parameter of [metabodecon::combine_peaks()]. If `ref` is supplied,
#' this is the maximum distance to the nearest reference peak. See 'Details'.
#'
#' @param ref
#' A single `align` or `decon2` object whose peaks define
#' the rows of the output matrix. If `NULL` (default), the
#' reference is auto-detected from `x` via `find_ref()`.
#'
#' @param drop_zero
#' If `TRUE`, rows where all values are zero are removed.
#'
#' @details
#' If `ref` is `NULL`, the transposed integral matrix is passed to
#' [metabodecon::combine_peaks()] with `range = maxCombine`. Partly-filled
#' neighbouring columns are merged when they share no common non-zero row.
#'
#' If `ref` is supplied, every peak is mapped to the nearest reference peak and
#' kept only if the distance is at most `maxCombine` datapoints. Areas of peaks
#' that map to the same reference peak are summed.
#'
#' Example for method 2 with ref peaks at indices 5, 9, 20 and
#' `maxCombine = 4`:
#'
#' ```txt
#' Step 1 – Build intervals [ref ± maxCombine]:
#'
#'   ref:     5              9                    20
#'            |              |                     |
#'   int:  [1 ····· 9]   [5 ···· 13]        [16 ···· 24]
#'             overlap!
#'
#' Step 2 – Shrink overlapping neighbours to midpoint.
#'          Refs 1 & 2 overlap → mid = floor((5+9)/2) = 7.
#'          Refs 2 & 3 don't  → keep maxCombine boundary.
#'
#'   ref:     5         9                         20
#'            |         |                          |
#'   int:  [1 ··· 7] [8 ·· 13]   gap        [16 ···· 24]
#'
#' Step 3 – Assign peaks to nearest ref; keep if ≤ maxCombine:
#'
#'   | Peak | Nearest | Dist | ≤ 4? | Action |
#'   |------|---------|------|------|--------|
#'   |    3 | ref 1   |    2 | yes  | keep   |
#'   |    6 | ref 1   |    1 | yes  | keep   |
#'   |    8 | ref 2   |    1 | yes  | keep   |
#'   |   11 | ref 2   |    2 | yes  | keep   |
#'   |   15 | ref 3   |    5 | no   | DROP   |
#'   |   21 | ref 3   |    1 | yes  | keep   |
#'
#' Step 4 – Sum areas (A × pi) per reference peak:
#'
#'   ref 1 ← A(3) + A(6)    (peaks at 3, 6)
#'   ref 2 ← A(8) + A(11)   (peaks at 8, 11)
#'   ref 3 ← A(21)          (peak 15 dropped)
#' ```
#'
#' @return
#' A numeric matrix with chemical shifts as rownames and
#' spectrum names as colnames.
#'
#' @author 2024-2025 Tobias Schmidt: initial version.
#'
#' @examples
#' if (interactive()) {
#'     decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
#'     aligns <- align(decons, maxCombine = 0)
#'     si_mat_0 <- get_si_mat(aligns)
#'     si_mat_1 <- get_si_mat(aligns, maxCombine = 20)
#'     si_mat_2 <- get_si_mat(aligns, maxCombine = 20, ref = aligns[[1]])
#' }
get_si_mat <- function(x, ref = NULL, drop_zero = FALSE, maxCombine = 0) {
    stopifnot(is_aligns(x))
    cs <- x[[1]]$cs
    if (maxCombine == 0) {
        mat <- sapply(x, function(xi) xi$sit$al)
        rownames(mat) <- cs
        if (drop_zero) {
            mat <- mat[rowSums(mat != 0) > 0, , drop = FALSE]
        }
    } else if (is.null(ref)) {
        mat <- sapply(x, function(xi) xi$sit$al)
        rownames(mat) <- cs
        tmat <- t(mat)
        cmat <- combine_peaks(tmat, range = maxCombine)
        mat <- t(cmat)
        if (drop_zero) {
            mat <- mat[rowSums(mat != 0) > 0, , drop = FALSE]
        }
    } else {
        get_idx <- function(vals) {
            idx <- match(vals, cs)
            if (anyNA(idx)) {
                idx <- round(convert_pos(vals, cs, seq_along(cs)))
            }
            pmin(length(cs), pmax(1L, as.integer(idx)))
        }
        ref_x0 <- if (is_align(ref)) ref$lcpar$x0_al else ref$lcpar$x0
        ref_idx <- get_idx(ref_x0)
        ns <- length(x)
        nc <- length(cs)
        tmat <- matrix(0, nrow = ns, ncol = nc)
        for (s in seq_len(ns)) {
            pi <- get_idx(x[[s]]$lcpar$x0_al)
            A <- x[[s]]$lcpar$A
            for (k in seq_along(pi)) tmat[s, pi[k]] <- tmat[s, pi[k]] + A[k] * base::pi
        }
        cmat <- combine_peaks(tmat, range = maxCombine, ref = ref_idx)
        mat <- t(cmat[, ref_idx, drop = FALSE])
        rownames(mat) <- ref_x0
        if (drop_zero) {
            mat <- mat[rowSums(mat != 0) > 0, , drop = FALSE]
        }
    }
    colnames(mat) <- get_names(x)
    mat
}

#' @noRd
#'
#' @title Combine peaks across aligned spectra
#' @description
#' Merge adjacent peak columns when no spectrum has non-zero values in both.
#' `shifted_mat[i, j]` = amplitude of spectrum i's peak at datapoint j (0 if
#' absent). With `ref = NULL`, uses greedy `merge_neighbors`; with `ref`
#' supplied, snaps each entry to the nearest reference column within `range`.
#' @param shifted_mat Numeric matrix (spectra x datapoints).
#' @param range Max merge/snap distance in datapoints. Default 5.
#' @param ref Optional integer vector of reference column indices.
#' @param lower_bound Min non-zero entries per column eligible for merging.
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
combine_peaks <- function(shifted_mat,
                          range = 5,
                          ref = NULL,
                          lower_bound = 1) {
    M <- replace(shifted_mat, is.na(shifted_mat), 0)
    if (!is.null(ref)) M <- snap_to_ref(M, ref, range)
    else               M <- merge_neighbors(M, range, lower_bound)
    M
}

#' @noRd
#' @title CluPA hierarchical-clustering alignment (Vu et al. 2011)
#' @description
#' Internal reimplementation of the CluPA algorithm used by `align()`.
#' `method = 1` delegates to `speaq::hClustAlign()`; `method = 2` (default)
#' uses metabodecon's built-in implementation.
#' @param X Numeric matrix (spectra x datapoints).
#' @param peakList List of integer peak-center index vectors, one per spectrum.
#' @param refInd Row index of the reference spectrum.
#' @param maxShift Max datapoints a peak may shift.
#' @param verbose Print progress messages.
#' @param method `1` = speaq, `2` = built-in.
#' @author
#' 2021-2024 Wolfram Gronwald: initial version.\cr
#' 2024-2025 Tobias Schmidt: refactored initial version.
dohCluster <- function(X,
                       peakList,
                       refInd = 0,
                       maxShift = 100,
                       verbose = TRUE,
                       method = 2) {
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    Y <- X
    peakListNew <- peakList
    startTime <- proc.time()
    nspec <- nrow(X)
    logf("Running dohCluster with maxShift = %d on %d spectra", maxShift, nspec)
    refSpec <- Y[refInd, ]
    for (tarInd in seq_len(nspec)) {
        if (tarInd != refInd) {
            logf("Aligning spectrum %d/%d", tarInd, nspec)
            targetSpec <- Y[tarInd, ]
            n_ref <- length(peakList[[refInd]])
            n_tar <- length(peakList[[tarInd]])
            myPeakList <- c(peakList[[refInd]], peakList[[tarInd]])
            myPeakLabel <- c(rep(1, n_ref), rep(0, n_tar))
            startP <- 1
            endP <- length(targetSpec)
            res <- if (method == 1) {
                speaq::hClustAlign(
                    refSpec, targetSpec,
                    myPeakList, myPeakLabel,
                    startP, endP,
                    maxShift = maxShift,
                    acceptLostPeak = FALSE
                )
            } else {
                hclust_align(
                    refSpec, targetSpec,
                    myPeakList, myPeakLabel,
                    startP, endP,
                    maxShift = maxShift
                )
            }
            Y[tarInd, ] <- res$tarSpec
            n_res <- length(res$peakList)
            n_end <- min(n_ref + n_tar, n_res)
            peakListNew[[tarInd]] <- res$peakList[
                (n_ref + 1):n_end
            ]
        }
    }
    elapsed <- (proc.time() - startTime)[3]
    logf("Finished dohCluster in %.1f s", elapsed)
    list("Y" = Y, "new_peakList" = peakListNew)
}

# Combine Peaks #####

#' @noRd
#'
#' @description
#' Calculates a "combine score" for a set of columns `nn` of a
#' Matrix `M`. The score describes how "beneficial" it is to merge column `n`
#' (from `nn`) into column `j`. A column `n` is considered "combinable" if there
#' is no row where columns `n` and `j` both have a non-zero element. If a column
#' is not combinable, its combine score is 0. If a column is combinable, its
#' combine score is the amount of non-zero elements in the column. This function
#' calculates the combine score for all columns `n` in `nn` and then returns the
#' calculated combine scores as vector.
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
merge_neighbors <- function(M, range, lower_bound = 1) {
    U <- M != 0
    uu <- colSums(U)
    nc <- ncol(M)
    for (i in (nrow(M) - 1):lower_bound) {
        for (j in which(uu == i)) {
            if (uu[j] == 0) next
            nn <- seq(max(1, j - range), min(nc, j + range))
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

snap_to_ref <- function(M, ref, range) {
    nr <- length(ref)
    mids <- if (nr >= 2) floor((ref[-nr] + ref[-1]) / 2) else integer()
    R <- matrix(0, nrow = nrow(M), ncol = ncol(M))
    rownames(R) <- rownames(M)
    for (i in seq_len(nrow(M))) {
        j <- which(M[i, ] != 0)
        if (length(j) == 0) next
        g <- if (nr == 1) rep.int(1L, length(j))
             else findInterval(j - 0.5, mids) + 1L
        keep <- abs(j - ref[g]) <= range
        if (!any(keep)) next
        sums <- rowsum(matrix(M[i, j[keep]], ncol = 1),
                       g[keep], reorder = FALSE)
        R[i, ref[as.integer(rownames(sums))]] <- sums[, 1]
    }
    R
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

# General Helpers #####

#' @noRd
#' @author 2024-2025 Tobias Schmidt: initial version.
get_peak_indices <- function(decon2) {
    x0 <- decon2$lcpar$x0
    cs <- decon2$cs
    convert_pos(x0, cs, seq_along(cs))
}

# Internal Helpers #####

#' @noRd
#' @description
#' Parse and validate the `...` arguments of [align()]. Issues warnings for
#' non-internal callers, errors on unknown args and invalid method. Applies
#' the verbose option as a side effect.
#' @param dots `list(...)` from the caller.
#' @param verbose Passed from `align()` to apply `local_options` side effect.
#' @return Named list with `$method` (integer) and `$full` (logical).
#' @author 2024-2026 Tobias Schmidt: initial version.
parse_align_dots <- function(dots, verbose) {
    internal <- isTRUE(dots$.internal)
    method   <- dots$method %||% 2L
    full     <- isTRUE(dots$full %||% TRUE)
    unknown  <- setdiff(names(dots), c(".internal", "method", "full", "install_deps"))
    if (length(unknown) > 0)
        stopf("Unknown args: %s", paste(unknown, collapse = ", "))
    if (!internal) {
        if ("method" %in% names(dots))
            warning("`method` is reserved for internal use. ",
                    "Pass `.internal = TRUE` to suppress.", call. = FALSE)
        if ("full" %in% names(dots))
            warning("`full` is reserved for internal use. ",
                    "Pass `.internal = TRUE` to suppress.", call. = FALSE)
        if ("install_deps" %in% names(dots))
            warning("`install_deps` is deprecated and has no effect.", call. = FALSE)
    }
    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())
    if (!is_int(method, 1) || !(method %in% 1:2))
        stop("`method` must be 1 or 2.", call. = FALSE)
    list(method = method, full = full)
}

#' @noRd
#' @description
#' Validate and coerce inputs for [align()]. Converts `ref` to decon2 if
#' needed and prepends it to `x`. Checks for uniform datapoint count and
#' rejects decons0/decons1. Checks Bioconductor deps for method 1.
#' @param x Input spectra (any decons type or aligns).
#' @param verbose Not used directly; present for symmetry with parse_align_dots.
#' @param ref Optional reference spectrum (align or decon2).
#' @param method 1 = speaq, 2 = built-in.
#' @return A `decons2` list with ref prepended when supplied.
#' @author 2024-2026 Tobias Schmidt: initial version.
chkargs_align <- function(x, verbose, ref, method) {
    if (!is.null(ref)) {
        if (is_align(ref)) class(ref) <- "decon2"
        x <- c(ref, x)
    }
    if (is_decons0(x) || is_decons1(x))
        stop("align() no longer accepts decons0 or decons1 objects. ",
             "Convert with as_decons2() first.", call. = FALSE)
    xx <- as_decons2(x)
    if (length(xx) < 2) stop("align() requires at least 2 spectra.", call. = FALSE)
    ndps <- sapply(xx, function(xi) length(xi$cs))
    if (length(unique(ndps)) != 1) {
        i <- which(ndps != ndps[1])[1]
        fmt <- paste("All spectra must have the same number of data points",
                     "but spectrum 1 has %d and spectrum %d has %d.")
        stop(sprintf(fmt, ndps[1], i, ndps[i]), call. = FALSE)
    }
    if (method == 1) {
        pkgvec <- c("MassSpecWavelet", "impute")
        inst <- sapply(pkgvec, requireNamespace, quietly = TRUE)
        if (any(!inst))
            stopf("method = 1 requires: %s. Install from Bioconductor or ",
                  "R-Universe, or use the default method = 2.",
                  paste(pkgvec[!inst], collapse = ", "))
    }
    xx
}

#' @noRd
#' @description
#' Run CluPA alignment and return continuous aligned peak-center indices.
#' Owns the `options("metabodecon.align_cache")` 1-slot cache: re-calls with
#' the same (xx, maxShift, ref) but a different maxCombine return cached pciaa.
#' @param xx A `decons2` list, with external ref prepended when applicable.
#' @param maxShift Max datapoints a peak may shift during alignment.
#' @param ref Original `ref` argument from `align()` (for cache key/refInd).
#' @param method 1 = speaq, 2 = built-in.
#' @param verbose Print progress messages.
#' @param nworkers Number of parallel workers.
#' @return List of numeric vectors (one per spectrum) of aligned peak-center indices.
#' @author 2024-2026 Tobias Schmidt: initial version.
clupa_decons <- function(xx, maxShift, ref, method, verbose, nworkers) {
    cache_key <- rlang::hash(list(xx, maxShift, ref))
    ac <- getOption("metabodecon.align_cache")
    if (!is.null(ac$key) && ac$key == cache_key) {
        logf("Skipping CluPA (cache hit)")
        return(ac$pciaa)
    }
    has_ext_ref <- !is.null(ref)
    peakList <- lapply(xx, get_peak_indices)
    find_ref_fn <- if (method == 1) speaq::findRef else find_ref
    refInd <- if (has_ext_ref) 1L else find_ref_fn(peakList)$refInd
    backend <- c("speaq", "built-in")[method]
    logf("Performing %s alignment with maxShift = %d", backend, maxShift)
    X <- do.call(rbind, lapply(xx, function(xi) xi$sit$sup))
    if (nworkers == 1) {
        obj <- dohCluster(
            X = X, peakList = peakList, refInd = refInd,
            maxShift = maxShift, verbose = verbose, method = method
        )
    } else {
        idx <- setdiff(seq_len(nrow(X)), refInd)
        k <- min(nworkers, length(idx))
        grp <- cut(seq_along(idx), breaks = k, labels = FALSE)
        chunks <- lapply(split(idx, grp), function(ids) c(refInd, ids))
        XX <- lapply(chunks, function(rows) X[rows, , drop = FALSE])
        peakLists <- lapply(chunks, function(rows) peakList[rows])
        nw_apply <- min(nworkers, length(chunks))
        objs <- mcmapply(
            nw_apply, dohCluster, XX, peakLists,
            refInd = 1, maxShift = maxShift,
            verbose = verbose, method = method
        )
        Y <- X
        new_peakList <- peakList
        for (j in seq_along(chunks)) {
            rows <- chunks[[j]]
            rows_no_ref <- rows[rows != refInd]
            if (length(rows_no_ref) == 0) next
            Y[rows_no_ref, ] <- objs[[j]]$Y[-1, , drop = FALSE]
            new_peakList[rows_no_ref] <- objs[[j]]$new_peakList[-1]
        }
        obj <- list(Y = Y, new_peakList = new_peakList)
    }
    pciaa <- obj$new_peakList
    options(metabodecon.align_cache = list(key = cache_key, pciaa = pciaa))
    pciaa
}

#' @noRd
#' @description
#' Convert continuous aligned peak-center indices (one numeric vector per
#' spectrum, as returned by `dohCluster()$new_peakList`) to a sparse signal
#' matrix. `smat[i, j] = A` if spectrum i has a peak with amplitude A at
#' datapoint j, else 0. Duplicate column indices within a row are resolved
#' by incrementing duplicates by 1 until all are distinct.
#' @param pciaa List of numeric vectors (one per spectrum).
#' @param xx `decons2` list providing `$cs` (for ndp) and `$lcpar$A`.
#' @return Numeric matrix (spectra x datapoints).
#' @author 2024-2026 Tobias Schmidt: initial version.
discretize_pciaa <- function(pciaa, xx) {
    logf("Discretizing peak center indices")
    nspec <- length(xx)
    ndp <- length(xx[[1]]$cs)
    smat <- matrix(0, nrow = nspec, ncol = ndp)
    for (i in seq_len(nspec)) {
        d <- round(pciaa[[i]])
        A <- xx[[i]]$lcpar$A
        if (length(d) == 0) stop(sprintf("No peaks found in spectrum %d", i))
        d <- pmin(ndp, pmax(1L, d))
        dups <- which(duplicated(d))
        while (length(dups) > 0) {
            fmt <- paste(
                "%d peak center indicies led to duplicate positions.",
                "Original positions: %s"
            )
            origstr <- paste(pciaa[[i]][dups], collapse = ", ")
            logf(fmt, length(dups), origstr)
            d[dups] <- d[dups] + 1
            dups <- which(duplicated(d))
        }
        smat[i, d] <- A
    }
    smat
}

#' @noRd
#' @description
#' Construct `align` objects from `decon2` objects using the combined signal
#' matrix `cmat`. Sets `lcpar$x0_al`, `sit$al`, and (when `full = TRUE`)
#' `sit$supal`. Returns an `aligns` collection.
#' @param xx `decons2` list.
#' @param cmat Signal matrix after `combine_peaks()` (equals `smat` when
#'   `maxCombine = 0`). `cmat[i, j]` = amplitude of spectrum i at datapoint j.
#' @param full If TRUE, also compute `sit$supal` (Lorentz superposition).
#' @return An `aligns` object.
#' @author 2024-2026 Tobias Schmidt: initial version.
build_aligns <- function(xx, cmat, full, has_ext_ref = FALSE) {
    cs <- xx[[1]]$cs
    n <- length(cs)
    for (i in seq_along(xx)) {
        al <- rep(0, n)
        pciac <- which(cmat[i, ] != 0)
        lcpar <- xx[[i]]$lcpar
        x0_al <- cs[pciac]
        al[pciac] <- cmat[i, pciac] * pi
        xx[[i]]$lcpar$x0_al <- x0_al
        xx[[i]]$sit$al <- al
        if (full) xx[[i]]$sit$supal <- lorentz_sup(cs, x0_al, lcpar$A, lcpar$lambda)
        class(xx[[i]]) <- "align"
    }
    aligns <- structure(xx, class = "aligns")
    if (has_ext_ref) {
        aligns <- aligns[-1]
        class(aligns) <- "aligns"
    }
    aligns
}

