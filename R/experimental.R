# Experimental #####
#
# Private, experimental functions kept for future work but intentionally
# NOT exported. They are reachable via `metabodecon:::<name>` for
# internal experiments and by the metabodecon-paper (`mdp`) package. None
# of these are part of the public v2.0 API and all are subject to change
# or removal without notice.
#
# Contents:
# - Needleman-Wunsch snap family: snap_nw(), snap_nw_blind(),
#   build_consensus() + helpers (normalize_A, snap_nw_lcpar). Alternative
#   to the public clupa() + snap_to_ref() alignment pipeline.
# - combine_peaks(): greedy post-CluPA column-merge snap + helpers
#   (combine_peaks_mat, combine_scores).
#
# All shared helpers they rely on (find_ref, dedupe_peaks, pci_on_cs,
# ensure_shared_cs) live in R/align.R.

# Needleman-Wunsch snap #####

#' @noRd
#' @title Pairwise Needleman-Wunsch snap of peak lists onto a reference
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

#' @noRd
#' @title Build a consensus peak-list reference for NW snapping
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

#' @noRd
#' @title Label-blind Needleman-Wunsch snap for fit_mdm
snap_nw_blind <- function(x, ref=NULL, maxCombine=20, w_A=0, ...) {
    stopifnot(inherits(x, "decons2"))
    cs <- ensure_shared_cs(x)
    spacing <- abs(stats::median(diff(cs)))
    gap_tol <- max(spacing, as.numeric(maxCombine) * spacing)
    pos_field <- if (!is.null(x[[1]]$lcpar$x0al)) "x0al" else "x0"
    if (is.null(ref)) {
        ref <- build_consensus(
            x, y=NULL, gap_tol=gap_tol, pos_field=pos_field
        )
    }
    out <- snap_nw(
        x, ref=ref, gap_tol=gap_tol, pos_field=pos_field, w_A=w_A
    )
    attr(out, "ref") <- ref
    out
}

# Greedy column-merge snap #####

#' @noRd
#' @title Greedy post-CluPA column-merge snap
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
