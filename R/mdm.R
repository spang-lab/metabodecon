
# API #####

#' @export
#'
#' @title Fit a Metabodecon Model
#'
#' @description
#' Performorms deconvolution, alignment and peak snapping for input spectea,
#' then fits a lasso model via [glmnet::cv.glmnet].
#'
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param y Factor vector with class labels for each spectrum.
#' @param sfr Signal free region. See [metabodecon::deconvolute()] for details.
#' @param nworkers Number of workers for parallel processing.
#' @param verbose Logical. Whether to print log messages.
#' @param use_rust Logical. Whether to use the Rust backend.
#' @param nfolds Number of folds for `cv.glmnet`. Default 10.
#' @param cadir Directory for caching deconvolution grid search results.
#' Defaults to [metabodecon::decon_cachedir()]. If `NULL`, caching is disabled.
#' Pass a custom path to use a different cache directory.
#' @param npmax Maximum number of peaks. See [metabodecon::deconvolute()] for details.
#' @param maxShift Maximum alignment shift. See [metabodecon::align()] for details.
#' @param maxSnap Maximum snap distance for feature matrix construction.
#' @param lambda Lambda selection strategy (`"lambda.1se"` or `"lambda.min"`).
#'
#' @return A list with class `mdm` and the following elements:
#' - `model`: The fitted `cv.glmnet` model object.
#' - `ref`: Reference `align` object used for snapping peaks.
#' - `args`: List of deconvolution/alignment/model parameters.
#'
#' @examples
#' \dontrun{
#'     m <- fit_mdm(spectra, y, sfr = c(11, -2))
#' }
fit_mdm <- function(
    spectra, y, sfr,
    use_rust = TRUE, npmax = 1000, maxShift = 100, maxSnap = 1,
    nworkers = 1, verbose = TRUE, nfolds = 10,
    cadir = decon_cachedir()
) {

    # When spectra carries an "mdm_hash" attribute (set by cv_mdm1),
    # reuse cached deconvolution/alignment results from a previous
    # fit_mdm call with the same spectra but different grid params.
    h <- attr(spectra, "mdm_hash")
    ca <- getOption("metabodecon.fit_mdm.cache")
    hit <- !is.null(h) && !is.null(ca) && identical(ca$hash, h)

    if (hit && identical(ca$npmax, npmax)) {
        logv("Reusing cached deconvolution (npmax=%d)", npmax)
        decons <- ca$decons
    } else {
        logf("Deconvoluting spectra with npmax=%d", npmax)
        decons <- deconvolute(
            x=spectra, sfr=sfr, verbose=verbose, use_rust=use_rust,
            npmax=npmax, nworkers=nworkers, cadir=cadir
        )
    }

    if (hit && identical(ca$npmax, npmax) && identical(ca$maxShift, maxShift)) {
        logv("Reusing cached alignment (maxShift=%d)", maxShift)
        aligns <- ca$aligns
    } else {
        logf("Aligning spectra with maxShift=%d", maxShift)
        aligns <- align(
            x=decons, maxShift=maxShift, maxCombine=0,
            verbose=verbose, nworkers=nworkers
        )
    }

    if (!is.null(h)) {
        options(metabodecon.fit_mdm.cache = list(
            hash=h, npmax=npmax, maxShift=maxShift,
            decons=decons, aligns=aligns
        ))
    }
    X <- {
        logv("Constructing feature matrix with maxSnap=%.1f", maxSnap)
        ref <- mdm_get_ref(aligns)
        t(get_si_mat(aligns, maxSnap=maxSnap, ref=ref))
    }
    model <- {
        logf("Fitting cv.glmnet")
        glmnet::cv.glmnet(x=X, y=y, family="binomial", alpha=1, nfolds=nfolds)
    }
    mdm <- {
        logf("Constructing and returning mdm object")
        meta <- list(
            model="lasso", snames=get_names(spectra), y=y,
            sfr=sfr, use_rust=use_rust, npmax=npmax,
            maxShift=maxShift, maxSnap=maxSnap,
        )
        structure(list(model=model, ref=ref, meta=meta), class="mdm")
    }
    mdm
}

#' @export
#' @title Deconvolute, Align and Classify Spectra in Cross-Validation
#' @description
#' Performs deconvolution, alignment and peak-snapping over a grid of parameters
#' and calls [glmnet::cv.glmnet()] on each corresponding matrix. Calculates the
#' corresponding `cv.glmnet` performances. The deconvolution, alignment and
#' snapping parameters that enable the best performances are considered optimal
#' and are used for fitting the final model with the corresponding optimal
#' lambda value.
#'
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param y Factor vector with class labels for each spectrum.
#' @param sfr Signal free region. See [metabodecon::deconvolute()] for details.
#' @param nworkers Number of workers for parallel grid deconvolution.
#' @param verbose Logical. Whether to print log messages.
#' @param use_rust Logical. Whether to use the Rust backend.
#' @param nfolds Number of folds for inner cv.glmnet. Default 10.
#' @param cadir Directory for caching deconvolution grid search results.
#' Defaults to [metabodecon::decon_cachedir()]. If `NULL`, caching is disabled.
#' Pass a custom path to use a different cache directory.
#'
#' @return
#' A object of class `mdm` as returned by [metabodecon::fit_mdm()], with two
#' additional elements `pgrid` and `ibest`, containing the performance grid and
#' index of the best parameter combination.
#'
#' @examples
#' \dontrun{
#'      aki <- read_aki_data()
#'      spectra <- aki$spectra
#'      y <- factor(aki$meta$type, levels=c("Control", "AKI"))
#'      nworkers <- ceiling(parallel::detectCores() / 2)
#'      stub(cv_mdm, spectra=spectra, y=y, sfr=c(11, -2), nworkers=nworkers)
#'      mdm <- cv_mdm(spectra, y, sfr=c(11, -2))
#' }
cv_mdm <- function(
    spectra, y, sfr, use_rust = TRUE,
    npmax = NULL, maxShift = seq(50, 250, 50), maxSnap = seq(0.5, 2.0, 0.5),
    nworkers = 1, verbose = TRUE, nfolds = 10,
    cadir = decon_cachedir()
) {

    logv("Starting cache initialization")
    gridlist <- grid_deconvolute_spectra(
        x=spectra, sfr=sfr, verbose=verbose, nw=nworkers,
        use_rust=use_rust, cadir=cadir
    )
    npmax <- mdm_get_npmaxs(gridlist)
    nw <- min(nworkers, length(npmax))

    attr(spectra, "mdm_hash") <- rlang::hash(spectra)

    logv("Starting grid search over npmax, maxShift, maxSnap and lambda")
    npmaxs_list <- split(npmax, cut2(npmax, nw))
    pgrids <- mcmapply(
        nw = nw, FUN = cv_mdm1, npmax = npmaxs_list,
        MoreArgs = list(
            spectra=spectra, y=y, sfr=sfr, use_rust=use_rust,
            maxShift=maxShift, maxSnap=maxSnap,
            nfolds=nfolds, verbose=verbose, cadir=cadir
        )
    )
    pgrid <- rbindlist(pgrids)
    ibest <- which.min(pgrid$cvm)
    best <- as.list(pgrid[ibest, ])
    fmt <- "Best: npmax=%d, maxShift=%d, maxSnap=%.1f, cvm=%.4f"
    logv(fmt, best$npmax, best$maxShift, best$maxSnap, best$cvm)

    logv("Fitting final model with best parameters")
    mdm <- fit_mdm(
        spectra=spectra, y=y, sfr=sfr, use_rust=use_rust,
        npmax=best$npmax, maxShift=best$maxShift, maxSnap=best$maxSnap,
        nworkers=nw, verbose=verbose, nfolds=nfolds, cadir=cadir
    )
    mdm$pgrid <- pgrid
    mdm
}

#' @noRd
#' @description
#' Strict single-core helper for [metabodecon::cv_mdm()].
#' Assumes pre-initialized deconvolution cache.
cv_mdm1 <- function(
    spectra, y, sfr, use_rust = TRUE,
    npmax = 1000, maxShift = seq(50, 250, 50), maxSnap = seq(0.5, 2.0, 0.5),
    nfolds = 10, verbose = TRUE, cadir = decon_cachedir()
) {

    logv("Building parameter grid")
    P <- expand.grid(
        npmax = npmax, maxShift = maxShift, maxSnap = maxSnap,
        lambda = c("lambda.1se", "lambda.min"), cvm = NA_real_,
        cvsd = NA_real_, KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    P <- P[order(P$npmax, P$maxShift, P$maxSnap), ]
    rownames(P) <- NULL
    np <- nrow(P)

    on.exit(options(metabodecon.fit_mdm.cache = NULL), add = TRUE)

    logv("Starting grid search (%d combinations)", np)
    for (i in seq_len(np)) {
        pc <- as.list(P[i, ])
        fmt <- "[%d/%d]: npmax=%d, maxShift=%d, maxSnap=%.1f"
        logv(fmt, i, np, pc$npmax, pc$maxShift, pc$maxSnap)
        mdm <- fit_mdm(
            spectra=spectra, y=y, sfr=sfr, use_rust=use_rust,
            npmax=pc$npmax, maxShift=pc$maxShift,
            maxSnap=pc$maxSnap, nworkers=1, verbose=FALSE,
            nfolds=nfolds, cadir=cadir
        )
        m <- mdm$model
        idx <- which(m$lambda == m[[pc$lambda]])
        P$cvm[i] <- m$cvm[idx]
        P$cvsd[i] <- m$cvsd[idx]
    }
    P
}

#' @export
#' @title Estimate MDM Performance via Nested Cross-Validation
#' @description
#' Splits the data into `nfo` outer folds. For each fold, calls
#' [metabodecon::cv_mdm()] on the training subset to select the best preprocessing
#' parameters and fit a model, then evaluates on the held-out test fold.
#' Grid search results are shared across all folds via a temporary cache
#' directory.
#' @details
#' Each outer worker keeps its own copy of `spectra` and `y`, so memory use
#' grows roughly with `nwo`, not with `nfo`. Large `nwo` values can therefore
#' still require substantial RAM for large input datasets.
#'
#' @inheritParams cv_mdm
#' @param nwo Number of workers for the outer cross-validation. Each worker gets
#' a full copy of `spectra` and `y`.
#' @param nwi Number of workers used inside each `cv_mdm()` fit and prediction.
#' @param nfo Number of outer folds. Default 3.
#' @param nfi Number of folds for inner cv.glmnet. Default 3.
#' @param seed Random seed for fold assignment.
#'
#' @return A list with elements:
#' - `models`: List of objects with classes `cv` and `mdm`, one per fold.
#' - `predictions`: Data frame with columns `fold`, `true`, `prob`, `pred`.
#'
#' @examples
#' \dontrun{
#'      aki <- read_aki_data()
#'      spectra <- aki$spectra
#'      y <- factor(aki$meta$type, levels = c("Control", "AKI"))
#'      stub(benchmark_mdm, spectra=spectra, y=y, sfr = c(11, -2), nwo = 5, nwi = 11)
#'      bm <- benchmark_mdm(spectra=spectra, y=y, sfr = c(11, -2), nwo = 5, nwi = 11)
#' }
benchmark_mdm <- function(
    spectra, y, sfr, nwo = 1, nwi = 1, verbose = TRUE, use_rust = TRUE,
    nfo = 5, nfi = 10, seed = 1, cadir = decon_cachedir()
) {

    # Pre-warm cache once for all spectra
    ns <- length(spectra)
    nw1 <- min(nwo * nwi, ns)
    logv("Initializing cache for %d spectra using %d workers", ns, nw1)
    grid_deconvolute_spectra(
        x = spectra, sfr = sfr, verbose = verbose,
        nw = nw1, use_rust = use_rust, cadir = cadir
    )
    logv("Finished cache initialization")

    # Assign folds
    set.seed(seed)
    ns <- length(spectra) # number of spectra
    fids <- seq_len(nfo) # fold indices (1:3 for nfo = 3)
    sids <- seq_len(ns) # spectra indices (1:9 for ns = 9)
    sfmap <- sample(rep_len(fids, ns)) # fid per spectrum (3 2 3 1 3 2 1 1 2)
    te_list <- split(sids, sfmap) # [(4 7 8), (2 6 9), (1 3 5)]

    # Fit cv_mdm on each fold's training set
    logv("Fitting cv_mdm for %d outer folds", nfo)
    nw2 <- min(nwo, nfo)
    fids_list <- if (nw2 > 1) cut(fids, nw2, labels = FALSE) else list(fids)
    ctx_list <- list(list(
        spectra=spectra, y=y, te_list=te_list, sfr=sfr, nworkers=nwi,
        verbose=verbose, use_rust=use_rust, nfolds=nfi, cadir=cadir
    ))
    models <- mcmapply(nw2, cv_mdm_on_folds, fids_list, ctx_list)
    models <- unlist(models, recursive=FALSE, use.names=FALSE)

    # Predict on each fold's test set
    logv("Predicting on held-out test folds")
    ctx_list[[1]]$models <- models
    preds <- mcmapply(nw2, predict_mdm_on_folds, fids_list, ctx_list)
    preds <- unlist(preds, recursive = FALSE, use.names = FALSE)

    Y <- do.call(rbind, preds)
    acc <- mean(Y$true == Y$pred)
    logv("Nested CV accuracy: %.2f%%", acc * 100)
    list(models = models, predictions = Y)
}

benchmark_aki <- function() {
    aki <- read_aki_data()
    spectra <- aki$spectra
    y <- factor(aki$meta$type, levels = c("Control", "AKI"))
    sfr <- c(11, -2)
    logv("Starting cache initialization")
    gridlist <- grid_deconvolute_spectra(
        x=spectra, sfr=sfr, verbose=TRUE, nw=6, use_rust=TRUE
    )
    bm <- benchmark_mdm(spectra=spectra, y=y, sfr=sfr, nwo = 1)
    bm
}

# Benchhmark helpers #####

#' @noRd
#' @title Fit cv_mdm models on outer folds
#'
#' @description
#' Internal worker helper for [metabodecon::benchmark_mdm()]. Given a vector of
#' outer-fold indices and a shared context object, it fits
#' [metabodecon::cv_mdm()] on the corresponding training subsets and
#' returns the per-fold results as a list.
#'
#' @details
#' This helper is designed only for use with [metabodecon::mcmapply()] inside
#' [metabodecon::benchmark_mdm()]. Each worker receives one chunk of fold ids
#' together with one shared context object and processes the chunk locally via
#' [base::lapply()]. This keeps the context explicit while still serializing
#' `spectra` only once per worker.
#'
#' @param fids Integer vector of outer-fold indices.
#' @param ctx Worker context.
#'
cv_mdm_on_folds <- function(fids, ctx) {
    lapply(fids, function(i) {
        te <- ctx$te_list[[i]]
        keep <- rep(TRUE, length(ctx$y))
        keep[te] <- FALSE
        cv_mdm(
            spectra = ctx$spectra[keep], y = ctx$y[keep],
            sfr = ctx$sfr, nworkers = ctx$nworkers,
            verbose = ctx$verbose, use_rust = ctx$use_rust,
            nfolds = ctx$nfolds, cadir = ctx$cadir
        )
    })
}

predict_mdm_on_folds <- function(fids, ctx) {
    lapply(fids, function(i) {
        te <- ctx$te_list[[i]]
        prob <- predict(
            object = ctx$models[[i]], newdata = ctx$spectra[te],
            type = "prob", nworkers = ctx$nworkers
        )
        pred <- ifelse(prob > 0.5, levels(ctx$y)[2], levels(ctx$y)[1])
        pred <- factor(pred, levels = levels(ctx$y))
        data.frame(fold = i, true = ctx$y[te], prob = prob, pred = pred)
    })
}

# Normalisation helpers #####

#' @title Build quantile reference
#' @description Computes the mean sorted profile across rows.
#' @param X Numeric feature matrix.
#' @return Numeric vector of reference quantiles.
#' @examples
#' X <- matrix(c(1, 3, 2, 4), nrow = 2)
#' mdm_get_quantile_reference(X)
#' @noRd
mdm_get_quantile_reference <- function(X) {
    mean_sorted <- apply(X, 1, sort)
    rowMeans(mean_sorted)
}

#' @title Apply quantile normalization
#' @description Maps each row to a shared reference distribution.
#' @param X Numeric feature matrix.
#' @param ref Quantile reference vector.
#' @return Quantile-normalized matrix.
#' @examples
#' X <- matrix(rnorm(20), nrow = 4)
#' r <- mdm_get_quantile_reference(X)
#' mdm_quantile_normalize(X, r)
#' @noRd
mdm_quantile_normalize <- function(X, ref) {
    Xn <- X
    for (i in seq_len(nrow(X))) {
        o <- order(X[i, ])
        Xn[i, o] <- ref
    }
    Xn
}

#' @title Build PQN reference
#' @description Computes feature-wise medians for probabilistic quotient normalization.
#' @param X Numeric feature matrix.
#' @return Numeric reference spectrum.
#' @examples
#' X <- matrix(rnorm(24), nrow = 6)
#' mdm_get_pqn_reference(X)
#' @noRd
mdm_get_pqn_reference <- function(X) {
    apply(X, 2, median)
}

#' @title Apply PQN normalization
#' @description Divides each row by its median quotient to the reference profile.
#' @param X Numeric feature matrix.
#' @param ref Reference spectrum.
#' @return PQN-normalized matrix.
#' @examples
#' X <- matrix(abs(rnorm(24)), nrow = 6)
#' r <- mdm_get_pqn_reference(X)
#' mdm_pqn_normalize(X, r)
#' @noRd
mdm_pqn_normalize <- function(X, ref) {
    Q <- sweep(X, 2, ref, "/")
    s <- apply(Q, 1, median)
    s[!is.finite(s) | s == 0] <- 1
    sweep(X, 1, s, "/")
}

#' @title Resolve creatinine bin
#' @description Finds creatinine bin index from input or bin-center metadata.
#' @param X Numeric feature matrix.
#' @param cre_idx Optional integer index to force a specific creatinine bin.
#' @return Integer scalar bin index.
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' attr(X, "bin_centers") <- seq(9.5, 0.5, length.out = ncol(X))
#' mdm_get_creatinine_ref(X)
#' @noRd
mdm_get_creatinine_ref <- function(X, cre_idx = NULL) {
    if (!is.null(cre_idx)) {
        return(as.integer(cre_idx))
    }
    cen <- attr(X, "bin_centers")
    if (is.null(cen)) {
        return(1L)
    }
    as.integer(which.min(abs(cen - 3.05)))
}

#' @title Fit normalizer metadata
#' @description Creates a method/ref pair for later normalization of new data.
#' @param X Numeric feature matrix.
#' @param method Normalization method.
#' @param cre_idx Optional creatinine bin index.
#' @return List with elements `method` and `ref`.
#' @examples
#' X <- matrix(abs(rnorm(30)), nrow = 6)
#' mdm_fit_normalizer(X, method = "median")
#' @noRd
mdm_fit_normalizer <- function(X,
                               method = c(
                                   "quantile", "pqn", "creatinine",
                                   "sum", "median", "none"
                               ),
                               cre_idx = NULL) {
    method <- match.arg(method)
    ref <- NULL
    if (method == "quantile") ref <- mdm_get_quantile_reference(X)
    if (method == "pqn") ref <- mdm_get_pqn_reference(X)
    if (method == "creatinine") ref <- mdm_get_creatinine_ref(X, cre_idx)
    list(method = method, ref = ref)
}

#' @title Apply normalization by method
#' @description Applies one of the supported normalization methods to matrix rows.
#' @param X Numeric feature matrix.
#' @param normalizer List with `method` and `ref`.
#' @return Normalized feature matrix.
#' @examples
#' X <- matrix(abs(rnorm(30)), nrow = 6)
#' nrm <- mdm_fit_normalizer(X, method = "sum")
#' mdm_apply_normalizer(X, nrm)
#' @noRd
mdm_apply_normalizer <- function(X, normalizer) {
    m <- normalizer$method
    r <- normalizer$ref
    if (m == "none") {
        return(X)
    }
    if (m == "quantile") {
        return(mdm_quantile_normalize(X, r))
    }
    if (m == "pqn") {
        return(mdm_pqn_normalize(X, r))
    }
    if (m == "creatinine") {
        s <- X[, r]
        s[!is.finite(s) | s == 0] <- median(s[s > 0], na.rm = TRUE)
        s[!is.finite(s) | s == 0] <- 1
        return(sweep(X, 1, s, "/"))
    }
    if (m == "sum") {
        s <- rowSums(X)
        s[!is.finite(s) | s == 0] <- 1
        return(sweep(X, 1, s, "/"))
    }
    if (m == "median") {
        s <- apply(X, 1, median)
        s[!is.finite(s) | s == 0] <- 1
        return(sweep(X, 1, s, "/"))
    }
    stop("Unsupported normalisation method: ", m)
}

# Classification helpers #####

#' @title Compute rank-based AUC
#' @description Computes area under the ROC curve using rank statistics.
#' @param y Binary labels coded as 0/1 or coercible to integer.
#' @param score Numeric prediction scores.
#' @return Numeric scalar AUC or `NA_real_` if one class is missing.
#' @examples
#' y <- c(0, 0, 1, 1)
#' s <- c(0.1, 0.3, 0.6, 0.8)
#' mdm_auc(y, s)
#' @noRd
mdm_auc <- function(y, score) {
    y <- mdm_as_binary01(y)
    pos <- y == 1
    n1 <- sum(pos)
    n0 <- sum(!pos)
    if (n1 == 0 || n0 == 0) {
        return(NA_real_)
    }
    r <- rank(score)
    (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @title Select top t-score features
#' @description Ranks columns by absolute two-class t-like score.
#' @param X Numeric matrix with samples in rows and features in columns.
#' @param y Binary labels coded as 0/1.
#' @param nfeat Number of top features to return.
#' @return Integer vector of selected feature indices.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40), nrow = 8)
#' y <- rep(0:1, each = 4)
#' mdm_select_top_features(X, y, nfeat = 3)
#' @noRd
mdm_select_top_features <- function(X, y, nfeat) {
    y <- mdm_as_binary01(y)
    i1 <- which(y == 1)
    i0 <- which(y == 0)
    m1 <- colMeans(X[i1, , drop = FALSE])
    m0 <- colMeans(X[i0, , drop = FALSE])
    v1 <- apply(X[i1, , drop = FALSE], 2, var)
    v0 <- apply(X[i0, , drop = FALSE], 2, var)
    s <- sqrt(v1 / length(i1) + v0 / length(i0))
    t <- abs((m1 - m0) / s)
    t[!is.finite(t)] <- 0
    order(t, decreasing = TRUE)[seq_len(min(nfeat, ncol(X)))]
}

# Utility helpers #####

mdm_as_binary01 <- function(y) {
    as_binary01(y)
}

#' @title Format mean and standard deviation
#' @description Creates compact strings of the form mean ± sd.
#' @param x Mean value.
#' @param y Standard-deviation value.
#' @param d Digits after decimal point.
#' @return Character scalar.
#' @examples
#' mdm_fmt_mean_sd(0.8123, 0.0345, d = 3)
#' @noRd
mdm_fmt_mean_sd <- function(x, y, d = 3) {
    sprintf(paste0("%.", d, "f \u00b1 %.", d, "f"), x, y)
}

mdm_get_ref <- function(x) {
    pl <- lapply(x, get_peak_indices)
    idx <- find_ref(pl)$refInd
    x[[idx]]
}

mdm_get_npmaxs <- function(gridlist) {
    np <- unname(unlist(lapply(gridlist, function(g) g$np)))
    np <- np[np > 0]
    if (length(np) == 0) {
        stop("All grid deconvolutions produced zero peaks.")
    }
    np_hi <- ceiling(max(np) / 100) * 100
    np_lo <- ceiling(min(np) / 100) * 100
    cands <- seq(np_lo, np_hi, by = 100)
    support <- vapply(cands, function(npm) {
        mean(vapply(gridlist, function(g) {
            any(g$np > 0 & g$np < npm)
        }, logical(1)))
    }, numeric(1))
    idx <- which(support >= 0.5)
    lower <- if (length(idx) == 0) np_hi else cands[idx[1]]
    seq(lower, np_hi, by = 100)
}

# S3 methods #####

#' @title Predict from mdm model
#' @description Predicts probabilities, classes, or link scores from an `mdm`
#' object. When `newdata` is a spectra object, the spectra are deconvoluted,
#' aligned and snapped to the reference stored in the model before prediction.
#' When `newdata` is a numeric matrix, it is used directly as the feature
#' matrix.
#' @param object Fitted `mdm` object.
#' @param newdata Spectra object or numeric feature matrix.
#' @param type Prediction type, one of `"prob"`, `"class"`, `"link"`.
#' @param s Regularization value for lasso predictions.
#' @param nworkers Number of workers used when `newdata` is a spectra object.
#' @param ... Additional arguments, currently unused.
#' @return Numeric vector of probabilities, classes, or link scores.
#' @examples
#' \dontrun{
#'   m <- cv_mdm(spectra, y, sfr = c(11, -2))
#'   predict(m, test_spectra, type = "prob")
#' }
predict.mdm <- function(object,
                        newdata,
                        type = c("prob", "class", "link"),
                        s = "lambda.min",
                        nworkers = 1,
                        ...) {
    type <- match.arg(type)

    if (is_spectra(newdata)) {
        meta <- object$meta
        decons <- deconvolute(
            x = newdata, sfr = meta$sfr,
            verbose = FALSE, use_rust = meta$use_rust,
            npmax = meta$npmax, nworkers = nworkers
        )
        als <- align(
            decons,
            maxShift = meta$maxShift,
            maxCombine = 0,
            verbose = FALSE,
            nworkers = nworkers,
            ref = object$ref
        )
        Xn <- t(get_si_mat(als,
            maxSnap = meta$maxSnap,
            ref = object$ref))
    } else {
        norm <- list(method = object$normalisation,
                     ref = object$ref)
        Xn <- mdm_apply_normalizer(newdata, norm)
    }

    if (object$meta$model == "svm") {
        idx <- object$meta$idx
        pred <- predict(
            object$model,
            Xn[, idx, drop = FALSE],
            decision.values = TRUE
        )
        score <- as.numeric(attr(pred, "decision.values"))
        if (isTRUE(object$meta$flip)) score <- -score
        score[!is.finite(score)] <- 0
        prob <- 1 / (1 + exp(-score))
    } else {
        score <- as.numeric(stats::predict(object$model,
            newx = Xn,
            s = s, type = "link"
        ))
        prob <- as.numeric(stats::predict(object$model,
            newx = Xn,
            s = s, type = "response"
        ))
    }

    if (type == "prob") {
        return(prob)
    }
    if (type == "link") {
        return(score)
    }
    as.integer(prob >= 0.5)
}

#' @title Print mdm object
#' @description Prints a compact model summary for `mdm` objects.
#' @param x Fitted `mdm` object.
#' @param ... Additional arguments, currently unused.
#' @return Invisibly returns `x`.
#' @examples
#' m <- structure(
#'   list(model = NULL, normalisation = "none", ref = NULL,
#'        meta = list(model = "lasso")),
#'   class = "mdm"
#' )
#' print(m)
print.mdm <- function(x, ...) {
    cat("metabodecon model (mdm)\n")
    cat("  model:         ", x$meta$model, "\n", sep = "")
    cat("  normalisation: ", x$normalisation, "\n", sep = "")
    if (x$meta$model == "svm") {
        cat("  nfeat:         ", x$meta$nfeat, "\n", sep = "")
        cat("  cost:          ", x$meta$cost, "\n", sep = "")
        cat("  gamma:         ", x$meta$gamma, "\n", sep = "")
    }
    invisible(x)
}

#' @title Extract mdm coefficients
#' @description Returns lasso coefficients for `mdm` objects fitted with lasso.
#' @param object Fitted `mdm` object.
#' @param ... Additional arguments passed to `stats::coef`.
#' @return Coefficient object for lasso, otherwise `NULL`.
#' @examples
#' \dontrun{
#'   m <- cv_mdm(spectra, y, sfr = c(11, -2))
#'   coef(m)
#' }
coef.mdm <- function(object, ...) {
    if (object$meta$model == "lasso") {
        return(stats::coef(object$model, s = "lambda.min", ...))
    }
    NULL
}

#' @title Plot mdm object
#' @description Plots lasso path for lasso models or support-vector count for SVM.
#' @param x Fitted `mdm` object.
#' @param ... Additional plotting arguments.
#' @return Invisibly returns `NULL`.
#' @examples
#' \dontrun{
#'   m <- cv_mdm(spectra, y, sfr = c(11, -2))
#'   plot(m)
#' }
plot.mdm <- function(x, ...) {
    if (x$meta$model == "lasso") {
        graphics::plot(x$model, ...)
    } else {
        nsv <- length(x$model$index)
        graphics::barplot(nsv,
            names.arg = "support vectors",
            main = "SVM model", ...
        )
    }
    invisible(NULL)
}

#' @title Summarize mdm object
#' @description Builds a compact summary list for mdm objects.
#' @param object Fitted `mdm` object.
#' @param ... Additional arguments, currently unused.
#' @return Object of class `summary.mdm`.
#' @examples
#' m <- structure(
#'   list(model = NULL, normalisation = "none", ref = NULL,
#'        meta = list(model = "svm")),
#'   class = "mdm"
#' )
#' summary(m)
summary.mdm <- function(object, ...) {
    out <- list(
        model = object$meta$model,
        normalisation = object$normalisation,
        ref_len = if (is.null(object$ref)) 0 else length(object$ref)
    )
    class(out) <- "summary.mdm"
    out
}

#' @title Print summary.mdm object
#' @description Prints formatted output for objects of class `summary.mdm`.
#' @param x Object of class `summary.mdm`.
#' @param ... Additional arguments, currently unused.
#' @return Invisibly returns `x`.
#' @examples
#' s <- structure(list(model = "svm", normalisation = "none", ref_len = 0),
#'   class = "summary.mdm"
#' )
#' print(s)
print.summary.mdm <- function(x, ...) {
    cat("Summary of mdm\n")
    cat("  model:         ", x$model, "\n", sep = "")
    cat("  normalisation: ", x$normalisation, "\n", sep = "")
    cat("  ref length:    ", x$ref_len, "\n", sep = "")
    invisible(x)
}
