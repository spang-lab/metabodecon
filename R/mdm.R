#' @title Metabodecon model module
#' @description Helper and API functions for model fitting, prediction, and
#' benchmark evaluation in the AKI reproduction workflow.
#'
#' + API functions:
#'   - fit_mdm: Fit an mdm (metabodecon model) with given parameters
#'   - cv_mdm: Fit an mdm with internal cv-based parameter selection
#'   - benchmark_cv_mdm: Use outer cv to estimate performance of cv_mdm
#' + Normalisation helpers:
#'   - mdm_get_quantile_reference: Build quantile-normalization reference from sample-wise sorted values.
#'   - mdm_quantile_normalize: Apply quantile normalization to each sample row with a shared reference.
#'   - mdm_get_pqn_reference: Build PQN reference as median profile across samples.
#'   - mdm_pqn_normalize: Apply probabilistic quotient normalization using a reference profile.
#'   - mdm_get_creatinine_ref: Resolve creatinine reference bin index from metadata or fallback.
#'   - mdm_fit_normalizer: Fit normalization metadata used during training and prediction.
#'   - mdm_apply_normalizer: Apply selected normalization to a feature matrix.
#' + Classification helpers:
#'   - mdm_auc: Compute AUC from binary labels and numeric scores using rank sums.
#'   - mdm_select_top_features: Rank features by Welch-like t-score and return top indices.
#' + Utility helpers:
#'   - mdm_as_binary01: Convert factor/integer labels to 0/1.
#'   - mdm_fmt_mean_sd: Format mean ± sd strings for result tables.
#'   - mdm_data_signature: Compute stable hash for matrix/label pairs.
#' + S3 methods:
#'   - predict.mdm: Predict probabilities, classes, or linear scores for mdm objects.
#'   - print.mdm: Print compact mdm summary.
#'   - coef.mdm: Extract lasso coefficients from mdm objects.
#'   - plot.mdm: Plot lasso paths or SVM support-vector count summaries.
#'   - summary.mdm: Build structured mdm summary object.
#'   - print.summary.mdm: Print formatted mdm summary.
#' @noRd
NULL

# API #####

#' @export
#'
#' @title Fit a Metabodecon Model
#'
#' @description Deconvolutes and aligns spectra, then fits a lasso model via
#' `cv.glmnet`.
#'
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param y Factor vector with class labels for each spectrum.
#' @param sfr Signal free region. See [metabodecon::deconvolute()] for details.
#' @param nworkers Number of workers for parallel processing.
#' @param verbose Logical. Whether to print log messages.
#' @param use_rust Logical. Whether to use the Rust backend.
#' @param nfolds Number of folds for `cv.glmnet`. Default 10.
#' @param cadir Directory for caching grid search results.
#' @param npmax Maximum number of peaks. See [metabodecon::deconvolute()] for details.
#' @param maxShift Maximum alignment shift. See [metabodecon::align()] for details.
#' @param maxSnap Maximum snap distance for feature matrix construction.
#' @param lambda Lambda selection strategy (`"lambda.1se"` or `"lambda.min"`).
#'
#' @return A list with class `mdm` and the following elements:
#' - `model`: The fitted `cv.glmnet` model object.
#' - `normalisation`: Normalization method used (currently always "none").
#' - `ref`: Reference metadata for normalization (currently always `NULL`).
#' - `meta`: List of metadata about the fitting process and parameters.
#'
#' @examples
#' \dontrun{
#'     m <- fit_mdm(spectra, y, sfr = c(11, -2))
#' }
fit_mdm <- function(
    spectra,
    y,
    sfr,
    nworkers = 4,
    verbose = TRUE,
    use_rust = TRUE,
    nfolds = 10,
    cadir = cachedir("grid_deconvolute_spectrum", FALSE),
    npmax = 1000,
    maxShift = 100,
    maxSnap = 1,
    lambda = "lambda.1se"
) {
    logf("Deconvoluting spectra (npmax=%d)", npmax)
    decons <- deconvolute(
        x = spectra, sfr = sfr, verbose = verbose,
        use_rust = use_rust, npmax = npmax,
        nworkers = nworkers, cachedir = cadir
    )
    logf("Aligning spectra (maxShift=%d)", maxShift)
    als <- align(
        decons,
        maxShift = maxShift,
        maxCombine = 0,
        verbose = verbose,
        nworkers = nworkers
    )
    X <- t(get_si_mat(als, maxSnap = maxSnap))
    nf <- min(nfolds, nrow(X))
    logf("Fitting cv.glmnet (nfolds=%d, lambda=%s)", nf, lambda)
    model <- glmnet::cv.glmnet(
        x = X, y = y,
        family = "binomial", alpha = 1, nfolds = nf
    )
    logf("Finished fit_mdm")
    structure(list(
        model = model,
        normalisation = "none",
        ref = NULL,
        meta = list(
            model = "lasso",
            lambda = lambda,
            npmax = npmax,
            maxShift = maxShift,
            maxSnap = maxSnap
        )
    ), class = "mdm")
}

#' @export
#' @title Deconvolute, Align and Classify Spectra in Cross-Validation
#' @description
#' Performs deconvolution+alignment over a grid of parameters and calculates the
#' corresponding cv.glmnet-lasso-lambda.min performance. The
#' deconvolution+alignment parameters that enable the best cv glmnet
#' performances are considered optimal and are used for fitting the final model
#' with the corresponding optimal lambda value.
#'
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param y Factor vector with class labels for each spectrum.
#' @param sfr Signal free region. See [metabodecon::deconvolute()] for details.
#' @param nworkers Number of workers for parallel grid deconvolution.
#' @param verbose Logical. Whether to print log messages.
#' @param use_rust Logical. Whether to use the Rust backend.
#' @param nfolds Number of folds for inner cv.glmnet. Default 10.
#' @param cadir Directory for caching grid search results. If `NULL`
#' (default), a temporary session-scoped directory is used. Pass a custom
#' path to share the cache across parallel calls.
#'
#' @return
#' A list (class `mdm`) with the following elements:
#' - `model`: The fitted `cv.glmnet` model object.
#' - `pgrid`: Data frame of parameter combinations and their cv.glmnet
#'   performances (columns: `npmax`, `maxShift`, `maxSnap`,
#'   `cvm`).
#' - `best`: Named list of the best parameter combination.
#'
#' @examples
#' \dontrun{
#'      aki <- read_aki_data()
#'      spectra <- aki$spectra
#'      y <- factor(aki$meta$type, levels = c("Control", "AKI"))
#'      mdm <- cv_mdm(spectra, y, sfr = c(11, -2))
#' }
cv_mdm <- function(
    spectra,
    y,
    sfr,
    nworkers = 4,
    verbose = TRUE,
    use_rust = TRUE,
    nfolds = 10,
    cadir = cachedir("grid_deconvolute_spectrum", FALSE)
) {

    # Pre-warm the grid cache so that all subsequent deconvolute() calls
    # with npmax >= 1 can reuse the cached per-spectrum grid results.
    logf("Starting cache initialization")
    gridlist <- grid_deconvolute_spectra(
        x=spectra, sfr=sfr, verbose=verbose, nw=nworkers,
        use_rust=use_rust, cachedir=cadir
    )
    logf("Finished cache initialization")

    # Build parameter grid from observed peak counts
    logf("Building parameter grid")
    np <- unname(unlist(lapply(gridlist, function(g) g$np)))
    np_hi <- ceiling(max(np) / 100) * 100
    np_lo <- ceiling(min(np[np != 0]) / 100) * 100
    npmaxs <- seq(np_lo, np_hi, by = 100)
    maxShifts <- seq(25, 250, by = 25)
    maxSnaps <- c(0.5, 1, 1.5, 2, 2.5, 3)
    P <- expand.grid(
        npmax = npmaxs, maxShift = maxShifts, maxSnap = maxSnaps
    )
    P$cvm <- NA_real_

    logf("Starting deconvolution + alignment + classification grid")
    row <- 1
    for (i in seq_along(npmaxs)) {
        logf("npmax=%d", npmaxs[i])
        decons <- deconvolute(
            x=spectra, sfr=sfr, verbose=verbose,
            use_rust=use_rust, npmax=npmaxs[i],
            nworkers=nworkers, cachedir=cadir
        )
        for (j in seq_along(maxShifts)) {
            logf("  maxShift=%d", maxShifts[j])
            als <- align(
                decons, maxShift=maxShifts[j], maxCombine=0,
                verbose=verbose, nworkers=nworkers
            )
            for (k in seq_along(maxSnaps)) {
                X <- t(get_si_mat(als, maxSnap=maxSnaps[k]))
                cvm <- glmnet::cv.glmnet(
                    x=X, y=y, family="binomial",
                    alpha=1, nfolds=nfolds
                )
                P$cvm[row] <- cvm$cvm[cvm$lambda == cvm$lambda.min]
                row <- row + 1
            }
        }
    }

    # Select best parameters and fit final model
    logf("Selecting best parameters")
    best_idx <- which.min(P$cvm)
    best <- as.list(P[best_idx, ])
    logf("Best: npmax=%d, maxShift=%d, maxSnap=%.1f, cvm=%.4f",
        best$npmax, best$maxShift, best$maxSnap, best$cvm)

    mdm <- fit_mdm(
        spectra = spectra, y = y, sfr = sfr,
        nworkers = nworkers, verbose = verbose,
        use_rust = use_rust, nfolds = nfolds,
        cadir = cadir, npmax = best$npmax,
        maxShift = best$maxShift, maxSnap = best$maxSnap
    )
    mdm$pgrid <- P
    mdm$best <- best
    logf("Finished cv_mdm")
    mdm
}

#' @export
#' @title Estimate MDM Performance via Nested Cross-Validation
#' @description
#' Splits the data into `kout` outer folds. For each fold, calls
#' [metabodecon::cv_mdm()] on the training subset to select the best preprocessing
#' parameters and fit a model, then evaluates on the held-out test fold.
#' Grid search results are shared across all folds via a temporary cache
#' directory.
#'
#' @inheritParams cv_mdm
#' @param kout Number of outer folds. Default 10.
#' @param seed Random seed for fold assignment.
#'
#' @return A data frame with columns `fold`, `true`, `prob`, `pred`.
#'
#' @examples
#' \dontrun{
#'      aki <- read_aki_data()
#'      spectra <- aki$spectra
#'      y <- factor(aki$meta$type, levels = c("Control", "AKI"))
#'      perf <- estimate_mdm_performance(spectra, y, sfr = c(11, -2),
#'          nworkers = 53, kout = 10, seed = 1)
#'      mean(perf$true == perf$pred)
#' }
benchmark_cv_mdm <- function(
    spectra, y, sfr, nworkers,
    verbose = TRUE, use_rust = TRUE,
    nfolds = 10, kout = 10, seed = 1
) {

    # Shared temp cache for all folds
    cd <- tempfile("grid_cache_")

    # Pre-warm cache once for all spectra
    logf("Pre-warming grid cache for %d spectra", length(spectra))
    grid_deconvolute_spectra(
        x = spectra, sfr = sfr, verbose = verbose,
        nw = nworkers, use_rust = use_rust, cachedir = cd
    )

    # Assign folds
    set.seed(seed)
    n <- length(spectra)
    foldid <- sample(rep(seq_len(kout), length.out = n))

    # Run outer CV
    results <- vector("list", kout)
    for (f in seq_len(kout)) {
        logf("Outer fold %d/%d", f, kout)
        tr <- which(foldid != f)
        te <- which(foldid == f)

        mdm <- cv_mdm(
            spectra = spectra[tr], y = y[tr], sfr = sfr,
            nworkers = nworkers, verbose = verbose,
            use_rust = use_rust, nfolds = nfolds, cadir = cd
        )

        # Deconvolute + align all spectra with best params
        all_decons <- deconvolute(
            x = spectra, sfr = sfr, verbose = verbose,
            use_rust = use_rust, npmax = mdm$best$npmax,
            nworkers = nworkers, cachedir = cd
        )
        all_als <- align(
            all_decons,
            maxShift = mdm$best$maxShift,
            maxCombine = 0,
            verbose = verbose,
            nworkers = nworkers
        )
        Xall <- t(get_si_mat(all_als, maxSnap = mdm$best$maxSnap))
        Xte <- Xall[te, , drop = FALSE]

        prob <- as.numeric(stats::predict(
            mdm$model, newx = Xte,
            s = "lambda.min", type = "response"
        ))
        pred <- factor(
            ifelse(prob > 0.5, levels(y)[2], levels(y)[1]),
            levels = levels(y)
        )
        results[[f]] <- data.frame(
            fold = f, true = y[te], prob = prob, pred = pred
        )
    }

    Y <- do.call(rbind, results)
    acc <- mean(Y$true == Y$pred)
    logf("Nested CV accuracy: %.2f%%", acc * 100)
    Y
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

#' @title Compute mdm data signature
#' @description Creates a stable short signature for matrix/label pairs.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @return Character hash string.
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' y <- rep(0:1, length.out = 5)
#' mdm_data_signature(X, y)
#' @noRd
mdm_data_signature <- function(X, y) {
    y <- mdm_as_binary01(y)
    payload <- list(X = X, y = y)
    txt <- serialize(payload, connection = NULL, ascii = TRUE, version = 2)
    paper_aki_cache_hash(rawToChar(txt))
}

# S3 methods #####

#' @title Predict from mdm model
#' @description Predicts probabilities, classes, or link scores from an `mdm` object.
#' @param object Fitted `mdm` object.
#' @param newdata Numeric feature matrix.
#' @param type Prediction type, one of `"prob"`, `"class"`, `"link"`.
#' @param s Regularization value for lasso predictions.
#' @param ... Additional arguments, currently unused.
#' @return Numeric vector of probabilities, classes, or link scores.
#' @examples
#' \dontrun{
#'   m <- cv_mdm(spectra, y, sfr = c(11, -2))
#'   predict(m, X, type = "prob")
#' }
predict.mdm <- function(object,
                        newdata,
                        type = c("prob", "class", "link"),
                        s = "lambda.min",
                        ...) {
    type <- match.arg(type)
    norm <- list(method = object$normalisation, ref = object$ref)
    Xn <- mdm_apply_normalizer(newdata, norm)

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
