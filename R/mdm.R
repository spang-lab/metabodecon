#' @title Metabodecon model module
#' @description Helper and API functions for model fitting, prediction, and
#' benchmark evaluation in the AKI reproduction workflow.
#'
#' + API functions:
#'   - fit_mdm: Fit an mdm object with svm or lasso and selected normalization.
#'   - mdm_score_grid: Score one train/test split for a hyperparameter setting grid.
#'   - mdm_find_best_params: Select best hyperparameters by inner CV AUC.
#'   - benchmark_cv_mdm: Run nested CV benchmark with cv-wise or global normalization timing.
#'   - benchmark_extension_table: Build summary table over model and normalization combinations.
#' + Helper functions:
#'   - mdm_auc: Compute AUC from binary labels and numeric scores using rank sums.
#'   - mdm_select_top_features: Rank features by Welch-like t-score and return top indices.
#'   - mdm_get_quantile_reference: Build quantile-normalization reference from sample-wise sorted values.
#'   - mdm_quantile_normalize: Apply quantile normalization to each sample row with a shared reference.
#'   - mdm_get_pqn_reference: Build PQN reference as median profile across samples.
#'   - mdm_pqn_normalize: Apply probabilistic quotient normalization using a reference profile.
#'   - mdm_get_creatinine_ref: Resolve creatinine reference bin index from metadata or fallback.
#'   - mdm_fit_normalizer: Fit normalization metadata used during training and prediction.
#'   - mdm_apply_normalizer: Apply selected normalization to a feature matrix.
#'   - mdm_fmt_mean_sd: Format mean ± sd strings for result tables.
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

mdm_as_binary01 <- function(y, arg_name = "y") {
    as_binary01(y, arg_name)
}

#' @title Fit metabodecon model
#' @description Fits an `mdm` object using SVM or lasso with selected normalization.
#' @param X Numeric feature matrix with samples in rows.
#' @param y Binary labels (0/1).
#' @param model Model type, either `"svm"` or `"lasso"`.
#' @param normalisation Normalization method name.
#' @param cost SVM cost parameter.
#' @param gamma SVM gamma parameter.
#' @param nfeat Number of selected features for SVM.
#' @param cre_idx Optional index for creatinine normalization.
#' @param nfolds_lasso Number of folds for `glmnet::cv.glmnet`.
#' @return Object of class `mdm`.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- fit_mdm(X, y, model = "svm", normalisation = "none", nfeat = 3)
#'   print(m)
#' }
fit_mdm <- function(X,
                    y,
                    model = c("svm", "lasso"),
                    normalisation = c(
                        "quantile", "pqn", "creatinine",
                        "sum", "median", "none"
                    ),
                    cost = 0.5,
                    gamma = 2^-10,
                    nfeat = 3,
                    cre_idx = NULL,
                    nfolds_lasso = 5) {
    model <- match.arg(model)
    normalisation <- match.arg(normalisation)
    y <- mdm_as_binary01(y)
    if (length(unique(y)) < 2) stop("training set must contain both classes")

    norm <- mdm_fit_normalizer(X, normalisation, cre_idx)
    Xn <- mdm_apply_normalizer(X, norm)

    meta <- list(model = model)
    if (model == "svm") {
        if (!requireNamespace("e1071", quietly = TRUE)) {
            stop("Package 'e1071' is required for svm")
        }
        idx <- mdm_select_top_features(Xn, y, nfeat)
        fit <- e1071::svm(
            x = Xn[, idx, drop = FALSE],
            y = factor(y, levels = c(0, 1)),
            type = "C-classification",
            kernel = "radial",
            cost = cost,
            gamma = gamma,
            probability = FALSE
        )
        pred <- predict(fit, Xn[, idx, drop = FALSE], decision.values = TRUE)
        score <- as.numeric(attr(pred, "decision.values"))
        m1 <- mean(score[y == 1], na.rm = TRUE)
        m0 <- mean(score[y == 0], na.rm = TRUE)
        flip <- is.finite(m1) && is.finite(m0) && m1 < m0
        meta$idx <- idx
        meta$flip <- flip
        meta$cost <- cost
        meta$gamma <- gamma
        meta$nfeat <- nfeat
        obj <- fit
    } else {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
            stop("Package 'glmnet' is required for lasso")
        }
        fit <- glmnet::cv.glmnet(
            x = Xn,
            y = y,
            family = "binomial",
            alpha = 1,
            nfolds = min(nfolds_lasso, nrow(Xn))
        )
        obj <- fit
    }

    structure(
        list(
            model = obj,
            normalisation = normalisation,
            ref = norm$ref,
            meta = meta
        ),
        class = "mdm"
    )
}

#' @title Score hyperparameter grid
#' @description Evaluates one train/test split for SVM or lasso.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param tr Integer indices for training rows.
#' @param te Integer indices for test rows.
#' @param model Model type.
#' @param normalisation Normalization method name.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats SVM feature-count candidates.
#' @return Data frame with columns `cost`, `gamma`, `nfeat`, `auc`, `acc`.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(100), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   tr <- 1:14
#'   te <- 15:20
#'   mdm_score_grid(X, y, tr, te, model = "svm", normalisation = "none")
#' }
#' @noRd
mdm_score_grid <- function(X,
                           y,
                           tr,
                           te,
                           model = c("svm", "lasso"),
                           normalisation = c(
                               "quantile", "pqn", "creatinine",
                               "sum", "median", "none"
                           ),
                           costs = c(0.2, 0.5, 0.8),
                           gammas = 2^c(-11, -10, -9),
                           nfeats = 2:3) {
    model <- match.arg(model)
    normalisation <- match.arg(normalisation)
    y <- mdm_as_binary01(y)
    if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) {
        return(data.frame(
            cost = NA, gamma = NA, nfeat = NA,
            auc = NA_real_, acc = NA_real_
        ))
    }

    if (model == "lasso") {
        m <- fit_mdm(X[tr, , drop = FALSE], y[tr], model, normalisation)
        p <- predict(m, X[te, , drop = FALSE], type = "prob")
        c <- predict(m, X[te, , drop = FALSE], type = "class")
        return(data.frame(
            cost = NA, gamma = NA, nfeat = NA,
            auc = mdm_auc(y[te], p),
            acc = mean(y[te] == c)
        ))
    }

    G <- expand.grid(cost = costs, gamma = gammas, nfeat = nfeats)
    G$auc <- NA_real_
    G$acc <- NA_real_
    for (i in seq_len(nrow(G))) {
        m <- fit_mdm(
            X[tr, , drop = FALSE],
            y[tr],
            model = "svm",
            normalisation = normalisation,
            cost = G$cost[i],
            gamma = G$gamma[i],
            nfeat = G$nfeat[i]
        )
        p <- predict(m, X[te, , drop = FALSE], type = "prob")
        c <- predict(m, X[te, , drop = FALSE], type = "class")
        G$auc[i] <- mdm_auc(y[te], p)
        G$acc[i] <- mean(y[te] == c)
    }
    G
}

#' @title Select best parameters by CV
#' @description Runs inner CV and chooses SVM parameters with maximal mean AUC.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param model Model type.
#' @param normalisation Normalization method name.
#' @param k Number of folds.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats SVM feature-count values.
#' @return List with chosen parameters and scored grid.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(120), nrow = 24)
#'   y <- rep(0:1, each = 12)
#'   mdm_find_best_params(X, y, model = "svm", normalisation = "none", k = 3)
#' }
#' @noRd
mdm_find_best_params <- function(X,
                                 y,
                                 model = c("svm", "lasso"),
                                 normalisation = c(
                                     "quantile", "pqn",
                                     "creatinine", "sum",
                                     "median", "none"
                                 ),
                                 k = 3,
                                 costs = c(0.2, 0.5, 0.8),
                                 gammas = 2^c(-11, -10, -9),
                                 nfeats = 2:3) {
    model <- match.arg(model)
    normalisation <- match.arg(normalisation)
    y <- mdm_as_binary01(y)

    te_ids <- get_test_ids(nfolds = k, nsamples = nrow(X), y = y)

    if (model == "lasso") {
        A <- sapply(te_ids, function(te) {
            tr <- setdiff(seq_len(nrow(X)), te)
            Gi <- mdm_score_grid(X, y, tr, te, model, normalisation)
            Gi$auc
        })
        return(list(
            cost = NA, gamma = NA, nfeat = NA,
            grid = data.frame(
                cost = NA, gamma = NA, nfeat = NA,
                auc = mean(A, na.rm = TRUE),
                acc = NA_real_
            )
        ))
    }

    G <- expand.grid(cost = costs, gamma = gammas, nfeat = nfeats)
    A <- matrix(NA_real_, nrow = nrow(G), ncol = length(te_ids))
    for (i in seq_along(te_ids)) {
        te <- te_ids[[i]]
        tr <- setdiff(seq_len(nrow(X)), te)
        Gi <- mdm_score_grid(
            X, y, tr, te, model, normalisation,
            costs, gammas, nfeats
        )
        A[, i] <- Gi$auc
    }
    G$auc <- rowMeans(A, na.rm = TRUE)
    G$auc[!is.finite(G$auc)] <- -Inf
    b <- G[which.max(G$auc), ]
    list(cost = b$cost, gamma = b$gamma, nfeat = b$nfeat, grid = G)
}

#' @title Benchmark nested CV for mdm
#' @description Evaluates one model/normalization setup under nested CV.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param model Model type.
#' @param normalisation Normalization method.
#' @param k Number of outer and inner folds.
#' @param norm_time Either `"cv"` or `"global"`.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats SVM feature-count values.
#' @return List with fold metrics and summary statistics.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(120), nrow = 24)
#'   y <- rep(0:1, each = 12)
#'   benchmark_cv_mdm(X, y, model = "svm", normalisation = "none", k = 3)
#' }
#' @noRd
benchmark_cv_mdm <- function(X,
                             y,
                             model = c("svm", "lasso"),
                             normalisation = c(
                                 "quantile", "pqn",
                                 "creatinine", "sum",
                                 "median", "none"
                             ),
                             k = 3,
                             norm_time = c("cv", "global"),
                             costs = c(0.2, 0.5, 0.8),
                             gammas = 2^c(-11, -10, -9),
                             nfeats = 2:3) {
    model <- match.arg(model)
    normalisation <- match.arg(normalisation)
    norm_time <- match.arg(norm_time)
    y <- mdm_as_binary01(y)

    key <- paper_aki_cache_key(
        "benchmark_cv_mdm",
        mdm_data_signature(X, y),
        model,
        normalisation,
        k,
        norm_time,
        costs,
        gammas,
        nfeats
    )
    cache <- paper_aki_get_cache_dir()
    hit <- paper_aki_cache_get(cache, key)
    if (!is.null(hit)) {
        return(hit)
    }

    X_use <- X
    norm_cv <- normalisation
    if (norm_time == "global" && normalisation != "none") {
        nobj <- mdm_fit_normalizer(X, normalisation)
        X_use <- mdm_apply_normalizer(X, nobj)
        norm_cv <- "none"
    }

    te_ids <- get_test_ids(nfolds = k, nsamples = nrow(X_use), y = y)
    rows <- lapply(te_ids, function(te) {
        tr <- setdiff(seq_len(nrow(X_use)), te)
        p <- mdm_find_best_params(
            X_use[tr, , drop = FALSE], y[tr],
            model, norm_cv, k,
            costs, gammas, nfeats
        )
        cost <- if (is.null(p$cost) || is.na(p$cost)) 0.5 else p$cost
        gamma <- if (is.null(p$gamma) || is.na(p$gamma)) 2^-10 else p$gamma
        nfeat <- if (is.null(p$nfeat) || is.na(p$nfeat)) 3 else p$nfeat
        m <- fit_mdm(X_use[tr, , drop = FALSE], y[tr],
            model = model,
            normalisation = norm_cv,
            cost = cost,
            gamma = gamma,
            nfeat = nfeat
        )
        prob <- predict(m, X_use[te, , drop = FALSE], type = "prob")
        cls <- predict(m, X_use[te, , drop = FALSE], type = "class")
        c(
            auc = mdm_auc(y[te], prob),
            acc = mean(y[te] == cls)
        )
    })

    M <- do.call(rbind, rows)
    out <- data.frame(
        mean_auc = mean(M[, "auc"], na.rm = TRUE),
        sd_auc = stats::sd(M[, "auc"], na.rm = TRUE),
        mean_acc = mean(M[, "acc"], na.rm = TRUE),
        sd_acc = stats::sd(M[, "acc"], na.rm = TRUE)
    )
    paper_aki_cache_set(cache, key, list(folds = M, summary = out))
}

#' @title Build extension benchmark table
#' @description Compares model and normalization pairs for CV-wise vs global normalization.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param models Character vector of model names.
#' @param norms Character vector of normalization names.
#' @param k Number of folds.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats SVM feature-count values.
#' @return Data frame with formatted AUC summaries.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(120), nrow = 24)
#'   y <- rep(0:1, each = 12)
#'   benchmark_extension_table(X, y, models = "svm", norms = c("none", "sum"))
#' }
#' @noRd
benchmark_extension_table <- function(X,
                                      y,
                                      models = c("svm", "lasso"),
                                      norms = c(
                                          "quantile", "pqn",
                                          "creatinine", "sum",
                                          "median", "none"
                                      ),
                                      k = 3,
                                      costs = c(0.2, 0.5, 0.8),
                                      gammas = 2^c(-11, -10, -9),
                                      nfeats = 2:3) {
    key <- paper_aki_cache_key(
        "benchmark_extension_table",
        mdm_data_signature(X, y),
        models,
        norms,
        k,
        costs,
        gammas,
        nfeats
    )
    cache <- paper_aki_get_cache_dir()
    hit <- paper_aki_cache_get(cache, key)
    if (!is.null(hit)) {
        return(hit)
    }

    out <- lapply(models, function(mdl) {
        lapply(norms, function(nrm) {
            cv <- benchmark_cv_mdm(
                X, y, mdl, nrm, k, "cv",
                costs, gammas, nfeats
            )
            gl <- benchmark_cv_mdm(
                X, y, mdl, nrm, k, "global",
                costs, gammas, nfeats
            )
            data.frame(
                model = mdl,
                normalisation = nrm,
                auc_cv_norm = mdm_fmt_mean_sd(
                    cv$summary$mean_auc,
                    cv$summary$sd_auc
                ),
                auc_global_norm = mdm_fmt_mean_sd(
                    gl$summary$mean_auc,
                    gl$summary$sd_auc
                )
            )
        })
    })
    tab <- do.call(rbind, unlist(out, recursive = FALSE))
    paper_aki_cache_set(cache, key, tab)
}

# Helpers #####

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
    sprintf(paste0("%.", d, "f ± %.", d, "f"), x, y)
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
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- fit_mdm(X, y, model = "svm", normalisation = "none")
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
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- fit_mdm(X, y, model = "svm", normalisation = "none")
#'   print(m)
#' }
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
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(100), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- fit_mdm(X, y, model = "lasso", normalisation = "none")
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
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- fit_mdm(X, y, model = "svm", normalisation = "none")
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
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- fit_mdm(X, y, model = "svm", normalisation = "none")
#'   summary(m)
#' }
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
