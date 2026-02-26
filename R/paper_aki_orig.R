#' @title AKI original helper module
#' @description Helper functions for AKI vignette preprocessing, modeling,
#' caching, fold generation, plotting, and lightweight parallel utilities.
#'
#' - bin_spectra: Bin spectra into fixed ppm bins and return feature matrix.
#' - get_quantile_reference: Build quantile normalization reference profile.
#' - quantile_normalize: Apply quantile normalization to feature matrix rows.
#' - get_pqn_reference: Build PQN reference spectrum from median feature values.
#' - pqn_normalize: Apply probabilistic quotient normalization.
#' - auc: Compute rank-based AUC for binary labels and numeric scores.
#' - select_top_features: Select features with highest two-class t-like scores.
#' - train_svm: Train radial SVM with optional quantile normalization.
#' - train_model: Train and cache radial SVM for repeated CV calls.
#' - predict_svm: Predict probabilities and classes from trained SVM wrapper.
#' - score_grid: Evaluate one split across SVM hyperparameter grid.
#' - find_best_params: Select best SVM hyperparameters by inner CV AUC.
#' - score_grid_cached: Cache score_grid outputs to disk as RDS files.
#' - find_best_params_cached: Cache find_best_params outputs to disk as RDS.
#' - paper_aki_cache_new: Create or ensure a cache directory exists.
#' - paper_aki_get_cache_dir: Get active cache directory from options.
#' - paper_aki_cache_key: Build string key from function inputs.
#' - paper_aki_cache_hash: Build stable short hash from cache key.
#' - paper_aki_cache_file: Resolve cache file path for a key.
#' - paper_aki_cache_get: Read cached R object from disk when available.
#' - paper_aki_cache_set: Save R object to disk cache as RDS.
#' - paper_aki_cache_size: Report total size of cached RDS files.
#' - paper_aki_cache_clear: Remove cached RDS files and return delete count.
#' - paper_aki_ncores: Compute worker count from available CPU cores.
#' - paper_aki_grid_apply: Apply function over grid with optional multicore.
#' - paper_aki_should_run_full_vignette: Decide if expensive vignette chunks run.
#' - get_test_ids: Create (stratified) fold test indices.
#' - plot_folds: Visualize outer or inner fold structure.
#' - plot_auc_grid: Plot AUC landscape across cost and gamma.
#' - plot_prob_scatter: Plot predicted probabilities against sample index.
#' @noRd
NULL

# Helpers #####

#' @title Bin spectra into features
#' @description Converts a list of spectra into a binned feature matrix.
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param regions Two-column matrix with ppm regions to include.
#' @param binwidth Bin width in ppm.
#' @return Numeric matrix with samples in rows and bins in columns.
#' @examples
#' sp <- list(
#'   list(cs = seq(9.5, 0.5, length.out = 100), si = rnorm(100)),
#'   list(cs = seq(9.5, 0.5, length.out = 100), si = rnorm(100))
#' )
#' X <- bin_spectra(sp)
#' dim(X)
#' @noRd
bin_spectra <- function(spectra,
                        regions = rbind(c(9.5, 6.5), c(4.5, 0.5)),
                        binwidth = 0.01) {
    edges <- unlist(lapply(seq_len(nrow(regions)), function(i) {
        seq(regions[i, 1], regions[i, 2], by = -binwidth)
    }))
    nb <- length(edges) - 1
    cen <- (edges[-1] + edges[-length(edges)]) / 2
    X <- matrix(0, nrow = length(spectra), ncol = nb)
    for (i in seq_along(spectra)) {
        cs <- spectra[[i]]$cs
        si <- spectra[[i]]$si
        for (j in seq_len(nb)) {
            hi <- edges[j]
            lo <- edges[j + 1]
            X[i, j] <- sum(si[cs <= hi & cs > lo])
        }
    }
    attr(X, "bin_centers") <- cen
    X
}

#' @title Build quantile reference
#' @description Computes mean sorted profile used for quantile normalization.
#' @param X Numeric feature matrix.
#' @return Numeric reference vector.
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' get_quantile_reference(X)
#' @noRd
get_quantile_reference <- function(X) {
    mean_sorted <- apply(X, 1, sort)
    rowMeans(mean_sorted)
}

#' @title Quantile-normalize matrix rows
#' @description Replaces each row with sorted reference values by rank.
#' @param X Numeric feature matrix.
#' @param ref Numeric quantile reference.
#' @return Numeric matrix with quantile-normalized rows.
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' ref <- get_quantile_reference(X)
#' quantile_normalize(X, ref)
#' @noRd
quantile_normalize <- function(X, ref) {
    Xn <- X
    for (i in seq_len(nrow(X))) {
        o <- order(X[i, ])
        Xn[i, o] <- ref
    }
    Xn
}

#' @title Build PQN reference spectrum
#' @description Computes median spectrum across samples for PQN.
#' @param X Numeric feature matrix.
#' @return Numeric reference spectrum.
#' @examples
#' X <- matrix(abs(rnorm(24)), nrow = 6)
#' get_pqn_reference(X)
#' @noRd
get_pqn_reference <- function(X) {
    apply(X, 2, median)
}

#' @title Apply PQN normalization
#' @description Normalizes rows by median quotient to a reference spectrum.
#' @param X Numeric feature matrix.
#' @param ref Numeric PQN reference spectrum.
#' @return PQN-normalized matrix.
#' @examples
#' X <- matrix(abs(rnorm(24)), nrow = 6)
#' ref <- get_pqn_reference(X)
#' pqn_normalize(X, ref)
#' @noRd
pqn_normalize <- function(X, ref) {
    Q <- sweep(X, 2, ref, "/")
    s <- apply(Q, 1, median)
    sweep(X, 1, s, "/")
}

#' @title Compute AUC
#' @description Computes rank-based AUC for binary labels.
#' @param y Binary labels coded as 0/1.
#' @param score Numeric scores.
#' @return Numeric scalar AUC or `NA_real_` when class is missing.
#' @examples
#' auc(c(0, 0, 1, 1), c(0.1, 0.2, 0.8, 0.9))
#' @noRd
auc <- function(y, score) {
    y <- as.integer(y)
    pos <- y == 1
    n1 <- sum(pos)
    n0 <- sum(!pos)
    if (n1 == 0 || n0 == 0) {
        return(NA_real_)
    }
    r <- rank(score)
    (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @title Select top scoring features
#' @description Uses a t-like score to rank features between two classes.
#' @param X Numeric feature matrix.
#' @param y Binary labels coded as 0/1.
#' @param nfeat Number of features to select.
#' @return Integer vector of selected column indices.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40), nrow = 8)
#' y <- rep(0:1, each = 4)
#' select_top_features(X, y, 3)
#' @noRd
select_top_features <- function(X, y, nfeat) {
    i1 <- which(y == 1)
    i0 <- which(y == 0)
    m1 <- colMeans(X[i1, , drop = FALSE])
    m0 <- colMeans(X[i0, , drop = FALSE])
    v1 <- apply(X[i1, , drop = FALSE], 2, var)
    v0 <- apply(X[i0, , drop = FALSE], 2, var)
    s <- sqrt(v1 / length(i1) + v0 / length(i0))
    t <- abs((m1 - m0) / s)
    t[!is.finite(t)] <- 0
    order(t, decreasing = TRUE)[seq_len(nfeat)]
}

#' @title Build fold test indices
#' @description Creates random folds, stratified when labels are supplied.
#' @param nfolds Number of folds.
#' @param nsamples Number of total samples.
#' @param seed RNG seed.
#' @param y Optional label vector for stratification.
#' @return List of integer vectors with test indices per fold.
#' @examples
#' get_test_ids(nfolds = 3, nsamples = 12, seed = 1)
#' y <- rep(0:1, each = 6)
#' get_test_ids(nfolds = 3, nsamples = 12, seed = 1, y = y)
#' @noRd
get_test_ids <- function(nfolds = 5, nsamples, seed = 1, y = NULL) {
    set.seed(seed)
    if (is.null(y)) {
        ids <- sample(seq_len(nsamples))
        grp <- split(ids, cut(seq_along(ids), nfolds, labels = FALSE))
        return(lapply(grp, sort))
    }

    y <- as.integer(y)
    levs <- sort(unique(y))
    out <- vector("list", nfolds)
    for (k in seq_len(nfolds)) out[[k]] <- integer(0)

    for (lev in levs) {
        ids <- sample(which(y == lev))
        grp <- split(ids, cut(seq_along(ids), nfolds, labels = FALSE))
        for (k in seq_len(nfolds)) {
            out[[k]] <- c(out[[k]], grp[[k]])
        }
    }

    lapply(out, sort)
}

# Model #####

#' @title Train radial SVM
#' @description Trains SVM with optional quantile normalization and feature selection.
#' @param X Numeric feature matrix.
#' @param y Binary labels coded as 0/1.
#' @param cost SVM cost.
#' @param gamma SVM gamma.
#' @param nfeat Number of selected features.
#' @param normalize Logical; apply quantile normalization before fitting.
#' @return List with model fit, normalization ref, selected features, and metadata.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- train_svm(X, y, cost = 0.5, gamma = 2^-10, nfeat = 3)
#'   names(m)
#' }
#' @noRd
train_svm <- function(X,
                      y,
                      cost,
                      gamma,
                      nfeat,
                      normalize = TRUE) {
    if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' is required.")
    }
    if (length(unique(y)) < 2) stop("training set must contain both classes")
    ref <- NULL
    Xn <- X
    if (normalize) {
        ref <- get_quantile_reference(X)
        Xn <- quantile_normalize(X, ref)
    }
    idx <- select_top_features(Xn, y, nfeat)
    fit <- e1071::svm(
        x = Xn[, idx, drop = FALSE],
        y = factor(y, levels = c(0, 1)),
        type = "C-classification",
        kernel = "radial",
        cost = cost,
        gamma = gamma,
        probability = FALSE
    )
    list(
        fit = fit,
        ref = ref,
        idx = idx,
        normalize = normalize,
        y_train = y
    )
}

#' @title Train cached model
#' @description Trains and caches radial SVM models for repeated CV calls.
#' @param X Numeric feature matrix.
#' @param y Binary labels coded as 0/1.
#' @param cost SVM cost.
#' @param gamma SVM gamma.
#' @param nfeat Number of selected features.
#' @param normalize Logical; apply quantile normalization before fitting.
#' @return List with model fit, normalization ref, selected features, metadata.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- train_model(X, y, cost = 0.5, gamma = 2^-10, nfeat = 3)
#'   names(m)
#' }
#' @noRd
train_model <- function(X,
                        y,
                        cost,
                        gamma,
                        nfeat,
                        normalize = TRUE) {
    payload <- list(X = X, y = y)
    txt <- serialize(payload, connection = NULL, ascii = TRUE, version = 2)
    sig <- paper_aki_cache_hash(rawToChar(txt))
    key <- paper_aki_cache_key("train_model", sig, cost, gamma, nfeat, normalize)
    cache <- paper_aki_get_cache_dir()
    hit <- paper_aki_cache_get(cache, key)
    if (!is.null(hit)) {
        return(hit)
    }
    out <- train_svm(X, y, cost, gamma, nfeat, normalize)
    paper_aki_cache_set(cache, key, out)
}

#' @title Predict with SVM wrapper
#' @description Predicts class probabilities and classes from `train_svm` output.
#' @param model Model list created by `train_svm`.
#' @param X Numeric feature matrix.
#' @return List with `prob` and `cls`.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   m <- train_svm(X, y, cost = 0.5, gamma = 2^-10, nfeat = 3)
#'   p <- predict_svm(m, X)
#'   names(p)
#' }
#' @noRd
predict_svm <- function(model, X) {
    if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("Package 'e1071' is required.")
    }
    Xn <- X
    if (model$normalize) Xn <- quantile_normalize(Xn, model$ref)
    pred <- stats::predict(
        model$fit,
        Xn[, model$idx, drop = FALSE],
        decision.values = TRUE
    )
    score <- as.numeric(attr(pred, "decision.values"))
    y0 <- as.integer(model$y_train)
    m1 <- mean(score[y0 == 1], na.rm = TRUE)
    m0 <- mean(score[y0 == 0], na.rm = TRUE)
    if (is.finite(m1) && is.finite(m0) && m1 < m0) score <- -score
    score[!is.finite(score)] <- 0
    p1 <- 1 / (1 + exp(-score))
    list(prob = p1, cls = as.integer(as.character(pred)))
}

#' @title Score one hyperparameter grid
#' @description Computes AUC/accuracy for each grid point on one split.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param tr Training indices.
#' @param te Test indices.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats Candidate selected-feature counts.
#' @param normalize Logical; apply quantile normalization in fitting.
#' @return Data frame with one row per parameter combination.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(100), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   g <- score_grid(X, y, 1:14, 15:20, c(0.5), c(2^-10), 3, TRUE)
#'   g[, c("cost", "gamma", "nfeat", "auc", "acc")]
#' }
#' @noRd
score_grid <- function(X,
                       y,
                       tr,
                       te,
                       costs,
                       gammas,
                       nfeats,
                       normalize = TRUE) {
    G <- expand.grid(cost = costs, gamma = gammas, nfeat = nfeats)
    G$auc <- NA_real_
    G$acc <- NA_real_
    if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) {
        return(G)
    }
    for (i in seq_len(nrow(G))) {
        f <- train_model(
            X[tr, , drop = FALSE],
            y[tr],
            cost = G$cost[i],
            gamma = G$gamma[i],
            nfeat = G$nfeat[i],
            normalize = normalize
        )
        p <- predict_svm(f, X[te, , drop = FALSE])
        G$auc[i] <- auc(y[te], p$prob)
        G$acc[i] <- mean(y[te] == p$cls)
    }
    G
}

#' @title Find best SVM parameters
#' @description Uses k-fold CV and chooses the setting with highest mean AUC.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats Candidate selected-feature counts.
#' @param normalize Logical; apply quantile normalization.
#' @param k Number of CV folds.
#' @return List with `cost`, `gamma`, `nfeat`, and scored `grid`.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(120), nrow = 24)
#'   y <- rep(0:1, each = 12)
#'   find_best_params(X, y, costs = c(0.5), gammas = c(2^-10), nfeats = 3, k = 3)
#' }
#' @noRd
find_best_params <- function(X,
                             y,
                             costs = (1:10) / 10,
                             gammas = 2^(-12:-5),
                             nfeats = 2:3,
                             normalize = TRUE,
                             k = 5) {
    te_ids <- get_test_ids(k, nrow(X), y = y)
    G <- expand.grid(cost = costs, gamma = gammas, nfeat = nfeats)
    A <- matrix(NA_real_, nrow = nrow(G), ncol = length(te_ids))
    for (i in seq_along(te_ids)) {
        te <- te_ids[[i]]
        tr <- setdiff(seq_len(nrow(X)), te)
        Gi <- score_grid(X, y, tr, te, costs, gammas, nfeats, normalize)
        A[, i] <- Gi$auc
    }
    G$auc <- rowMeans(A, na.rm = TRUE)
    G$auc[!is.finite(G$auc)] <- -Inf
    b <- G[which.max(G$auc), ]
    list(cost = b$cost, gamma = b$gamma, nfeat = b$nfeat, grid = G)
}

# Cache #####

#' @title Create cache directory
#' @description Ensures a disk cache directory exists and returns its path.
#' @param path Optional directory path.
#' @return Character scalar path to cache directory.
#' @examples
#' d <- paper_aki_cache_new(tempfile("aki-cache-"))
#' dir.exists(d)
#' @noRd
paper_aki_cache_new <- function(path = NULL) {
    if (is.null(path)) {
        path <- file.path(tempdir(), "metabodecon-aki-cache")
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

#' @title Get active cache directory
#' @description Returns cache directory from options and initializes it when needed.
#' @return Character scalar path to cache directory.
#' @examples
#' d <- paper_aki_get_cache_dir()
#' dir.exists(d)
#' @noRd
paper_aki_get_cache_dir <- function() {
    pdir <- file.path(datadir_persistent(), "cache")
    if (dir.exists(pdir)) {
        options(metabodecon.aki_cache = pdir)
        return(pdir)
    }

    cache_dir <- getOption("metabodecon.aki_cache")
    if (is.null(cache_dir) || !is.character(cache_dir) || length(cache_dir) != 1) {
        cache_dir <- paper_aki_cache_new()
        options(metabodecon.aki_cache = cache_dir)
    }
    if (!dir.exists(cache_dir)) {
        cache_dir <- paper_aki_cache_new(cache_dir)
        options(metabodecon.aki_cache = cache_dir)
    }
    cache_dir
}

#' @title Build cache key
#' @description Concatenates inputs into a stable key string.
#' @param ... Objects used to identify a cached result.
#' @return Character key.
#' @examples
#' paper_aki_cache_key("score_grid", 1:3, c(0.2, 0.5))
#' @noRd
paper_aki_cache_key <- function(...) {
    x <- list(...)
    paste(vapply(x, function(z) paste(z, collapse = ","), ""), collapse = "|")
}

#' @title Hash cache key
#' @description Creates a short deterministic hash from a key string.
#' @param key Character cache key.
#' @return Character hash string.
#' @examples
#' paper_aki_cache_hash("abc|123")
#' @noRd
paper_aki_cache_hash <- function(key) {
    r <- utf8ToInt(key)
    h <- 0
    m <- 2147483647
    for (i in seq_along(r)) {
        h <- (h * 131 + r[i]) %% m
    }
    paste0("k", nchar(key), "_", sprintf("%010d", h))
}

#' @title Resolve cache file path
#' @description Builds full RDS file path for cache key in cache directory.
#' @param cache_dir Cache directory path.
#' @param key Cache key.
#' @return Character path to cache file.
#' @examples
#' d <- paper_aki_cache_new(tempfile("aki-cache-"))
#' key <- paper_aki_cache_key("x", 1:3)
#' paper_aki_cache_file(d, key)
#' @noRd
paper_aki_cache_file <- function(cache_dir, key) {
    file.path(cache_dir, paste0(paper_aki_cache_hash(key), ".rds"))
}

#' @title Get cached object from disk
#' @description Reads cached RDS object when available, otherwise returns NULL.
#' @param cache_dir Cache directory path.
#' @param key Cache key.
#' @return Cached object or `NULL`.
#' @examples
#' d <- paper_aki_cache_new(tempfile("aki-cache-"))
#' key <- paper_aki_cache_key("demo")
#' paper_aki_cache_get(d, key)
#' @noRd
paper_aki_cache_get <- function(cache_dir, key) {
    file <- paper_aki_cache_file(cache_dir, key)
    if (!file.exists(file)) {
        return(NULL)
    }
    out <- tryCatch(readRDS(file), error = function(e) NULL)
    out
}

#' @title Store cached object on disk
#' @description Saves object as RDS file under hashed cache key.
#' @param cache_dir Cache directory path.
#' @param key Cache key.
#' @param value Object to cache.
#' @return Invisibly returns `value`.
#' @examples
#' d <- paper_aki_cache_new(tempfile("aki-cache-"))
#' key <- paper_aki_cache_key("demo")
#' paper_aki_cache_set(d, key, list(a = 1))
#' @noRd
paper_aki_cache_set <- function(cache_dir, key, value) {
    file <- paper_aki_cache_file(cache_dir, key)
    saveRDS(value, file)
    value
}

#' @title Compute cache size
#' @description Sums bytes of cached `.rds` files and returns a compact summary.
#' @param cache_dir Cache directory path.
#' @return List with `bytes`, `files`, and `human`.
#' @examples
#' d <- paper_aki_cache_new(tempfile("aki-cache-"))
#' paper_aki_cache_size(d)
#' @noRd
paper_aki_cache_size <- function(cache_dir = paper_aki_get_cache_dir()) {
    files <- list.files(cache_dir,
        pattern = "\\.rds$",
        full.names = TRUE,
        ignore.case = TRUE
    )
    if (length(files) == 0) {
        return(list(bytes = 0, files = 0L, human = "0 B"))
    }
    sz <- sum(file.info(files)$size, na.rm = TRUE)
    units <- c("B", "KB", "MB", "GB", "TB")
    idx <- 1L
    val <- as.numeric(sz)
    while (val >= 1024 && idx < length(units)) {
        val <- val / 1024
        idx <- idx + 1L
    }
    list(
        bytes = as.numeric(sz),
        files = length(files),
        human = sprintf("%.2f %s", val, units[idx])
    )
}

#' @title Clear disk cache
#' @description Deletes cached `.rds` files and returns number of removed files.
#' @param cache_dir Cache directory path.
#' @return Integer count of deleted files.
#' @examples
#' d <- paper_aki_cache_new(tempfile("aki-cache-"))
#' paper_aki_cache_clear(d)
#' @noRd
paper_aki_cache_clear <- function(cache_dir = paper_aki_get_cache_dir()) {
    files <- list.files(cache_dir,
        pattern = "\\.rds$",
        full.names = TRUE,
        ignore.case = TRUE
    )
    if (length(files) == 0) {
        return(0L)
    }
    as.integer(sum(file.remove(files)))
}

#' @title Cached score_grid
#' @description Runs `score_grid` and stores/retrieves results from disk cache.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param tr Training indices.
#' @param te Test indices.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats Candidate selected-feature counts.
#' @param normalize Logical; apply quantile normalization.
#' @return Data frame from `score_grid`.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(80), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   score_grid_cached(X, y, 1:14, 15:20, 0.5, 2^-10, 3, TRUE)
#' }
#' @noRd
score_grid_cached <- function(X,
                              y,
                              tr,
                              te,
                              costs,
                              gammas,
                              nfeats,
                              normalize = TRUE) {
    key <- paper_aki_cache_key(
        "score_grid", tr, te, costs, gammas, nfeats, normalize
    )
    cache <- paper_aki_get_cache_dir()
    hit <- paper_aki_cache_get(cache, key)
    if (!is.null(hit)) {
        return(hit)
    }
    out <- score_grid(X, y, tr, te, costs, gammas, nfeats, normalize)
    paper_aki_cache_set(cache, key, out)
}

#' @title Cached parameter search
#' @description Runs `find_best_params` and stores/retrieves from disk cache.
#' @param X Numeric feature matrix.
#' @param y Binary labels.
#' @param costs SVM cost values.
#' @param gammas SVM gamma values.
#' @param nfeats Candidate selected-feature counts.
#' @param normalize Logical; apply quantile normalization.
#' @param k Number of CV folds.
#' @return List from `find_best_params`.
#' @examples
#' if (requireNamespace("e1071", quietly = TRUE)) {
#'   set.seed(1)
#'   X <- matrix(rnorm(100), nrow = 20)
#'   y <- rep(0:1, each = 10)
#'   find_best_params_cached(X, y, 0.5, 2^-10, 3, TRUE, 3)
#' }
#' @noRd
find_best_params_cached <- function(X,
                                    y,
                                    costs,
                                    gammas,
                                    nfeats,
                                    normalize = TRUE,
                                    k = 3) {
    key <- paper_aki_cache_key(
        "find_best_params", y, costs, gammas, nfeats, normalize, k
    )
    cache <- paper_aki_get_cache_dir()
    hit <- paper_aki_cache_get(cache, key)
    if (!is.null(hit)) {
        return(hit)
    }
    out <- find_best_params(X, y, costs, gammas, nfeats, normalize, k)
    paper_aki_cache_set(cache, key, out)
}

# Runtime #####

#' @title Detect number of worker cores
#' @description Returns core count scaled by `frac` and bounded below by one.
#' @param frac Fraction of detected logical cores.
#' @return Integer worker count.
#' @examples
#' paper_aki_ncores(0.5)
#' @noRd
paper_aki_ncores <- function(frac = 0.5) {
    n <- parallel::detectCores(logical = TRUE)
    max(1L, floor(n * frac))
}

#' @title Apply function over grid
#' @description Evaluates a function for each row of a grid, optionally multicore.
#' @param grid Data frame of parameter settings.
#' @param fun Function taking `(row_df, i)`.
#' @param ncores Number of cores.
#' @return List of results with `ncores` attribute.
#' @examples
#' g <- expand.grid(a = 1:2, b = 3:4)
#' paper_aki_grid_apply(g, function(x, i) c(i = i, s = x$a + x$b), ncores = 1)
#' @noRd
paper_aki_grid_apply <- function(grid,
                                 fun,
                                 ncores = paper_aki_ncores(0.5)) {
    stopifnot(is.function(fun))
    idx <- seq_len(nrow(grid))
    run_one <- function(i) fun(grid[i, , drop = FALSE], i)
    if (.Platform$OS.type == "unix" && ncores > 1) {
        out <- parallel::mclapply(idx, run_one, mc.cores = ncores)
    } else {
        out <- lapply(idx, run_one)
    }
    attr(out, "ncores") <- ncores
    out
}

#' @title Decide full vignette execution
#' @description Returns TRUE on GitHub Actions or when explicit local option/env is set.
#' @return Logical scalar.
#' @examples
#' paper_aki_should_run_full_vignette()
#' @noRd
paper_aki_should_run_full_vignette <- function() {
    is_github <- identical(Sys.getenv("GITHUB_ACTIONS"), "true")
    opt <- isTRUE(getOption("metabodecon.full_vignette", FALSE))
    env <- identical(Sys.getenv("METABODECON_FULL_VIGNETTE"), "true")
    is_github || opt || env
}

# Plots #####

#' @title Plot fold structure
#' @description Plots outer folds or outer+inner validation split layout.
#' @param test_ids List of outer test indices.
#' @param nsamples Total number of samples.
#' @param val_ids Optional list of inner validation indices.
#' @param main Plot title.
#' @return Invisible integer matrix used for plotting.
#' @examples
#' ids <- list(c(1, 4), c(2, 5), c(3, 6))
#' tmp <- tempfile(fileext = ".pdf")
#' grDevices::pdf(tmp)
#' plot_folds(ids, nsamples = 6)
#' grDevices::dev.off()
#' @noRd
plot_folds <- function(test_ids,
                       nsamples,
                       val_ids = NULL,
                       main = "Cross-validation folds") {
    n <- length(test_ids)
    all_ids <- seq_len(nsamples)

    if (is.null(val_ids)) {
        z <- matrix(1L, n, nsamples)
        for (i in seq_len(n)) z[i, test_ids[[i]]] <- 3L
        ylab <- paste0("Fold ", seq_len(n))
    } else {
        m <- length(val_ids)
        z <- matrix(1L, m, nsamples)
        te <- sort(unique(unlist(test_ids)))
        for (i in seq_len(m)) {
            z[i, te] <- 3L
            z[i, val_ids[[i]]] <- 2L
            tr <- setdiff(all_ids, c(te, val_ids[[i]]))
            z[i, tr] <- 1L
        }
        ylab <- paste0("Inner ", seq_len(m))
    }

    cols <- c("#a1d99b", "#9ecae1", "#fdd870")
    withr::with_par(list(mar = c(4, 6, 2, 1)), {
        image(
            x = seq_len(nsamples),
            y = seq_len(nrow(z)),
            z = t(z[rev(seq_len(nrow(z))), ]),
            col = cols,
            axes = FALSE,
            xlab = "Samples",
            ylab = "",
            main = main
        )
        axis(1)
        axis(2, at = seq_len(nrow(z)), labels = rev(ylab), las = 1)
        box()
        legend(
            "topright",
            fill = cols,
            legend = c("train", "validation", "test"),
            cex = 0.85,
            bg = "white"
        )
    })

    invisible(z)
}

#' @title Plot AUC parameter grid
#' @description Plots AUC values over cost/gamma for one selected `nfeat`.
#' @param grid Data frame with `cost`, `gamma`, `nfeat`, `auc` columns.
#' @param nfeat Feature count to visualize.
#' @param zlim Optional numeric limits for color scale.
#' @param main Optional title.
#' @return Invisible matrix of plotted AUC values.
#' @examples
#' g <- expand.grid(cost = c(0.2, 0.5), gamma = c(2^-11, 2^-10), nfeat = 2:3)
#' g$auc <- seq(0.6, 0.9, length.out = nrow(g))
#' tmp <- tempfile(fileext = ".pdf")
#' grDevices::pdf(tmp)
#' plot_auc_grid(g, nfeat = 2)
#' grDevices::dev.off()
#' @noRd
plot_auc_grid <- function(grid,
                          nfeat,
                          zlim = NULL,
                          main = NULL) {
    g <- grid[grid$nfeat == nfeat, c("cost", "gamma", "auc")]
    cs <- sort(unique(g$cost))
    gs <- sort(unique(g$gamma))
    z <- matrix(NA_real_, nrow = length(cs), ncol = length(gs))

    for (i in seq_len(nrow(g))) {
        c <- match(g$cost[i], cs)
        r <- match(g$gamma[i], gs)
        z[c, r] <- g$auc[i]
    }

    if (is.null(zlim)) zlim <- range(grid$auc, na.rm = TRUE)
    withr::with_par(list(mar = c(4, 4, 2, 1)), {
        image(
            x = cs,
            y = gs,
            z = z,
            col = hcl.colors(30, "YlOrRd", rev = FALSE),
            zlim = zlim,
            xlab = "cost",
            ylab = "gamma",
            main = if (is.null(main)) paste("nfeat =", nfeat) else main
        )
        contour(
            x = cs,
            y = gs,
            z = z,
            add = TRUE,
            drawlabels = FALSE,
            nlevels = 8,
            col = "grey25"
        )
    })

    invisible(z)
}

#' @title Plot predicted probabilities
#' @description Draws probability scatter with class-specific colors.
#' @param prob Numeric vector of predicted probabilities.
#' @param y Binary class labels.
#' @param main Plot title.
#' @return `NULL` invisibly.
#' @examples
#' tmp <- tempfile(fileext = ".pdf")
#' grDevices::pdf(tmp)
#' plot_prob_scatter(c(0.1, 0.7, 0.4), c(0, 1, 0), "Example")
#' grDevices::dev.off()
#' @noRd
plot_prob_scatter <- function(prob, y, main = "Predicted probabilities") {
    cols <- ifelse(y == 1, "#d95f02", "#1b9e77")
    withr::with_par(list(mar = c(4, 4, 2, 1)), {
        plot(
            x = seq_along(prob),
            y = prob,
            col = cols,
            pch = 19,
            xlab = "sample",
            ylab = "P(AKI)",
            main = main,
            ylim = c(0, 1)
        )
        abline(h = 0.5, lty = 2)
        legend(
            "topright",
            pch = 19,
            col = c("#1b9e77", "#d95f02"),
            legend = c("no AKI", "AKI"),
            bg = "white",
            cex = 0.85
        )
    })
}
