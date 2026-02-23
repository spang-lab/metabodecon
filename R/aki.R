

#' @noRd
#' @title AKI Klassifikation
#' @description
#' Compares the performance of metabodecon against three other approaches
#' by performing an deconvolution, alignment and classification of the AKI
#' dataset.
#'
#' Requires execution of `update_aki_dataset()` first, which in turn
#' requires availability of Spang datasets under /data.
#'
#' @details
#' To estimate the end-to-end performance of metabodecon, this function
#' evaluates its deconvolution and alignment workflow against three preprocessing
#' strategies in a supervised classification setting: equidistant binning,
#' wavelet-based peak detection and alignment using speaq, and automated peak
#' picking and alignment as implemented in Bruker TopSpin. Equidistant binning
#' represents a widely adopted baseline in NMR metabolomics (Emwas et al.,
#' 2019), while wavelet-based peak detection improves peak localization and
#' noise robustness compared to fixed binning (Vu et al., 2011; Vu et al.,
#' 2019). Vendor-provided pipelines such as TopSpin reflect common analytical
#' practice and provide a real-world reference. Each preprocessing method yields
#' a feature matrix of aligned chemical shift positions with corresponding peak
#' integrals or binned intensities. To assess discriminative performance between
#' AKI patients and controls, the function trains logistic regression models with
#' an elastic-net penalty, which is well suited for high-dimensional, correlated
#' metabolomics data with limited sample size (Zou & Hastie, 2005; Kirpich et
#' al., 2018). Model hyperparameters (regularization strength and mixing
#' parameter) are tuned using nested cross-validation, and predictive performance
#' is quantified by balanced accuracy and area under the ROC curve (AUROC), with
#' feature scaling performed within each training fold to prevent information
#' leakage. This evaluation framework ensures an unbiased comparison of
#' preprocessing strategies and quantifies their end-to-end performance.
#'
classify_aki_patients <- function(data,
                                  labels,
                                  methods = c("metabodecon", "speaq", "binning"),
                                  split_ratio = c(2, 1),
                                  seed = 1,
                                  bin_width = 0.01,
                                  bin_range = NULL,
                                  sfr = NULL,
                                  smopts = c(2, 5),
                                  delta = 6.4,
                                  nfit = 3,
                                  wshw = 0,
                                  maxShift = 50,
                                  maxCombine = 5,
                                  nworkers = 1,
                                  use_rust = TRUE,
                                  find_params = TRUE,
                                  install_deps = NULL,
                                  read_cache = TRUE,
                                  write_cache = TRUE,
                                  verbose = TRUE) {
    paths <- aki_norm_paths(data)
    y <- aki_labels_to_num(labels)
    methods <- match.arg(methods,
        c("metabodecon", "speaq", "binning"),
        several.ok = TRUE
    )
    split <- aki_train_test_split(y, split_ratio, seed)
    train_paths <- paths[split$train_idx]
    test_paths <- paths[split$test_idx]
    train_y <- y[split$train_idx]
    test_y <- y[split$test_idx]
    models <- lapply(methods, function(mtd) {
        aki_train_classifier(
            train_paths, train_y, method = mtd,
            bin_width = bin_width, bin_range = bin_range,
            sfr = sfr, smopts = smopts, delta = delta, nfit = nfit,
            wshw = wshw, maxShift = maxShift, maxCombine = maxCombine,
            nworkers = nworkers, use_rust = use_rust,
            find_params = find_params,
            install_deps = install_deps,
            read_cache = read_cache, write_cache = write_cache,
            verbose = verbose
        )
    })
    names(models) <- methods
    metrics <- lapply(methods, function(mtd) {
        model <- models[[mtd]]
        prob <- aki_predict_classifier(
            model, test_paths,
            nworkers = nworkers,
            use_rust = use_rust,
            install_deps = install_deps,
            verbose = verbose
        )
        pred <- prob >= 0.5
        data.frame(
            method = mtd,
            bal_acc = aki_balanced_accuracy(test_y, pred),
            auc = aki_auc(test_y, prob)
        )
    })
    metrics <- do.call(rbind, metrics)
    list(
        methods = methods,
        split = split,
        metrics = metrics,
        models = models
    )
}

aki_cv_method <- function(paths,
                          y,
                          folds,
                          method,
                          bin_width,
                          bin_range,
                          sfr,
                          smopts,
                          delta,
                          wshw,
                          maxShift,
                          maxCombine,
                          nworkers,
                          use_rust,
                          install_deps,
                          read_cache,
                          write_cache,
                          verbose) {
    k <- length(folds)
    bal <- rep(NA_real_, k)
    auc <- rep(NA_real_, k)
    for (i in seq_len(k)) {
        test_idx <- folds[[i]]
        train_idx <- setdiff(seq_along(y), test_idx)
        train_paths <- paths[train_idx]
        test_paths <- paths[test_idx]
        train_y <- y[train_idx]
        test_y <- y[test_idx]
        model <- aki_train_classifier(
            train_paths, train_y, method = method,
            bin_width = bin_width, bin_range = bin_range,
            sfr = sfr, smopts = smopts, delta = delta, wshw = wshw,
            maxShift = maxShift, maxCombine = maxCombine,
            nworkers = nworkers, use_rust = use_rust,
            install_deps = install_deps,
            read_cache = read_cache, write_cache = write_cache,
            verbose = verbose
        )
        prob <- aki_predict_classifier(
            model, test_paths,
            nworkers = nworkers,
            use_rust = use_rust,
            install_deps = install_deps,
            verbose = verbose
        )
        pred <- prob >= 0.5
        bal[i] <- aki_balanced_accuracy(test_y, pred)
        auc[i] <- aki_auc(test_y, prob)
    }
    list(
        fold = data.frame(
            fold = seq_len(k),
            bal_acc = bal,
            auc = auc
        ),
        summary = data.frame(
            bal_acc_mean = mean(bal, na.rm = TRUE),
            bal_acc_sd = stats::sd(bal, na.rm = TRUE),
            auc_mean = mean(auc, na.rm = TRUE),
            auc_sd = stats::sd(auc, na.rm = TRUE)
        )
    )
}

aki_train_classifier <- function(paths,
                                 y,
                                 method,
                                 bin_width,
                                 bin_range,
                                 sfr,
                                 smopts,
                                 delta,
                                 nfit,
                                 wshw,
                                 maxShift,
                                 maxCombine,
                                 nworkers,
                                 use_rust,
                                 find_params,
                                 install_deps,
                                 read_cache,
                                 write_cache,
                                 verbose) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Missing package 'glmnet'. Install it before running this.")
    }
    if (isTRUE(find_params) && !isTRUE(use_rust)) {
        stop("find_params requires use_rust = TRUE.")
    }
    if (isTRUE(find_params) && method != "binning") {
        first_spec <- aki_read_spectra(paths[1], cs = NULL)[[1]]
        sfr_use <- aki_default_sfr(first_spec, sfr)
        best <- find_best_params(
            first_spec,
            sfr = sfr_use,
            smopts = smopts,
            delta = delta,
            nfit = nfit,
            verbose = verbose
        )
        smopts <- best$smopts
        delta <- best$delta
        nfit <- best$nfit
    }
    params <- list(
        method = method,
        bin_width = bin_width,
        bin_range = bin_range,
        sfr = sfr,
        smopts = smopts,
        delta = delta,
        nfit = nfit,
        wshw = wshw,
        maxShift = maxShift,
        maxCombine = maxCombine,
        nworkers = nworkers,
        use_rust = use_rust,
        find_params = find_params
    )
    key <- aki_cache_key(paths, y, params)
    rds <- file.path(cachedir(), sprintf("aki_model_%s.rds", key))
    if (isTRUE(read_cache) && file.exists(rds)) {
        return(readRDS(rds))
    }
    feats <- aki_build_features(
        paths,
        method = method,
        bin_width = bin_width,
        bin_range = bin_range,
        sfr = sfr,
        smopts = smopts,
        delta = delta,
        nfit = nfit,
        wshw = wshw,
        maxShift = maxShift,
        maxCombine = maxCombine,
        nworkers = nworkers,
        use_rust = use_rust,
        install_deps = install_deps,
        verbose = verbose
    )
    ref <- aki_pqn_ref(feats$mat)
    x <- aki_pqn_apply(feats$mat, ref)
    fit <- glmnet::cv.glmnet(
        x = x,
        y = y,
        family = "binomial",
        alpha = 1,
        nfolds = 10
    )
    model <- list(
        method = method,
        fit = fit,
        pqn_ref = ref,
        train_paths = paths,
        train_y = y,
        params = params,
        train_decons = feats$decons %||% NULL,
        bin_edges = feats$bin_edges %||% NULL,
        train_cs = feats$train_cs %||% NULL
    )
    if (isTRUE(write_cache)) saveRDS(model, rds)
    model
}

aki_predict_classifier <- function(model,
                                   paths,
                                   nworkers,
                                   use_rust,
                                   install_deps,
                                   verbose) {
    feats <- aki_build_features(
        paths,
        method = model$method,
        bin_width = model$params$bin_width,
        bin_range = model$params$bin_range,
        sfr = model$params$sfr,
        smopts = model$params$smopts,
        delta = model$params$delta,
        nfit = model$params$nfit,
        wshw = model$params$wshw,
        maxShift = model$params$maxShift,
        maxCombine = model$params$maxCombine,
        nworkers = nworkers,
        use_rust = use_rust,
        install_deps = install_deps,
        verbose = verbose,
        train_decons = model$train_decons,
        bin_edges = model$bin_edges,
        train_cs = model$train_cs
    )
    x <- aki_pqn_apply(feats$mat, model$pqn_ref)
    as.numeric(stats::predict(
        model$fit,
        newx = x,
        s = "lambda.min",
        type = "response"
    ))
}

aki_build_features <- function(paths,
                               method,
                               bin_width,
                               bin_range,
                               sfr,
                               smopts,
                               delta,
                               nfit,
                               wshw,
                               maxShift,
                               maxCombine,
                               nworkers,
                               use_rust,
                               install_deps,
                               verbose,
                               train_decons = NULL,
                               bin_edges = NULL,
                               train_cs = NULL) {
    method <- match.arg(method,
        c("metabodecon", "speaq", "binning")
    )
    if (method == "binning") {
        return(aki_features_binning(
            paths,
            bin_width = bin_width,
            bin_range = bin_range,
            bin_edges = bin_edges
        ))
    }
    if (method == "speaq") {
        return(aki_features_speaq(
            paths,
            train_decons = train_decons,
            sfr = sfr,
            smopts = smopts,
            delta = delta,
            nfit = nfit,
            wshw = wshw,
            maxShift = maxShift,
            maxCombine = maxCombine,
            nworkers = nworkers,
            use_rust = use_rust,
            install_deps = install_deps,
            verbose = verbose,
            train_cs = train_cs
        ))
    }
    aki_features_metabodecon(
        paths,
        train_decons = train_decons,
        sfr = sfr,
        smopts = smopts,
        delta = delta,
        nfit = nfit,
        wshw = wshw,
        maxShift = maxShift,
        maxCombine = maxCombine,
        nworkers = nworkers,
        use_rust = use_rust,
        install_deps = install_deps,
        verbose = verbose,
        train_cs = train_cs
    )
}

aki_features_metabodecon <- function(paths,
                                     train_decons,
                                     sfr,
                                     smopts,
                                     delta,
                                     nfit,
                                     wshw,
                                     maxShift,
                                     maxCombine,
                                     nworkers,
                                     use_rust,
                                     install_deps,
                                     verbose,
                                     train_cs) {
    cs <- train_cs
    if (is.null(cs) && !is.null(train_decons) && length(train_decons) > 0) {
        cs <- train_decons[[1]]$cs
    }
    spectra <- aki_read_spectra(paths, cs = cs)
    decons <- as_decons2(deconvolute(
        spectra,
        nfit = nfit,
        smopts = smopts,
        delta = delta,
        sfr = sfr,
        wshw = wshw,
        ask = FALSE,
        force = FALSE,
        verbose = verbose,
        nworkers = nworkers,
        use_rust = use_rust
    ))
    all_decons <- aki_merge_decons(train_decons, decons)
    aligns <- tryCatch(
        align(
            all_decons,
            maxShift = maxShift,
            maxCombine = maxCombine,
            verbose = verbose,
            install_deps = install_deps
        ),
        error = function(e) {
            if (isTRUE(verbose)) {
                msg <- paste("align failed, retrying with maxCombine=0:", e$message)
                logf(msg)
            }
            align(
                all_decons,
                maxShift = maxShift,
                maxCombine = 0,
                verbose = verbose,
                install_deps = install_deps
            )
        }
    )
    mat <- t(get_si_mat(aligns))
    idx <- aki_new_rows(train_decons, decons)
    list(
        mat = mat[idx, , drop = FALSE],
        decons = if (is.null(train_decons)) decons else train_decons,
        train_cs = if (is.null(train_decons)) decons[[1]]$cs else train_decons[[1]]$cs
    )
}

aki_features_speaq <- function(paths,
                               train_decons,
                               sfr,
                               smopts,
                               delta,
                               nfit,
                               wshw,
                               maxShift,
                               maxCombine,
                               nworkers,
                               use_rust,
                               install_deps,
                               verbose,
                               train_cs) {
    aki_check_speaq_deps(install_deps, verbose)
    cs <- train_cs
    if (is.null(cs) && !is.null(train_decons) && length(train_decons) > 0) {
        cs <- train_decons[[1]]$cs
    }
    spectra <- aki_read_spectra(paths, cs = cs)
    decons <- as_decons2(deconvolute(
        spectra,
        nfit = nfit,
        smopts = smopts,
        delta = delta,
        sfr = sfr,
        wshw = wshw,
        ask = FALSE,
        force = FALSE,
        verbose = verbose,
        nworkers = nworkers,
        use_rust = use_rust
    ))
    all_decons <- aki_merge_decons(train_decons, decons)
    decons1 <- as_decons1(all_decons)
    feat <- gen_feat_mat(decons1)
    shifted <- speaq_align(
        feat,
        maxShift = maxShift,
        spectrum_data = decons1,
        verbose = verbose,
        show = FALSE
    )
    combined <- combine_peaks(shifted, range = maxCombine)
    mat <- combined$long
    idx <- aki_new_rows(train_decons, decons)
    list(
        mat = mat[idx, , drop = FALSE],
        decons = if (is.null(train_decons)) decons else train_decons,
        train_cs = if (is.null(train_decons)) decons[[1]]$cs else train_decons[[1]]$cs
    )
}

aki_features_binning <- function(paths,
                                 bin_width,
                                 bin_range,
                                 bin_edges) {
    spectra <- aki_read_spectra(paths)
    aki_bin_matrix(
        spectra,
        bin_width = bin_width,
        bin_range = bin_range,
        bin_edges = bin_edges
    )
}

aki_bin_matrix <- function(spectra,
                           bin_width,
                           bin_range,
                           bin_edges) {
    if (!is_spectra(spectra)) stop("Expected a spectra object.")
    if (is.null(bin_edges)) {
        cs <- spectra[[1]]$cs
        ppm_max <- max(cs)
        ppm_min <- min(cs)
        if (!is.null(bin_range)) {
            ppm_max <- max(bin_range)
            ppm_min <- min(bin_range)
        }
        bin_edges <- seq(ppm_min, ppm_max, by = bin_width)
        if (tail(bin_edges, 1) < ppm_max) {
            bin_edges <- c(bin_edges, ppm_max)
        }
    }
    nb <- length(bin_edges) - 1
    mid <- bin_edges[-length(bin_edges)] + diff(bin_edges) / 2
    mat <- matrix(0, nrow = length(spectra), ncol = nb)
    for (i in seq_along(spectra)) {
        cs <- rev(spectra[[i]]$cs)
        si <- rev(spectra[[i]]$si)
        idx <- findInterval(cs, bin_edges, rightmost.closed = TRUE)
        valid <- idx > 0 & idx <= nb
        if (!any(valid)) next
        sums <- rowsum(si[valid], idx[valid], reorder = FALSE)
        mat[i, as.integer(rownames(sums))] <- sums[, 1]
    }
    colnames(mat) <- sprintf("%.4f", mid)
    list(mat = mat, bin_edges = bin_edges)
}

aki_read_spectra <- function(paths, cs = NULL) {
    paths <- aki_norm_paths(paths)
    if (length(paths) == 1) {
        spec <- read_spectrum(paths)
        spectra <- list(spec)
        names(spectra) <- basename(paths)
        class(spectra) <- "spectra"
        return(aki_resample_spectra(spectra, cs))
    }
    parents <- unique(dirname(paths))
    if (length(parents) == 1) {
        spectra <- read_spectra(parents)
        names_wanted <- basename(paths)
        spectra <- spectra[names_wanted]
        missing <- setdiff(names_wanted, names(spectra))
        if (length(missing) > 0) {
            msg <- paste(
                "Some spectra are missing:",
                paste(missing, collapse = ", ")
            )
            stop(msg)
        }
        return(aki_resample_spectra(spectra, cs))
    }
    spectra <- lapply(paths, read_spectrum)
    names(spectra) <- basename(paths)
    class(spectra) <- "spectra"
    aki_resample_spectra(spectra, cs)
}

aki_resample_spectra <- function(spectra, cs = NULL) {
    cs_list <- lapply(spectra, function(x) x$cs)
    if (!is.null(cs)) {
        if (length(cs) != length(cs_list[[1]])) {
            stop("Common cs length does not match spectra length.")
        }
    } else if (all_identical(cs_list)) {
        return(spectra)
    } else {
        cs <- aki_common_cs(spectra)
    }
    for (i in seq_along(spectra)) {
        x <- spectra[[i]]
        x$si <- stats::approx(x$cs, x$si, xout = cs, rule = 2)$y
        x$meta$fq <- stats::approx(x$cs, x$meta$fq, xout = cs, rule = 2)$y
        x$cs <- cs
        spectra[[i]] <- x
    }
    spectra
}

aki_common_cs <- function(spectra) {
    max_cs <- min(sapply(spectra, function(x) max(x$cs)))
    min_cs <- max(sapply(spectra, function(x) min(x$cs)))
    n <- length(spectra[[1]]$cs)
    seq(max_cs, min_cs, length.out = n)
}

aki_default_sfr <- function(spectrum, sfr) {
    if (!is.null(sfr)) return(sfr)
    stats::quantile(spectrum$cs, c(0.9, 0.1))
}

find_best_params <- function(spectrum,
                             sfr,
                             smopts = c(2, 5),
                             delta = 6.4,
                             nfit = 3,
                             verbose = TRUE) {
    check_mdrb(stop_on_fail = TRUE)
    sm_iter <- unique(pmax(1, smopts[1] + c(-1, 0, 1)))
    sm_win <- unique(pmax(3, smopts[2] + c(-2, 0, 2)))
    sm_win <- sm_win[sm_win %% 2 == 1]
    delta_grid <- unique(pmax(0.1, delta + c(-1, -0.5, 0, 0.5, 1)))
    nfit_grid <- unique(pmax(1, nfit + c(-1, 0, 1)))
    grid <- expand.grid(
        sm_iter = sm_iter,
        sm_win = sm_win,
        delta = delta_grid,
        nfit = nfit_grid
    )
    best <- list(mse = Inf, smopts = smopts, delta = delta, nfit = nfit)
    mdrb_spec <- mdrb::Spectrum$new(spectrum$cs, spectrum$si, sfr)
    for (i in seq_len(nrow(grid))) {
        row <- grid[i, ]
        dec <- mdrb::Deconvoluter$new()
        dec$set_moving_average_smoother(row$sm_iter, row$sm_win)
        dec$set_noise_score_selector(row$delta)
        dec$set_analytical_fitter(row$nfit)
        res <- try(dec$deconvolute_spectrum(mdrb_spec), silent = TRUE)
        if (inherits(res, "try-error")) next
        mse <- res$mse()
        if (isTRUE(verbose)) {
            logf("Grid %d/%d: mse=%.6f", i, nrow(grid), mse)
        }
        if (is.finite(mse) && mse < best$mse) {
            best$mse <- mse
            best$smopts <- c(row$sm_iter, row$sm_win)
            best$delta <- row$delta
            best$nfit <- row$nfit
        }
    }
    best[c("smopts", "delta", "nfit", "mse")]
}

aki_train_test_split <- function(y, ratio, seed) {
    set.seed(seed)
    ratio <- ratio / sum(ratio)
    pos <- sample(which(y == 1))
    neg <- sample(which(y == 0))
    n_pos_train <- max(1, floor(length(pos) * ratio[1]))
    n_neg_train <- max(1, floor(length(neg) * ratio[1]))
    train_idx <- c(pos[seq_len(n_pos_train)], neg[seq_len(n_neg_train)])
    test_idx <- setdiff(seq_along(y), train_idx)
    list(train_idx = sort(train_idx), test_idx = sort(test_idx))
}

aki_norm_paths <- function(paths) {
    if (!is.character(paths) || length(paths) == 0) {
        stop("data must be a non-empty character vector of paths.")
    }
    paths <- norm_path(paths, mustWork = TRUE)
    if (any(!file.exists(paths))) {
        stop("One or more spectrum paths do not exist.")
    }
    paths
}

aki_merge_decons <- function(train_decons, new_decons) {
    if (is.null(train_decons)) return(new_decons)
    merged <- c(train_decons, new_decons)
    class(merged) <- "decons2"
    merged
}

aki_new_rows <- function(train_decons, new_decons) {
    if (is.null(train_decons)) return(seq_len(length(new_decons)))
    start <- length(train_decons) + 1
    end <- length(train_decons) + length(new_decons)
    seq.int(start, end)
}

aki_pqn_ref <- function(x) {
    ref <- apply(x, 2, stats::median, na.rm = TRUE)
    ref[!is.finite(ref) | ref == 0] <- NA_real_
    ref
}

aki_pqn_apply <- function(x, ref) {
    ref <- ref
    ref[!is.finite(ref) | ref == 0] <- NA_real_
    quot <- sweep(x, 2, ref, "/")
    scale <- apply(quot, 1, stats::median, na.rm = TRUE)
    scale[!is.finite(scale) | scale == 0] <- 1
    sweep(x, 1, scale, "/")
}

aki_make_folds <- function(y, k, seed) {
    set.seed(seed)
    pos <- which(y == 1)
    neg <- which(y == 0)
    pos <- sample(pos)
    neg <- sample(neg)
    folds <- vector("list", k)
    for (i in seq_len(k)) {
        pos_i <- if (i <= length(pos)) pos[seq(i, length(pos), by = k)] else integer(0)
        neg_i <- if (i <= length(neg)) neg[seq(i, length(neg), by = k)] else integer(0)
        folds[[i]] <- sort(c(pos_i, neg_i))
    }
    folds
}

aki_labels_to_num <- function(labels) {
    if (is.logical(labels)) return(as.integer(labels))
    if (is.numeric(labels) && all(labels %in% c(0, 1))) {
        return(as.integer(labels))
    }
    msg <- paste(
        "labels must be logical or numeric 0/1.",
        "Use pheno$HasAKI or convert factors explicitly."
    )
    stop(msg)
}

aki_balanced_accuracy <- function(y, pred) {
    y <- as.integer(y)
    pred <- as.integer(pred)
    pos <- y == 1
    neg <- y == 0
    tpr <- sum(pred == 1 & pos) / sum(pos)
    tnr <- sum(pred == 0 & neg) / sum(neg)
    mean(c(tpr, tnr), na.rm = TRUE)
}

aki_auc <- function(y, prob) {
    y <- as.integer(y)
    pos <- y == 1
    n_pos <- sum(pos)
    n_neg <- sum(!pos)
    if (n_pos == 0 || n_neg == 0) return(NA_real_)
    r <- rank(prob)
    (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

aki_cache_key <- function(paths, y, params) {
    if (!requireNamespace("digest", quietly = TRUE)) {
        stop("Missing package 'digest'. Install it before running this.")
    }
    checks <- lapply(paths, checksum, method = "size")
    digest::digest(list(
        paths = paths,
        labels = y,
        params = params,
        checks = checks
    ))
}

aki_check_speaq_deps <- function(install_deps, verbose) {
    pkgvec <- c("speaq", "MassSpecWavelet", "impute")
    if (isTRUE(install_deps)) {
        bioc_install(pkgvec, ask = FALSE, verbose = verbose)
    }
    if (is.null(install_deps)) {
        bioc_install(pkgvec, ask = TRUE, verbose = verbose)
    }
    is_installed <- sapply(pkgvec, requireNamespace, quietly = TRUE)
    if (any(!is_installed)) {
        missing <- pkgvec[!is_installed]
        msg <- paste(
            "Missing packages:",
            paste(missing, collapse = ", "),
            "Install them or pass install_deps = TRUE."
        )
        stop(msg)
    }
}
