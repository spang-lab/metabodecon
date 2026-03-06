
# API #####

run_aki_benchmark_rest <- function() {
    tscores <- compute_t_scores(X_qn, y)
    pvals <- 2 * pt(-abs(tscores), df = nrow(X_qn) - 2)
    withr::with_par(
        new = list(mfrow = c(1, 2)),
        code = {
            hist(tscores, breaks = 100, main = "T-scores for AKI vs Control", xlab = "T-score")
            hist(pvals, breaks = 100, main = "P-values for AKI vs Control", xlab = "P-value")
        }
    )
    df <- data.frame(feature = colnames(X_qn), tscore = tscores, pval = pvals)
    df <- df[order(df$pval), ]
    head(df, 20)
    # top2_fixed <- df$feature[1:2]

    set.seed(1)
    te <- get_test_ids(seq_len(nrow(X_qn)), nfolds = 5, seed = 1, y = y)
    Y <- estimate_svm_performance(
        X = X_qn,
        y = y,
        test_ids = te,
        costs = 2^(0:5),
        gammas = 2^(-8:-3),
        nfeats = 2:3,
        inner_nfolds = 5,
        ncores = 5,
        verbose = FALSE
    )
    Yord <- Y[order(Y$true, Y$prob), ]
    plot_empty(
        xlim = c(1, nrow(Yord)), ylim = c(0, 1),
        axes = TRUE, ylab = "Predicted score",
        xlab = "Row Nr in Data.frame",
        main = "Lasso CV Predictions"
    )
    text(
        x = seq_along(Yord$prob),
        y = Yord$prob,
        labels = as.character(Yord$sample),
        col = as.integer(Yord$true == "AKI") + 1
    )

    # Show accuracy and AUC
    accuracy <- mean(Y$true == Y$pred)
    auc <- pROC::auc(Y$true, Y$prob)
    fmt <- "Accuracy: %.2f%%\nAUC: %.3f\n"
    table(Y$true, Y$pred)
    cat(sprintf(fmt, accuracy * 100, auc))
}

# Plotting #####

#' @noRd
#' @title Plot selected spectra panel
#' @description Draws selected spectra in a multi-panel layout.
#' @param spectra Spectra object.
#' @param ids Integer indices to plot.
#' @param meta Optional data frame with `sample` and `type` columns.
#' @param foc_rgn Numeric length-2 vector for focus region in ppm.
#' @param mfrow Integer length-2 panel layout.
#' @return Invisibly returns plotted indices.
#' @examples
#' # See AKI vignette for a full example.
plot_spectra_panel <- function(spectra,
                               meta = NULL,
                               ids = c(30, 90),
                               foc_rgn = c(4.5, 1.5),
                               mfrow = NULL) {
    spectra <- spectra[ids]

    if (is.null(mfrow)) {
        n <- length(spectra)
        s <- sqrt(n)
        lo <- floor(s)
        hi <- ceiling(s)

        if (lo * lo >= n) {
            mfrow <- c(lo, lo)
        } else if (lo * hi >= n) {
            mfrow <- c(lo, hi)
        } else if (hi * hi >= n) {
            mfrow <- c(hi, hi)
        } else {
            stop("could not derive mfrow from n; possible rounding issue")
        }
    }

    withr::with_par(list(
        mfrow = mfrow,
        font.main = 1,
        cex.main = par("cex")
    ), {
        for (j in seq_along(spectra)) {
            i <- ids[j]
            main <- NULL
            if (!is.null(meta)) {
                main <- paste0(meta$sample[i], " (", meta$type[i], ")")
            }
            plot_spectrum(
                spectra[[j]],
                foc_rgn = foc_rgn,
                sub1 = list(main = main, bt_axis = TRUE)
            )
        }
    })

    invisible(ids)
}

#' @noRd
#' @title Draw base heatmap
#' @description Draws a base heatmap with robust colors.
#' @param M Numeric matrix with samples in rows and bins in columns.
#' @param y Unused. Kept for API compatibility.
#' @param title Optional plot title. Set `NULL` to hide title.
#' @param x_axis X-axis title.
#' @param mar Plot margins in the order bottom, left, top, right.
#' @param mai Plot margins in inches as bottom, left, top, right.
draw_complex_heatmap <- function(M,
                                 y = NULL,
                                 title = NULL,
                                 x_axis = "ppm",
                                 mar = c(4, 1, 1, 1),
                                 mai = NULL) {
    M <- as.matrix(M)
    nr <- nrow(M)
    nc <- ncol(M)
    if (nr == 0 || nc == 0) {
        plot.new()
        return(invisible(NULL))
    }

    vals <- as.numeric(M)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) {
        plot.new()
        return(invisible(NULL))
    }
    q <- as.numeric(stats::quantile(
        vals,
        probs = c(0.01, 0.10, 0.50, 0.75, 0.90, 0.99),
        na.rm = TRUE,
        type = 8
    ))
    lo <- q[1]
    hi <- q[length(q)]
    M_plot <- pmin(pmax(M, lo), hi)

    br <- as.numeric(stats::quantile(
        vals,
        probs = seq(0, 1, length.out = 257),
        na.rm = TRUE,
        type = 8
    ))
    if (length(unique(br)) < 2) {
        br <- seq(lo - 0.5, hi + 0.5, length.out = 257)
    }

    z <- matrix(as.integer(cut(
        M_plot,
        breaks = br,
        include.lowest = TRUE,
        labels = FALSE
    )), nrow = nr, ncol = nc)
    z[!is.finite(z)] <- 1L

    ppm <- attr(M, "bin_centers")
    if (is.null(ppm) || length(ppm) != nc) {
        ppm <- suppressWarnings(as.numeric(colnames(M)))
    }
    if (length(ppm) != nc || anyNA(ppm)) {
        ppm <- seq_len(nc)
    }
    n_ticks <- min(8L, nc)
    x_at <- unique(round(seq(1, nc, length.out = n_ticks)))
    x_labs <- format(ppm[x_at], digits = 3)

    pal <- grDevices::colorRampPalette(c(
        "#08306b", "#2171b5", "#6baed6",
        "#fee090", "#fdae61", "#d73027"
    ))(256)

    popts <- list(mar = mar)
    if (!is.null(mai)) {
        popts$mai <- mai
    }

    withr::with_par(popts, {
        image(
            x = seq_len(nc),
            y = seq_len(nr),
            z = t(z[rev(seq_len(nr)), , drop = FALSE]),
            col = pal,
            zlim = c(1, length(pal)),
            xaxs = "i",
            yaxs = "i",
            axes = FALSE,
            xlab = "",
            ylab = "",
            main = title
        )
        axis(1, at = x_at, labels = x_labs, las = 1, cex.axis = 0.8)
        title(xlab = x_axis)
        box()
    })

    invisible(NULL)
}

#' @noRd
#' @title Signed log2 transform
#' @description
#' Applies a signed log2 transform defined as sign(x) * log2(1 + abs(x)).
#' @param x Numeric vector or matrix.
symlog2 <- function(x) {
    sign(x) * log2(1 + abs(x))
}

#' @noRd
#' @title Inverse signed log2 transform
#' @description
#' Inverse of `symlog2()`, defined as sign(x) * (2^abs(x) - 1).
#' @param x Numeric vector or matrix on symlog2 scale.
inv_symlog2 <- function(x) {
    sign(x) * (2^abs(x) - 1)
}

#' @noRd
#' @title Plot sample-wise boxplots for full intensities
#' @description
#' Draws one boxplot per sample using the full intensity range.
#' @param M Numeric matrix with samples in rows and bins in columns.
#' @param y Optional class labels with values "Control" and "AKI".
#' @param tick_every Show every k-th sample index on y-axis.
#' @param y_axis X-axis label for the boxplot values.
#' @param x_transform Transform on x-values: "symlog2" or "none".
#' @param mar Plot margins in the order bottom, left, top, right.
#' @param mai Plot margins in inches as bottom, left, top, right.
plot_sample_boxplots <- function(M,
                                 y = NULL,
                                 tick_every = 20,
                                 y_axis = "Intensity",
                                 x_transform = c("symlog2", "none"),
                                 mar = c(5, 4, 0, 0),
                                 mai = NULL) {
    M <- as.matrix(M)
    x_transform <- match.arg(x_transform)

    if (x_transform == "symlog2") M <- symlog2(M)
    n <- nrow(M)
    at <- seq(1, n, by = tick_every)
    if (tail(at, 1) != n) at <- c(at, n)
    x_full <- lapply(seq_len(n), function(i) {
        x <- as.numeric(M[i, ])
        x[is.finite(x)]
    })
    empty <- vapply(x_full, length, integer(1)) == 0L
    if (any(empty)) {
        x_full[empty] <- list(NA_real_)
    }
    cols <- rep("#bdbdbd", n)
    if (!is.null(y) && length(y) == n) {
        cols[y == "Control"] <- "#33a02c"
        cols[y == "AKI"] <- "#ff7f00"
    }
    popts <- list(mar = mar)
    if (!is.null(mai)) {
        popts$mai <- mai
    }

    withr::with_par(popts, {
        pos <- rev(seq_len(n))
        x_all <- unlist(x_full, use.names = FALSE)
        x_all <- x_all[is.finite(x_all)]
        if (length(x_all) == 0) {
            x_all <- 0
        }

        q <- stats::quantile(x_all, probs = c(0.01, 0.99), na.rm = TRUE,
            names = FALSE)
        x_rng <- as.numeric(q)
        if (!all(is.finite(x_rng)) || x_rng[1] >= x_rng[2]) {
            x_rng <- range(x_all, finite = TRUE)
        }

        dx <- diff(x_rng)
        if (is.finite(dx) && dx > 0) {
            x_rng <- x_rng + c(-0.05, 0.05) * dx
        } else {
            x_rng <- x_rng + c(-0.5, 0.5)
        }

        x_ticks <- pretty(x_rng, n = 6)

        if (x_transform == "symlog2") {
            x_labs <- format(inv_symlog2(x_ticks), digits = 3, trim = TRUE)
        } else {
            x_labs <- format(x_ticks, digits = 3, trim = TRUE)
        }

        boxplot(
            x_full,
            at = pos,
            horizontal = TRUE,
            outline = TRUE,
            yaxt = "n",
            xaxt = "n",
            col = cols,
            border = "#4d4d4d",
            xlab = y_axis,
            ylab = "",
            main = "",
            xlim = x_rng,
            ylim = c(0.5, n + 0.5),
            yaxs = "i"
        )
        axis(1, at = x_ticks, labels = x_labs, las = 1, cex.axis = 0.8)
        axis(2, at = n - at + 1, labels = at, las = 1, cex.axis = 0.7)
    })

    invisible(NULL)
}

#' @noRd
#' @title Visualize feature matrix in a 1x2 panel
#' @description
#' Draws sample-wise boxplots and a heatmap of raw feature intensities.
#' @param X Numeric matrix with samples in rows and bins in columns.
#' @param y Optional class labels with values "Control" and "AKI".
#' @return Invisible list with `raw` matrix.
visualize_feature_matrix <- function(X,
                                     y = NULL) {
    X <- as.matrix(X)
    withr::local_par(list(mfrow = c(1, 2)))
    plot_sample_boxplots(
        X,
        y = y,
        tick_every = 10,
        y_axis = "Intensity (symlog2)",
        mar = c(5, 4, 0, 0)
    )
    draw_complex_heatmap(
        X,
        title = "raw",
        x_axis = "ppm",
        mar = c(5, 4, 0, 0)
    )

    invisible(list(raw = X))
}

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
#' @param sample_ids Optional sample IDs used for x-axis annotations.
#' @param main Plot title.
#' @return `NULL` invisibly.
#' @examples
#' tmp <- tempfile(fileext = ".pdf")
#' grDevices::pdf(tmp)
#' plot_prob_scatter(c(0.1, 0.7, 0.4), c(0, 1, 0), "Example")
#' grDevices::dev.off()
#' @noRd
plot_prob_scatter <- function(prob,
                              y,
                              main = "Predicted probabilities",
                              sample_ids = NULL) {
    y01 <- as_binary01(y)
    o <- order(y01, prob)
    prob <- prob[o]
    y01 <- y01[o]

    if (!is.null(sample_ids)) {
        sample_ids <- sample_ids[o]
    }

    cols <- ifelse(y01 == 1, "#d95f02", "#1b9e77")

    x <- seq_along(prob)
    x_labs <- x
    if (!is.null(sample_ids)) {
        sid <- as.character(sample_ids)
        parts <- strsplit(sid, "_", fixed = TRUE)
        num <- vapply(parts, function(z) {
            if (length(z) >= 4) z[[4]] else NA_character_
        }, character(1))
        bad <- is.na(num) | !grepl("^[0-9]+$", num)
        num[bad] <- sid[bad]
        x_labs <- num
    }

    withr::with_par(list(mar = c(7, 4, 2, 1)), {
        plot(
            x = x,
            y = prob,
            col = cols,
            pch = 19,
            xaxt = "n",
            xlab = if (is.null(sample_ids)) "sample" else "sample number",
            ylab = "P(AKI)",
            main = main,
            ylim = c(0, 1)
        )
        axis(1, at = x, labels = x_labs, las = 2, cex.axis = 0.6)
        abline(h = 0.5, lty = 2)
        legend(
            "topright",
            pch = 19,
            col = c("#1b9e77", "#d95f02"),
            legend = c("Control", "AKI"),
            bg = "white",
            cex = 0.85
        )
    })
}

# Modeling #####

compute_t_scores <- function(X, y)
{
    y <- as_binary01(y)
    i1 <- which(y == 1)
    i0 <- which(y == 0)
    m1 <- colMeans(X[i1, , drop = FALSE])
    m0 <- colMeans(X[i0, , drop = FALSE])
    v1 <- apply(X[i1, , drop = FALSE], 2, var)
    v0 <- apply(X[i0, , drop = FALSE], 2, var)
    s <- sqrt(v1 / length(i1) + v0 / length(i0))
    t <- abs((m1 - m0) / s)
    t[!is.finite(t)] <- 0
    t
}

#' @title Estimate SVM Performance
#' @description Evaluates radial SVM using nested cross-validation.
#' @param X Feature matrix with samples in rows.
#' @param y Binary labels (factor or 0/1).
#' @param test_ids List of test index vectors for outer folds.
#' @param costs Cost parameter candidates for grid search.
#' @param gammas Gamma parameter candidates for grid search.
#' @param nfeats Feature count candidates for grid search.
#' @param inner_nfolds Number of inner CV folds for hyperparameter tuning.
#' @param verbose Print progress messages.
#' @param ncores Number of CPU cores for parallel grid search.
#' @param sample_ids Optional sample IDs. Defaults to row indices.
#' @return Data frame with columns `sample`, `fold`, `cost`, `gamma`,
#'   `nfeat`, `true`, `prob`, `pred`.
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(200), 20)
#' y <- factor(rep(c("Control", "AKI"), each = 10))
#' te <- list(1:5, 6:10, 11:15, 16:20)
#' Y <- estimate_svm_performance(X, y, te,
#'   costs = 2^(0:2), gammas = 2^(-5:-3), nfeats = 2:3)
#' head(Y)
#' }
#' @export
estimate_svm_performance <- function(X, y, test_ids, costs, gammas, nfeats,
                                     inner_nfolds = 20, verbose = TRUE, ncores = 8L,
                                     sample_ids = NULL) {
    if (is.null(sample_ids)) sample_ids <- seq_len(nrow(X))
    if (is.null(rownames(X))) rownames(X) <- sample_ids

    eval_fold <- function(te, fold_idx) {
        if (verbose) {
            ts <- format(Sys.time(), "%H:%M:%S")
            cat(sprintf("[%s] Processing fold %d/%d (n_test=%d)\n", ts, fold_idx, length(test_ids), length(te)))
        }
        tr <- setdiff(seq_len(nrow(X)), te)
        p <- find_best_svm_params(
            X = X[tr, ],
            y = y[tr],
            costs = costs,
            gammas = gammas,
            nfeats = nfeats,
            normalize = TRUE,
            k = inner_nfolds,
            verbose = verbose,
            fold_prefix = sprintf("[Fold %d] ", fold_idx),
            ncores = ncores
        )
        svm <- train_svm(X[tr, ], y[tr], p$cost, p$gamma, p$nfeat, TRUE)
        yhat <- predict_svm(svm, X[te, ])
        cls <- factor(yhat$cls, levels = c(0, 1), labels = levels(y))
        fold_auc <- auc(y[te], yhat$prob)
        fold_acc <- mean(y[te] == cls)

        if (verbose) {
            ts <- format(Sys.time(), "%H:%M:%S")
            cat(sprintf("[%s] Completed: AUC=%.3f, ACC=%.3f\n\n", ts, fold_auc, fold_acc))
        }

        list(
            ids = te,
            fold = fold_idx,
            cost = p$cost,
            gamma = p$gamma,
            nfeat = p$nfeat,
            true = as.character(y[te]),
            prob = yhat$prob,
            pred = as.character(cls),
            auc = fold_auc,
            acc = fold_acc
        )
    }

    if (verbose) {
        ts <- format(Sys.time(), "%H:%M:%S")
        cat(sprintf("\n[%s] Starting nested CV: %d outer, %d inner folds\n", ts, length(test_ids), inner_nfolds))
        cat(sprintf("[%s] Grid: %d costs * %d gammas * %d nfeats = %d combinations\n",
            ts, length(costs), length(gammas), length(nfeats), length(costs) * length(gammas) * length(nfeats)))
        cat(sprintf("[%s] Using %d cores\n\n", ts, ncores))
    }

    res <- lapply(seq_along(test_ids), function(i) eval_fold(test_ids[[i]], i))

    Y <- do.call(rbind, lapply(res, function(r) {
        data.frame(
            sample = sample_ids[r$ids],
            fold = r$fold,
            cost = r$cost,
            gamma = r$gamma,
            nfeat = r$nfeat,
            true = r$true,
            prob = r$prob,
            pred = r$pred,
            stringsAsFactors = FALSE
        )
    }))
    rownames(Y) <- NULL

    # Pooled metrics (combine all fold predictions)
    prob_all <- rep(NA_real_, nrow(X))
    cls_all <- rep(NA_integer_, nrow(X))
    for (r in res) {
        prob_all[r$ids] <- r$prob
        cls_all[r$ids] <- as_binary01(factor(r$pred, levels = levels(y)))
    }

    auc_native <- auc(y, prob_all)
    acc_native <- mean(cls_all == as_binary01(y))
    acc_thresh <- mean((prob_all >= 0.5) == as_binary01(y))
    mean_fold_auc <- mean(vapply(res, function(r) r$auc, numeric(1)),
        na.rm = TRUE)
    mean_fold_acc <- mean(vapply(res, function(r) r$acc, numeric(1)),
        na.rm = TRUE)

    if (verbose) {
        ts <- format(Sys.time(), "%H:%M:%S")
        cat(sprintf("[%s] Pooled AUC: %.3f (mean fold: %.3f)\n",
            ts, auc_native, mean_fold_auc))
        cat(sprintf("[%s] Pooled ACC (native): %.3f (mean fold: %.3f)\n",
            ts, acc_native, mean_fold_acc))
        cat(sprintf("[%s] Pooled ACC (thresh): %.3f\n\n", ts,
            acc_thresh))
    }

    Y
}

#' @title Train SVM
#' @description Trains radial SVM with optional quantile normalization and feature selection.
#' @param X Feature matrix.
#' @param y Binary labels (0/1 or factor).
#' @param cost SVM cost parameter.
#' @param gamma SVM gamma parameter.
#' @param nfeat Number of top features to select.
#' @param normalize Apply quantile normalization.
#' @return List with fit, ref, idx, normalize flag, y_train.
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(80), 20)
#' y <- rep(0:1, each = 10)
#' m <- train_svm(X, y, cost = 1, gamma = 0.01, nfeat = 3)
#' }
#' @export
train_svm <- function(X, y, cost, gamma, nfeat, normalize = TRUE) {
    stopifnot(requireNamespace("e1071", quietly = TRUE))
    y <- as_binary01(y)
    stopifnot(length(unique(y)) == 2)
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
        gamma = gamma
    )
    list(fit = fit, ref = ref, idx = idx, normalize = normalize, y_train = y)
}

#' @title Predict SVM
#' @description Predicts probabilities and classes from trained SVM.
#' @param model Model from train_svm().
#' @param X Feature matrix.
#' @return List with prob (probabilities) and cls (predicted classes 0/1).
#' @examples
#' \dontrun{
#' m <- train_svm(X, y, 1, 0.01, 3)
#' p <- predict_svm(m, X_test)
#' }
#' @export
predict_svm <- function(model, X) {
    stopifnot(requireNamespace("e1071", quietly = TRUE))
    Xn <- X
    if (model$normalize) Xn <- quantile_normalize(Xn, model$ref)
    pred <- stats::predict(model$fit, Xn[, model$idx, drop = FALSE], decision.values = TRUE)
    score <- as.numeric(attr(pred, "decision.values"))
    cls <- as.integer(as.character(pred))

    # Ensure probabilities align with class predictions
    # If most class-1 predictions have negative scores, flip direction
    if (sum(cls == 1) > 0 && sum(cls == 0) > 0) {
        m1 <- mean(score[cls == 1], na.rm = TRUE)
        m0 <- mean(score[cls == 0], na.rm = TRUE)
        if (is.finite(m1) && is.finite(m0) && m1 < m0) score <- -score
    }

    score[!is.finite(score)] <- 0
    p1 <- 1 / (1 + exp(-score))
    list(prob = p1, cls = cls)
}

#' @title Score Hyperparameter Grid
#' @description Evaluates all hyperparameter combinations on train/test split.
#' @param X Feature matrix.
#' @param y Binary labels.
#' @param tr Training indices.
#' @param te Test indices.
#' @param costs Cost values to try.
#' @param gammas Gamma values to try.
#' @param nfeats Feature counts to try.
#' @param normalize Apply quantile normalization.
#' @param ncores  Number of CPU cores for parallel evaluation.
#' @return Data frame with cost, gamma, nfeat, auc, acc columns.
#' @export
score_grid <- function(X, y, tr, te, costs, gammas, nfeats, normalize = TRUE, ncores = 8L) {
    y <- as_binary01(y)
    G <- expand.grid(cost = costs, gamma = gammas, nfeat = nfeats)
    G$auc <- NA_real_
    G$acc <- NA_real_
    if (length(unique(y[tr])) < 2 || length(unique(y[te])) < 2) return(G)
    eval_one <- function(i) {
        f <- train_svm(X[tr, , drop = FALSE], y[tr], G$cost[i], G$gamma[i], G$nfeat[i], normalize)
        p <- predict_svm(f, X[te, , drop = FALSE])
        list(auc = auc(y[te], p$prob), acc = mean(y[te] == p$cls))
    }
    results <- parallel::mclapply(seq_len(nrow(G)), eval_one, mc.cores = ncores, mc.preschedule = TRUE)
    for (i in seq_len(nrow(G))) {
        G$auc[i] <- results[[i]]$auc
        G$acc[i] <- results[[i]]$acc
    }
    G
}

#' @title Find Best Hyperparameters
#' @description Uses k-fold CV to select optimal cost, gamma, nfeat by max mean AUC.
#' @param X Feature matrix.
#' @param y Binary labels.
#' @param costs Cost candidates.
#' @param gammas Gamma candidates.
#' @param nfeats Feature count candidates.
#' @param normalize Apply quantile normalization.
#' @param k Number of inner CV folds.
#' @param verbose Print progress.
#' @param fold_prefix Prefix for progress messages.
#' @param ncores Number of cores.
#' @return List with cost, gamma, nfeat, grid.
#' @export
find_best_svm_params <- function(X, y, costs = 2^(0:5), gammas = 2^(-8:-3), nfeats = 2:3,
                                 normalize = TRUE, k = 20, verbose = FALSE, fold_prefix = "", ncores = 8L) {
    te_ids <- get_test_ids(seq_len(nrow(X)), nfolds = k, y = y)
    G <- expand.grid(cost = costs, gamma = gammas, nfeat = nfeats)
    A <- matrix(NA_real_, nrow = nrow(G), ncol = length(te_ids))
    for (i in seq_along(te_ids)) {
        if (verbose) {
            ts <- format(Sys.time(), "%H:%M:%S")
            cat(sprintf(
                "[%s] %sInner CV: fold %d/%d (using %d cores)\n",
                ts, fold_prefix, i, length(te_ids), ncores
            ))
        }
        te <- te_ids[[i]]
        tr <- setdiff(seq_len(nrow(X)), te)
        Gi <- score_grid(X, y, tr, te, costs, gammas, nfeats, normalize, ncores)
        A[, i] <- Gi$auc
    }
    G$auc <- rowMeans(A, na.rm = TRUE)
    G$auc[!is.finite(G$auc)] <- -Inf
    b <- G[which.max(G$auc), ]
    list(cost = b$cost, gamma = b$gamma, nfeat = b$nfeat, grid = G)
}

# Helpers #####

#' @noRd
#' @title Coerce binary labels to 0/1
#' @description
#' Validates binary labels and maps them to integer 0/1 codes.
#' If `y` is a factor, arbitrary two levels are supported. The first level maps
#' to 0 and the second level maps to 1.
#' @param y Label vector.
#' @return Integer vector with values 0/1.
#' @examples
#' as_binary01(factor(c("Control", "AKI", "Control")))
as_binary01 <- function(y) {
    arg_name <- "y"
    if (is.factor(y)) {
        if (nlevels(y) != 2) {
            stop(sprintf("%s must have exactly 2 levels.", arg_name))
        }
        out <- as.integer(y) - 1L
        if (anyNA(out)) {
            stop(sprintf("%s must not contain NA.", arg_name))
        }
        return(out)
    }

    vals <- y[!is.na(y)]
    if (length(vals) == 0) {
        stop(sprintf("%s must contain non-missing values.", arg_name))
    }

    levs <- sort(unique(vals))
    if (length(levs) != 2) {
        stop(sprintf("%s must have exactly 2 classes.", arg_name))
    }

    if (is.logical(y)) {
        out <- as.integer(y)
        if (anyNA(out)) {
            stop(sprintf("%s must not contain NA.", arg_name))
        }
        return(out)
    }

    if ((is.integer(y) || is.numeric(y)) && all(levs == c(0, 1))) {
        out <- as.integer(y)
        if (anyNA(out)) {
            stop(sprintf("%s must not contain NA.", arg_name))
        }
        return(out)
    }

    fac <- factor(y)
    if (nlevels(fac) != 2) {
        stop(sprintf("%s must have exactly 2 classes.", arg_name))
    }
    out <- as.integer(fac) - 1L
    if (anyNA(out)) {
        stop(sprintf("%s must not contain NA.", arg_name))
    }
    out
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
    y <- as_binary01(y)
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
    y <- as_binary01(y)
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

    y <- as_binary01(y)
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


# Old #####

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
        best <- aki_find_best_params(
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
