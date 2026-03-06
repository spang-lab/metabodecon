
# API #####

#' @noRd
#' @title Run AKI benchmark
#'
#' @description
#' Runs the full AKI benchmark pipeline from data loading to result summarization.
#'
#' @param seed Random seed for reproducibility.
#' @param kout Number of outer folds for cross-validation.
#' @param kin Number of inner folds for cross-validation.
#' @param fx Feature extraction method. Options are "bin", "deconvolute", "speaq".
#' @param mt Model type. Options are "lasso", "svm".
#' @param norm
#' Normalization method to apply on top of creatinine normalization, which is
#' always applied. Options are "none", "quantile", "pqn".
#'
#' @return
#' Invisibly returns a data frame with pooled predictions and true labels.
#'
run_aki_benchmark <- function(  seed = 1,
                                kout = 5,
                                kin = 5,
                                fx = "bin", # "deconvolute", "speaq"
                                mt = "lasso", # "svm"
                                norm = "none" # "quantile", "pqn"
                                ) {

    # Read in data
    xds_path <- metabodecon::download_example_datasets()
    aki_path <- file.path(xds_path, "bruker", "aki")
    meta <- read_aki_metadata(aki_path)
    spectra <- metabodecon::read_spectra(aki_path)
    stopifnot(all.equal(names(spectra), meta$sid))

    # Normalize spectra
    spectra_normed <- creatinine_normalize(spectra)

    # Prepare IDs for CV
    n <- length(spectra)
    folds <- seq_len(kout)
    test <- get_test_ids(kout, n, seed)
    train <- lapply(test, function(x) setdiff(seq_len(n), x))

    # Prepare y labels for CV
    y <- factor(meta$type, levels = c("Control", "AKI"))
    ytrain <- lapply(train, function(tr) y[tr])
    ytest <- lapply(test, function(te) y[te])

    # Prepare X matrices for CV
    if (fx == "bin") {
        X <- bin_spectra(spectra_normed)
    } else if (fx == "deconvolute") {
        decons <- deconvolute(spectra_normed, use_rust = TRUE)
        aligns <- align(decons)
        Xall <- t(get_si_mat(aligns))
        X <- Xall[, which(colSums(Xall) > 0), drop = FALSE]
    } else {
        stop("Invalid feature extraction method.")
    }
    Xtrain <- lapply(train, function(tr) X[tr, , drop = FALSE])
    Xtest <- lapply(test, function(te) X[te, , drop = FALSE])

    # Train models with hyperparameters selected from inner CV
    set.seed(seed)
    cvms <- map(train_cv_lasso, Xtrain, ytrain, kin)
    Y <- eval_cv_models(cvms, Xtest, ytest, folds)

    # Summarize results
    summarize_aki_benchmark_results(aki_path, meta, spectra, X, y, Y)
    invisible(Y)
}

# Helpers #####

map <- function(FUN, ..., MoreArgs = NULL) {
    FUN <- match.fun(FUN)
    dots <- list(...)
    .Internal(mapply(FUN, dots, MoreArgs))
}

read_aki_metadata <- function(aki_path) {
    meta <- read.csv(file.path(aki_path, "s_MTBLS24.txt"), sep = "\t")
    meta <- meta[, c("Sample.Name", "Factor.Value.Acute.Kidney.Injury.")]
    meta <- setNames(meta, c("sid", "type"))
    meta$type <- ifelse(meta$type == "Acute Kidney Injury", "AKI", "Control")
    # Fix samples with wrong date in their sample
    meta$sid[meta$sid == "AKI_8_24_105_110812"] <- "AKI_8_24_105_110816"
    meta$sid[meta$sid == "AKI_8_24_106_110812"] <- "AKI_8_24_106_110816"
    meta$sid[meta$sid == "AKI_8_24_107_110812"] <- "AKI_8_24_107_110816"
    meta$sid[meta$sid == "AKI_8_24_108_110812"] <- "AKI_8_24_108_110816"
    meta$sid[meta$sid == "AKI_8_24_109_110812"] <- "AKI_8_24_109_110816"
    meta$sid[meta$sid == "AKI_8_24_110_110812"] <- "AKI_8_24_110_110816"
    # Sort in the same order as the spectra on disk
    rownames(meta) <- meta$sid
    spectra_dirs <- grep("AKI_8_24", dir(aki_path), value = TRUE)
    meta <- meta[spectra_dirs, , drop = FALSE]
    meta
}

#' @noRd
#' @title Bin spectra into features
#'
#' @description Converts a list of spectra into a binned feature matrix.
#'
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param regions Two-column matrix with ppm regions to include.
#' @param binwidth Bin width in ppm.
#'
#' @details
#' The implementation uses one outer loop over spectra and computes all bin
#' sums within a spectrum by vectorized endpoint arithmetic.
#'
#' Nomenclature used in the code:
#'
#' - Region endpoints are named `reg_hi` and `reg_lo` and are initialized
#'   as `reg_hi = regions[, 1]` and `reg_lo = regions[, 2]`. If `regions[, 1]`
#'   is not greater than `regions[, 2]` for any row, an error is raised.
#'
#' - Bin endpoints are named `lo` and `hi` and built from the region endpoints.
#'
#' Algorithm (for each spectrum):
#'
#' 1. Sort points by chemical shift (`cs`) ascending and align intensities (`si`).
#' 2. Build bin intervals `(lo, hi]` from region bounds and `binwidth`.
#' 3. Compute cumulative sums `cum <- cumsum(si)`.
#' 4. For all bins at once, locate endpoints with `findInterval()`:
#'    - `i_hi`: number of points with `cs <= hi`
#'    - `i_lo`: number of points with `cs <= lo`
#' 5. Bin sums from cumulative differences:
#'    `bin_sum = S(i_hi) - S(i_lo)` with `S(k) = 0` for `k = 0` and
#'    `S(k) = cum[k]` for `k > 0`, which equals
#'    `sum(si[cs <= hi & cs > lo])`.
#'
#' @return Numeric matrix with samples in rows and bins in columns.
#'
#' @examples
#'
#' ## Toy Example
#' cs1 <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
#' si1 <- c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' si2 <- c( 2, 4, 1, 3, 5, 2, 6, 1, 4, 2)
#' ##                        ##
#' ##                  ##    ##
#' ##         ##       ##    ##    ##
#' ##         ##    ## ##    ##    ##
#' ##      ## ##    ## ## ## ##    ## ##
#' ##      ## ## ## ## ## ## ## ## ## ##
#' ##        |-----|-----|  |-----|
#' ##          4+1   3+5      6+1
#' spec1 <- structure(class = "spectrum", .Data = list(cs = cs1, si = si1))
#' spec2 <- structure(class = "spectrum", .Data = list(cs = cs1, si = si2))
#' spectra <- as_spectra(list(spec1, spec2))
#' regions <- rbind(c(9.5, 5.5), c(4.5, 2.5))
#' binwidth <- 2
#' bin_spectra(spectra, regions, binwidth)
#' ##            8.5 6.5 3.5
#' ## spectrum_1   3   7   13
#' ## spectrum_2   5   8   7
#'
#' ## Real Example
#' spectra <- sim
#' regions <- rbind(c(3.50, 3.45), c(3.35, 3.30))
#' binwidth <- 0.01
#' bin_spectra(sim, regions, binwidth = 0.01)
bin_spectra <- function(spectra,
                        regions = rbind(c(9.5, 6.5), c(4.5, 0.5)),
                        binwidth = 0.01) {

    # Region matrix convention: column 1 = region high, column 2 = region low.
    # Therefore every row must satisfy region high > region low.
    reg_hi <- regions[, 1]
    reg_lo <- regions[, 2]
    stopifnot(all(reg_hi > reg_lo))
    intwidths <- reg_hi - reg_lo
    nbins <- intwidths / binwidth
    stopifnot(all(is_int(nbins)))

    # Build global sorted edges from all regions.
    edges <- apply(regions, 1, function(x) seq(x[2], x[1], by = binwidth), simplify = FALSE)
    hi <- sort(unlist(lapply(edges, tail, -1)), decreasing = TRUE)
    lo <- hi - binwidth
    cen <- (lo + hi) / 2

    # Initialize output matrix
    nb <- length(cen)
    X <- matrix(0, nrow = length(spectra), ncol = nb)

    # Compute bin sums for each spectrum by vectorized endpoint arithmetic.
    for (i in seq_along(spectra)) {

        # Sort points by chemical shift in ascending order, as required by
        # findInterval().
        cs <- spectra[[i]]$cs
        si <- spectra[[i]]$si
        ord <- order(cs)
        cs <- cs[ord]
        si <- si[ord]

        # Locate bin endpoints and compute sums by cumulative differences.
        cum <- cumsum(si)
        idx_hi <- findInterval(hi, cs)
        idx_lo <- findInterval(lo, cs)
        sum_hi <- ifelse(idx_hi > 0L, cum[idx_hi], 0)
        sum_lo <- ifelse(idx_lo > 0L, cum[idx_lo], 0)
        X[i, ] <- sum_hi - sum_lo
    }

    colnames(X) <- cen
    rownames(X) <- get_names(spectra)
    attr(X, "bin_centers") <- cen
    X
}

creatinine_normalize <- function(spectra, cr = c(3.053, 3.011)) {
    stopifnot(length(cr) == 2, cr[1] > cr[2])
    spectra_normed <- spectra # copy spectra to new object for normalization
    for (i in seq_along(spectra)) {
        s <- spectra[[i]]
        idx <- which(s$cs >= cr[2] & s$cs <= cr[1])
        ci <- sum(s$si[idx]) # creatinine intensity
        spectra_normed[[i]]$si <- s$si / ci
    }
    spectra_normed
}

train_cv_lasso <- function(X, y, nfolds = 10) {
    glmnet::cv.glmnet(X, y, alpha = 1, nfolds = nfolds, family = "binomial")
}

eval_cv_model <- function(cv_model, Xte, yte, fold) {
    data.frame(
        sid = rownames(Xte),
        snr = sapply(strsplit(rownames(Xte), "_"), "[", 4),
        fold = fold,
        true = yte,
        prob = predict(cv_model, Xte, s = "lambda.min", type = "response")[, 1],
        pred = predict(cv_model, Xte, s = "lambda.min", type = "class")[, 1]
    )
}

eval_cv_models <- function(cvms, Xtest, ytest, folds) {
    Ys <- map(eval_cv_model, cvms, Xtest, ytest, folds)
    Y <- rbindlist(Ys)
    Y <- Y[order(Y$true, Y$prob), ]
    Y
}

summarize_aki_benchmark_results <- function(aki_path, meta, spectra, X, y, Y) {
    # metabodecon::tree_preview(aki_path)
    # print(headtail(meta, 3))
    # print(object.size(spectra), units = "Mb")
    # summary(head(spectra, 4))
    # plot_spectra_panel(spectra, meta)
    # draw_complex_heatmap(X, title = "Creatine Normed")
    # plot_probabilities(Y)
    summarize_performance(Y)
}

plot_probabilities <- function(Y) {
    plot_empty(
        xlim = c(1, nrow(Y)), ylim = c(0, 1),
        axes = TRUE, ylab = "Predicted probability",
        xlab = "Row Nr in Data.frame",
        main = "Lasso CV Predictions"
    )
    text(seq_along(Y$prob), Y$prob, labels = Y$snr, col = as.integer(Y$true))
}

summarize_performance <- function(Y) {
    accuracy <- mean(Y$true == Y$pred)
    auc <- pROC::auc(Y$true, Y$prob)
    fmt <- "Accuracy: %.2f%%\nAUC: %.3f\n"
    table(Y$true, Y$pred)
    cat(sprintf(fmt, accuracy * 100, auc))
}
