
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
#' @param nfit Number of peaks to fit per spectrum in deconvolution.
#' @param smit Smoothing iterations for deconvolution.
#' @param smws Smoothing window size for deconvolution.
#' @param delta Delta parameter for deconvolution.
#' @param npmax
#' Maximum number of peaks to allow in deconvolution.
#' Seeting npmax to any value greater 0 causes the deconvolution
#' to ignore the given nfit, smit, smws and delta parameters and
#' instead performs a grid search over a predefined grid of deconvolution
#' parameters and selects the parameter combination we lowest MAE and
#' less than npmax peaks.
#' @param maxShift Maximum shift allowed in alignment.
#' @param maxCombine Maximum number of peaks to combine in alignment.
#' @param verbose Whether to print verbose output during deconvolution and alignment.
#'
#' @return
#' Invisibly returns a data frame with pooled predictions and true labels.
#'
run_aki_benchmark <- function(  aki = read_aki_data(),
                                seed = 1,
                                kout = 3,
                                kin = 3,
                                fx = "deconvolute",
                                mt = "lasso",
                                norm = "none",
                                nfit = 3,
                                smit = 2,
                                smws = 5,
                                delta = 6.4,
                                npmax = 0,
                                maxShift = 50,
                                maxCombine = 5,
                                verbose = TRUE,
                                nworkers = 1,
                                use_rust = TRUE,
                                sfr = NULL,
                                wshw = 0) {

    if (isFALSE(verbose)) local_options(toscutil.logf.file = nullfile())

    logf("Preparing IDs for CV")
    spectra <- aki$spectra
    n <- length(spectra)
    folds <- seq_len(kout)
    test <- get_test_ids(kout, n, seed)
    train <- lapply(test, function(x) setdiff(seq_len(n), x))

    logf("Preparing labels for CV")
    meta <- aki$meta
    y <- factor(meta$type, levels = c("Control", "AKI"))
    ytrain <- lapply(train, function(tr) y[tr])
    ytest <- lapply(test, function(te) y[te])

    logf("Preparing matrices for CV")
    if (fx == "bin") {
        X <- bin_spectra(spectra)
    } else if (fx == "deconvolute") {
        decons <- deconvolute(
            spectra, nfit=nfit, smopts=c(smit, smws), delta=delta,
            sfr=sfr, wshw=wshw, ask=FALSE, force=FALSE, verbose=verbose,
            nworkers=nworkers, use_rust=use_rust, npmax=npmax
        )
        aligns <- align(
            decons, maxShift=maxShift, verbose=verbose,
            nworkers=nworkers
        )
        Xall <- t(get_si_mat(aligns))
        X <- Xall[, which(colSums(Xall) > 0), drop=FALSE]
    } else {
        stop("Invalid feature extraction method.")
    }
    Xtrain <- lapply(train, function(tr) X[tr, , drop=FALSE])
    Xtest <- lapply(test, function(te) X[te, , drop=FALSE])

    logf("Training models with hyperparameters selected from inner CV")
    set.seed(seed)
    cvms <- map(train_cv_lasso, Xtrain, ytrain, kin)
    Y <- eval_cv_models(cvms, Xtest, ytest, folds)

    logf("Summarizing results")
    acc <- mean(Y$true == Y$pred)
    auc <- auc(Y$true, Y$prob)
    logf("Accuracy: %.2f%%, AUC: %.3f", acc * 100, auc)

    logf("Finished AKI benchmark")
    invisible(list(acc = acc, auc = auc))
}

try_aki_benchmark <- function(  aki = read_aki_data(),
                                seed = 1,
                                kout = 3,
                                kin = 3,
                                fx = "deconvolute",
                                mt = "lasso",
                                norm = "none",
                                nfit = 3,
                                smit = 2,
                                smws = 5,
                                delta = 6.4,
                                npmax = 0,
                                maxShift = 50,
                                maxCombine = 5,
                                verbose = TRUE,
                                nworkers = 1,
                                use_rust = TRUE,
                                sfr = NULL,
                                wshw = 0) {
     sfrstr <- if (is.null(sfr)) "NULL" else paste(sfr, collapse = "-")
     fmt <- paste(
        "Starting: seed=%d, nfit=%d, smit=%d, smws=%d, delta=%.1f, npmax=%d,",
        "maxShift=%d, maxCombine=%d, sfr=%s, wshw=%d"
    )
    logf(fmt, seed, nfit, smit, smws, delta, npmax, maxShift, maxCombine, sfrstr, wshw)
    x <- tryCatch(
        expr = {
            run_aki_benchmark(
                aki, seed, kout, kin, fx, mt, norm, nfit, smit, smws, delta,
                npmax, maxShift, maxCombine, verbose, nworkers, use_rust,
                sfr, wshw
            )
        },
        error = function(e) {
            message("Error in run_aki_benchmark: ", conditionMessage(e))
            invisible(list(acc = NA, auc = NA))
        }
    )
    fmt <- paste(
        "Done: seed=%d, nfit=%d, smit=%d, smws=%d, delta=%.1f, npmax=%d,",
        "maxShift=%d, maxCombine=%d, sfr=%s, wshw=%s. Acc=%.2f%%, AUC=%.3f."
    )
    logf(
        fmt, seed, nfit, smit, smws, delta, npmax, maxShift, maxCombine,
        sfrstr, wshw, x$acc * 100, x$auc
    )
    x
}

#' @noRd
#' @title Run AKI benchmark over parameter grid
#'
#' @description
#' Iterates over a grid of deconvolution and alignment parameters and calls
#' `run_aki_benchmark()` for each parameter combination. This is statistically
#' invalid because of two reasons:
#'
#' 1. We are doing the preprocessing (deconvolution and alignment) on the full
#'    dataset instead of only on the training folds.
#' 2. We are evaluating thousands of models on the same test set. I.e., even if
#'    the different preprocssing parameters have no effect on the performance,
#'    just by chance some parameter combinations will perform better than others
#'    on the test folds, due to random assignment of samples to folds. I.e., it
#'    will look like a certain preprocessing parameter combination is better than
#'    the rest, even if in reality there is no difference.
#'
#' The correct way to do this, would be to do the deconvolution and alignment as
#' part of the inner CV loop, and then only apply the best parameter combination
#' to the outer test fold. This will be implemented in [metabodecon::cv_mdm_grid()].
#'
#' However, because this will multiply runtime by kout*kin and also requires
#' additional implementation effort, we will first run this invalid version to
#' get a rough idea of whether the preprocessing parameters have any effect at
#' all. We can counter the "multiple evaluations problem" by repeat the
#' deconvolution with default parameter multiple times with different random
#' seeds to estimate the variance of the performance estimates and then check
#' whether some preprocessings lead to a performance increase that is far
#' outside the expected range.
#'
run_aki_benchmarks_invalid <- function() {

    # Prepare input data
    aki <- read_aki_data()

    # Run one benchmark with npmax = 1000. This will trigger a grid search to
    # find the best deconvolution parameters producing less than 1000 peaks.
    # Grid search results are cached automatically.
    perfMC <- run_aki_benchmark(
        aki=aki, seed=1, kout=10, kin=10, fx="deconvolute", mt="lasso",
        norm="none", nfit=3, smit=2, smws=5, delta=6.4, npmax=1000,
        maxShift=50, maxCombine=16, verbose=TRUE,
        nworkers=53
    )

    # Estimate variance of performances for base model with Rust Backend
    n_base <- 50
    nw <- min(half_cores(), n_base)
    aki_list <- list(aki)
    start <- Sys.time()
    logf("Running %d benchmarks with default params", n_base)
    baseperfs <- mcmapply(nw, run_aki_benchmark,
        aki=list(aki), seed=seq_len(n_base), kout=10, kin=10, fx="deconvolute",
        mt="lasso", norm="none", nfit=3, smit=2, smws=5,
        delta=6.4, npmax=0, maxShift=50, maxCombine=5,
        verbose=TRUE, nworkers=1,
        use_rust=TRUE
    )
    duration <- format(Sys.time() - start, units = "secs")
    logf("Finished %d benchmarks in %s seconds", n_base, duration)
    baseperfs <- rbindlist(baseperfs)
    # mean(baseperfs$acc)  # 0.73
    # range(baseperfs$acc) # 0.67 - 0.78

    # Repeat for R backend (use less repeats because of runtime)
    n_base <- 20
    nw <- min(half_cores(), n_base)
    aki_list <- list(aki)
    start <- Sys.time()
    logf("Running %d benchmarks with default params", n_base)
    baseperfsR <- mcmapply(1, run_aki_benchmark,
        aki=list(aki), seed=seq_len(n_base), kout=10, kin=10, fx="deconvolute",
        mt="lasso", norm="none", nfit=3, smit=3, smws=9,
        delta=6.4, npmax=0, maxShift=200, maxCombine=100,
        verbose=TRUE, nworkers=1,
        use_rust=FALSE
    )
    duration <- format(Sys.time() - start, units = "secs")
    logf("Finished %d benchmarks in %s seconds", n_base, duration)
    baseperfsR_acc <- sapply(baseperfsR, function(x) x$acc)
    mean(baseperfsR_acc) # 0.682
    range(baseperfsR_acc) # 0.641 - 0.707

    # First check which deconvolutions produce valid peaklists.
    D1 <- expand.grid(smit=c(2,3), smws=c(3,5,7), delta=2:8, nfit=5, npmax=0)
    D2 <- data.frame(smit=0, smws=0, delta=0, nfit=0, npmax=seq(500,2000,100))
    G <- as.data.frame(rbind(D1, D2))
    g <- nrow(G)
    nw <- 53
    smit <- G$smit; smws <- G$smws; delta <- G$delta; nfit <- G$nfit;
    npmax <- G$npmax;
    for (i in seq_len(nrow(G))) {
        logf("Checking deconvolution %d/%d: smit=%d, smws=%d, delta=%.1f, nfit=%d, npmax=%d",
            i, g, smit[i], smws[i], delta[i], nfit[i], npmax[]
        )
        decons <- try(deconvolute(
            aki$spectra, nfit=nfit[i], smopts=c(smit[i], smws[i]), delta=delta[i],
            sfr=NULL, wshw=0, ask=FALSE, force=FALSE, verbose=FALSE,
            nworkers=53, use_rust=TRUE, npmax=npmax[i]
        ))
        if (inherits(decons, "try-error")) next
        peaks <- sapply(decons, function(d) nrow(d$peak))
        mses <- sapply(decons, function(d) d$mse$raw)
        G[i, "npminobs"] <- min(peaks)
        G[i, "npmaxobs"] <- max(peaks)
        G[i, "npzero"] <- sum(peaks == 0)
        G[i, "idzero"] <- paste(which(peaks == 0), collapse = ",")
        G[i, "msemin"] <- min(mses)
        G[i, "msemax"] <- max(mses)
        G[i, "mseavg"] <- mean(mses)
        G[i, "msesd"] <- sd(mses)
    }
    Gsaved <- G

    # Run benchmarks for all parameter combinations in grid
    NZ <- G[G$npzero == 0, 1:5]
    A <- expand.grid(maxShift=2^(4:8), maxCombine=2^(4:8))
    P <- merge(NZ, A, by = NULL)
    g <- nrow(P)
    nw <- min(half_cores(), g)
    smit <- P$smit; smws <- P$smws; delta <- P$delta; nfit <- P$nfit;
    npmax <- P$npmax; maxShift <- P$maxShift; maxCombine <- P$maxCombine
    gridperfs <- mcmapply(10, try_aki_benchmark,
        aki=list(aki), seed=1, kout=10, kin=10, fx="deconvolute",
        mt="lasso", norm="none", nfit=nfit, smit=smit, smws=smws,
        delta=delta, npmax=npmax, maxShift=maxShift, maxCombine=maxCombine,
        verbose=FALSE, nworkers=1
    )
    acc <- sapply(gridperfs, function(x) x$acc)
    auc <- sapply(gridperfs, function(x) x$auc)
    P2 <- data.frame(P, acc = acc, auc = auc)
    bestacc <- which.max(P2$acc)
    bestauc <- which.max(P2$auc)
    bestavg <- which.max(P2$acc + P2$auc)
    Pbest <- P2[c(bestacc, bestauc, bestavg), ]
    #        smit smws delta nfit npmax maxShift maxCombine       acc       auc
    # 122       0    0     0    0  1000       64         16 0.8207547 0.8198529 (best complex model)
    # 1014      3    3     3    5     0      128        256 0.8018868 0.8864379 (best simple model)
    # 1014      3    3     3    5     0      128        256 0.8018868 0.8864379
    write.csv(P2, "gridperfs.csv", row.names = FALSE)

    # Estimate variance of performances of best "simple" model
    n_base <- 20
    nw <- min(half_cores(), n_base)
    aki_list <- list(aki)
    start <- Sys.time()
    logf("Running %d benchmarks with default params", n_base)
    bestperfs <- mcmapply(nw, run_aki_benchmark,
        aki=list(aki), seed=seq_len(n_base), kout=10, kin=10, fx="deconvolute",
        mt="lasso", norm="none", nfit=5, smit=3, smws=3,
        delta=3, npmax=0, maxShift=128, maxCombine=256,
        verbose=TRUE
    )
    duration <- format(Sys.time() - start, units = "secs")
    logf("Finished %d benchmarks in %s seconds", n_base, duration)
    bestaccs <- sapply(bestperfs, function(x) x$acc)
    bestaucs <- sapply(bestperfs, function(x) x$auc)
    mean(bestaccs)  # 0.797
    range(bestaccs) # 0.754 - 0.830

    # Estimate variance of performances of best "complex" model
    n_base <- 20
    nw <- min(half_cores(), n_base)
    aki_list <- list(aki)
    start <- Sys.time()
    logf("Running %d benchmarks with default params", n_base)
    bestperfs <- mcmapply(1, run_aki_benchmark,
        aki=list(aki), seed=seq_len(n_base), kout=10, kin=10, fx="deconvolute",
        mt="lasso", norm="none", nfit=5, smit=0, smws=0,
        delta=0, npmax=1000, maxShift=64, maxCombine=16,
        verbose=TRUE
    )
    duration <- format(Sys.time() - start, units = "secs")
    logf("Finished %d benchmarks in %s seconds", n_base, duration)
    bestaccs <- sapply(bestperfs, function(x) x$acc)
    bestaucs <- sapply(bestperfs, function(x) x$auc)
    mean(bestaccs)  # 0.784
    range(bestaccs) # 0.745 - 0.820
}

plot_deconvolution_metrics <- function() {
    aki <- read_aki_data()
    spectra <- aki$spectra
    G <- expand.grid(delta=2:8, smit=c(2,3), smws=c(3,5,7))
    rownames(G) <- sprintf("smit%d_smws%d_delta%d", G$smit, G$smws, G$delta)
    smopts_list <- lapply(seq_len(nrow(G)), function(i) c(G$smit[i], G$smws[i]))
    delta_list <- as.list(G$delta)
    spectra_list <- list(spectra)
    decons_list <- mcmapply(
        nrow(G), deconvolute, x=spectra_list, smopts=smopts_list, delta=delta_list,
        MoreArgs = list(
            nfit=5, sfr=NULL, wshw=0, ask=FALSE, force=FALSE, verbose=TRUE,
            nworkers=1, use_rust=TRUE, npmax=0
        )
    )
    decons_list_bak <- decons_list
    NP <- as.data.frame(matrix(NA, nrow(G), length(spectra)))
    colnames(NP) <- paste0("A", sapply(strsplit(get_names(spectra), "_"), "[", 4))
    rownames(NP) <- rownames(G)
    AR <- NP
    for (i in seq_len(nrow(G))) {
        decons <- decons_list[[i]]
        for (j in seq_along(decons)) {
            decon <- decons[[j]]
            y <- decon$si
            yhat <- decon$sit$sup
            AR[i, j] <- sum(abs(y - yhat)) / sum(abs(y))
            NP[i, j] <- nrow(decon$peak)
        }
    }
    meta <- aki$meta
    stopifnot(all(meta$sid == names(spectra)))
    ord <- order(meta$type)
    AR <- AR[, ord]
    NP <- NP[, ord]
    # hm <- function(M, file, title) {
    #     Z <- as.matrix(M)
    #     col <- grDevices::colorRampPalette(c("#FFF7BC", "#225EA8"))(100)
    #     grDevices::pdf(file, 11, 7)
    #     ht <- ComplexHeatmap::Heatmap(
    #         Z,
    #         name = title,
    #         col = col,
    #         cluster_rows = FALSE,
    #         cluster_columns = FALSE,
    #         row_names_gp = grid::gpar(fontsize = 8),
    #         column_names_gp = grid::gpar(fontsize = 6),
    #         column_names_rot = 90
    #     )
    #     ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
    #     grDevices::dev.off()
    # }
    # hm(AR, "AR_heatmap.pdf", "Area Ratio")
    # hm(NP, "NP_heatmap.pdf", "Number of Peaks")
}

plot_deconvolution_metrics_R <- function() {
    aki <- read_aki_data()
    spectra <- aki$spectra
    G <- expand.grid(delta=2:8, smit=c(2,3), smws=c(3,5,7))
    rownames(G) <- sprintf("smit%d_smws%d_delta%d", G$smit, G$smws, G$delta)
    smopts_list <- lapply(seq_len(nrow(G)), function(i) c(G$smit[i], G$smws[i]))
    delta_list <- as.list(G$delta)
    spectra_list <- list(spectra)
    decons_list <- mcmapply(
        nrow(G), deconvolute, x=spectra_list, smopts=smopts_list, delta=delta_list,
        MoreArgs = list(
            nfit=5, sfr=NULL, wshw=0, ask=FALSE, force=FALSE, verbose=TRUE,
            nworkers=1, use_rust=FALSE, npmax=0
        )
    )
    NP <- as.data.frame(matrix(NA, nrow(G), length(spectra)))
    colnames(NP) <- paste0("A", sapply(strsplit(get_names(spectra), "_"), "[", 4))
    rownames(NP) <- rownames(G)
    AR <- NP
    for (i in seq_len(nrow(G))) {
        decons <- decons_list[[i]]
        for (j in seq_along(decons)) {
            decon <- decons[[j]]
            y <- decon$si
            yhat <- decon$sit$sup
            AR[i, j] <- sum(abs(y - yhat)) / sum(abs(y))
            NP[i, j] <- nrow(decon$peak)
        }
    }
    meta <- aki$meta
    stopifnot(all(meta$sid == names(spectra)))
    ord <- order(meta$type)
    AR <- AR[, ord]
    NP <- NP[, ord]
    # hm <- function(M, file, title) {
    #     Z <- as.matrix(M)
    #     col <- grDevices::colorRampPalette(c("#FFF7BC", "#225EA8"))(100)
    #     grDevices::pdf(file, 11, 7)
    #     ht <- ComplexHeatmap::Heatmap(
    #         Z,
    #         name = title,
    #         col = col,
    #         cluster_rows = FALSE,
    #         cluster_columns = FALSE,
    #         row_names_gp = grid::gpar(fontsize = 8),
    #         column_names_gp = grid::gpar(fontsize = 6),
    #         column_names_rot = 90
    #     )
    #     ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
    #     grDevices::dev.off()
    # }
    write.csv(AR, "AR_R.csv")
    write.csv(NP, "NP_R.csv")
    # hm(AR, "AR_R_heatmap.pdf", "Area Ratio")
    # hm(NP, "NP_R_heatmap.pdf", "Number of Peaks")
}


# Helpers #####

map <- function(FUN, ..., MoreArgs = NULL) {
    FUN <- match.fun(FUN)
    mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = FALSE)
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

plot_probabilities <- function(Y) {
    plot_empty(
        xlim = c(1, nrow(Y)), ylim = c(0, 1),
        axes = TRUE, ylab = "Predicted probability",
        xlab = "Row Nr in Data.frame",
        main = "Lasso CV Predictions"
    )
    text(seq_along(Y$prob), Y$prob, labels = Y$snr, col = as.integer(Y$true))
}

#' @noRd
#' @title Plot Grid Performance Tree
#' @description
#' Fits and plots a simple decision tree to explain accuracy from tuning
#' parameters. Rows with `npmax != 0` are excluded.
plot_gridperf_tree <- function(
    infile = "gridperfs.csv",
    minsplit = 20,
    cp = 0.01
) {
    req <- c("acc", "smit", "smws", "delta", "nfit",
        "npmax", "maxShift", "maxCombine")
    P <- read.csv(infile, check.names = FALSE)
    miss <- setdiff(req, names(P))
    if (length(miss) > 0) {
        stop(sprintf("Missing columns: %s", paste(miss, collapse = ", ")))
    }
    P <- P[P$npmax == 0 & !is.na(P$acc), req, drop = FALSE]
    if (nrow(P) < 2) stop("Need at least two rows with npmax == 0")

    x <- P[, setdiff(req, c("acc", "npmax")), drop = FALSE]
    for (nm in names(x)) x[[nm]] <- as.factor(x[[nm]])
    df <- data.frame(acc = P$acc, x, check.names = FALSE)

    fit <- rpart::rpart(
        acc ~ .,
        data = df,
        method = "anova",
        control = rpart::rpart.control(minsplit = minsplit, cp = cp)
    )
    graphics::plot(fit, uniform = TRUE, margin = 0.05)
    graphics::text(fit, use.n = TRUE, cex = 0.7)
    invisible(fit)
}

# End-to-end Deconvolution Training #####

#' @noRd
#' @description
#' Fits Lasso models over a grid of lambda and decpar values,
#' with decpar is a string encoding all 'deconvolution parameters'
#' passed on to `decconvolute()`.
#' @return Returns a matrix with LASSO
grid_train_cv_lasso <- function() {
    aki <- read_aki_data()
    spectra_normed <- creatinine_normalize(aki$spectra)
    outdir <- file.path(datadir(), "grid_train_cv_lasso")
    grid_deconvolute_spectra(spectra_normed, outdir)
}

read_aki_data <- function() {
    # Read in data
    aki_path <- datadir("example_datasets/bruker/aki")
    if (!dir.exists(aki_path)) metabodecon::download_example_datasets()
    meta <- read_aki_metadata(aki_path)
    spectra_raw <- metabodecon::read_spectra(aki_path)
    stopifnot(all.equal(names(spectra_raw), meta$sid))

    # Normalize spectra
    spectra <- creatinine_normalize(spectra_raw)

    list(spectra = spectra, meta = meta)
}

build_param_grid <- function() {

    # Deconvolution grid (58 x 5)
    D1 <- expand.grid(smit=c(2,3), smws=c(3,5,7), delta=2:8, nfit=5, npmax=0)
    D2 <- data.frame(smit=0, smws=0, delta=0, nfit=0, npmax=seq(500,2000,100))
    D <- as.data.frame(rbind(D1, D2))

    # Alignment grid (25 x 2)
    A <- expand.grid(maxShift=2^(4:8), maxCombine=2^(4:8))

    # Combined grid (1450 x 7)
    merge(D, A, by = NULL)
}
