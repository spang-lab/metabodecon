
# API #####

#' @export
#' @name mdm
#' @rdname mdm
#'
#' @title Metabodecon Models
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' **WARNING: These functions are experimental and must not be used in
#' production. Their API is very likely to change in non-backwards-compatible
#' ways over the next few weeks.**
#'
#' Utilities for fitting, tuning and benchmarking 'metabodecon models' (mdm).
#'
#' A mdm is a essentially a glmnet lasso model, fitted on a feature matrix X,
#' obtained by deconvoluting and aligning spectra and snapping their peaks to a
#' shared reference spectrum. Lambda is chosen via cross-validation on X.
#' Deconvolution parameters (`npmax` or `nfit`/`smit`/`smws`/`delta`),
#' alignment parameters (`maxShift`, `maxCombine`) and the snapping parameter
#' (`maxSnap`) are tunable hyperparameters.
#'
#' When `npmax > 0`, deconvolution runs an internal grid search over smoothing
#' and fitting parameters and keeps at most `npmax` peaks. The explicit
#' `nfit`/`smit`/`smws`/`delta` values are ignored in that case. When
#' `npmax == 0`, the explicit parameters are used directly.
#'
#' [metabodecon::fit_mdm()] deconvolutes spectra, aligns detected peaks, snaps
#' them to a shared reference and fits one lasso model via
#' [glmnet::cv.glmnet()]. Its reported cross-validation error covers only the
#' model-fitting step on the already constructed feature matrix, so it is
#' slightly over-optimistic. Use [metabodecon::benchmark_mdm()] to get a
#' realistic performance estimate via nested cross-validation.
#'
#' [metabodecon::fit_mdm_grid()] calls [metabodecon::fit_mdm()] over a grid of
#' preprocessing parameter combinations and returns the model with the best
#' estimated performance.
#'
#' [metabodecon::benchmark_mdm()] wraps [metabodecon::fit_mdm_grid()] in an
#' outer cross-validation loop to estimate end-to-end predictive performance.
#' It returns the per-fold models and held-out predictions.
#'
#' @details
#'
#' ## Disk caching across processes
#'
#' [metabodecon::fit_mdm_grid()] and [metabodecon::benchmark_mdm()] can share
#' deconvolution results through directory `cadir`. This on-disk cache lets
#' parallel workers reuse expensive preprocessing results while keeping RAM use
#' manageable. Disk caching can be disabled by setting `cadir =
#' NULL`, but this will increase runtime drastically (from a few minutes/hours
#' to several days).
#'
#' ## RAM caching within processes
#'
#' [metabodecon::fit_mdm()] caches the most recent deconvolution and alignment
#' results in RAM when the input `spectra` carries a `"hash"` attribute. This
#' avoids redundant preprocessing when [metabodecon::fit_mdm_grid()] evaluates
#' several settings on the same spectra.
#'
#' ## Memory usage by parallel workers
#'
#' Each outer worker in [metabodecon::benchmark_mdm()] keeps its own copy of
#' `spectra` and `y`, so memory use grows roughly with `nwo`. Large `nwo` values
#' can therefore require substantial RAM for large input datasets.
#'
#' Inner workers (`nwi`) are used inside each outer worker for deconvolution,
#' alignment and prediction on held-out spectra. They process only the spectra
#' needed for the current task, so their memory impact is much smaller. In
#' RAM-constrained settings, use `nwo = 1` and increase `nwi` instead.
#'
#' The total amount of processes spawned is `nwo * nwi.`
#'
#' @param spectra List-like spectra object with `cs` and `si` vectors.
#' @param y Factor vector with class labels for each spectrum.
#' @param sfr Signal free region. See [metabodecon::deconvolute()] for details.
#' @param use_rust Logical. Whether to use the Rust backend.
#' @param npmax
#' Maximum number of peaks to retain. When `npmax > 0`, deconvolution runs an
#' internal grid search and the explicit `nfit`/`smit`/`smws`/`delta` are
#' ignored. Set to 0 to use explicit deconvolution parameters instead.
#' @param nfit Number of Lorentz-curve fitting iterations (used when `npmax == 0`).
#' @param smit Number of smoothing iterations (used when `npmax == 0`).
#' @param smws Smoothing window size (used when `npmax == 0`).
#' @param delta Peak-filter threshold (used when `npmax == 0`).
#' @param maxShift Maximum alignment shift.
#' @param maxCombine Maximum number of peaks to combine during alignment.
#' @param maxSnap Maximum snap distance in datapoints for feature-matrix construction.
#' @param nworkers
#' Number of workers used by [metabodecon::fit_mdm()] and
#' [metabodecon::fit_mdm_grid()] for deconvolution and alignment.
#' @param verbosity Integer. Verbosity level; each nested call
#' decrements by 1. Messages print when `verbosity >= 1`.
#' @param nfolds
#' Number of folds for the `cv.glmnet()` call in [metabodecon::fit_mdm()]
#' and [metabodecon::fit_mdm_grid()]. Default 10.
#' @param nfp
#' Number of folds used for repeated cross-validated performance estimation
#' inside [metabodecon::fit_mdm()]. Default 10.
#' @param cadir
#' Directory used to cache deconvolution results across grid points and
#' processes. Defaults to [metabodecon::decon_cachedir()].
#' @param nwo
#' Number of workers for the outer cross-validation in
#' [metabodecon::benchmark_mdm()]. Each outer worker holds a full copy of
#' `spectra` and `y`.
#' @param nwi
#' Number of workers used inside each outer fold for
#' [metabodecon::fit_mdm_grid()] and [metabodecon::predict.mdm()].
#' @param nfo Number of outer folds in [metabodecon::benchmark_mdm()].
#' @param nfl
#' Number of folds for the `cv.glmnet()` call inside
#' [metabodecon::benchmark_mdm()]. Passed as `nfolds` to
#' [metabodecon::fit_mdm()].
#' @param seed
#' Random seed used for fold assignment in [metabodecon::benchmark_mdm()] and
#' for inner cross-validation splits in [metabodecon::fit_mdm_grid()] and
#' [metabodecon::fit_mdm()].
#' @param check Logical. Whether to validate inputs at function entry.
#' @param pgrid Data frame of preprocessing parameter combinations as returned
#'   by [metabodecon::get_pgrid()].
#' @param ignore Optional integer vector of sample indices to exclude.
#' @param warm_cache Logical. Whether to pre-populate the disk cache before
#'   starting the grid search.
#' @param conf Character string selecting a predefined parameter grid
#'   configuration.
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with elements
#' `model`, `ref`, `meta`, `perf` (list with `by_seed` and `pooled`
#' data frames reporting accuracy and AUC) and `preds` (data frame with
#' columns `seed`, `fold`, `true`, `link`, `prob`, `pred`).
#'
#' [metabodecon::fit_mdm_grid()] returns an `mdm` object with an additional element
#' `pgrid` containing the performance grid.
#'
#' [metabodecon::benchmark_mdm()] returns a list with elements:
#' - `models`: List of `mdm` objects, one per outer fold.
#' - `predictions`: Data frame with columns `fold`, `true`, `link`, `prob`,
#'   `pred`.
#'
#' @examples
#' \dontrun{
#'
#'
#'
#'      # -~-~-~ Inputs -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#'
#'      aki <- read_aki_data()
#'      spectra <- aki$spectra
#'      attr(spectra, "hash") <- rlang::hash(spectra)
#'      y <- factor(aki$meta$type, levels = c("Control", "AKI"))
#'      names(y) <- rownames(aki$meta)
#'      sfr <- c(11, -2)
#'      cadir <- cachedir("deconvs", persistent = TRUE)
#'
#'
#'
#'      # -~-~-~-~-~-~-~-~-~ Single -~-~-~-~-~-~-~-~-~
#'
#'      # Best "simple" model from aki.R (0.797 [0.75 - 0.83])
#'      mdm <- fit_mdm(
#'          spectra, y, sfr=NULL,
#'          nfit=5, smit=3, smws=3, delta=3, npmax=0,
#'          maxShift=128, maxCombine=256, maxSnap=0,
#'          nworkers=half_cores(), cadir=cadir
#'      )
#'      # Per-seed acc: [77.36%, 82.08%]
#'      # Per-seed AUC: [0.7884, 0.8656]
#'      # Pooled acc: 79.34% +- 1.75%
#'      # Pooled AUC: 0.8235 +- 0.0244
#'
#'      # Best "complex" model from aki.R (0.784, [0.74 - 0.82])
#'      mdm <- fit_mdm(
#'          spectra, y, sfr = NULL,
#'          nfit=5, smit=0, smws=0, delta=0, npmax=1000,
#'          maxShift=64, maxCombine=16, maxSnap=0,
#'          nworkers=half_cores(), cadir=cadir
#'      )
#'      # Pooled acc: 74.81% +- 2.56%
#'      # Pooled AUC: 0.7764 +- 0.023
#'
#'
#'
#'      # -~-~-~-~-~-~-~-~-~ Search -~-~-~-~-~-~-~-~-~
#'
#'      mdm_grid_stat1 <- fit_mdm_grid(
#'          spectra, y, pgrid=get_pgrid("static1"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'      saveRDS(mdm_grid_stat1, "tmp/mdm_grid_stat1.rds")
#'
#'
#'      mdm_grid_stat2 <- fit_mdm_grid(
#'          spectra, y, pgrid=get_pgrid("static2"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'      saveRDS(mdm_grid_stat2, "tmp/mdm_grid_stat2.rds")
#'      #
#'      # Best static2: acc=0.82, auc=0.89
#'      # smit=2, smws=3, delta=6, nfit=7, npmax=0, maxShift=100, maxCombine=0, maxSnap=30
#'
#'
#'      mdm_grid_stat3 <- fit_mdm_grid(
#'          spectra, y, pgrid=get_pgrid("static3"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'      saveRDS(mdm_grid_stat3, "tmp/mdm_grid_stat3.rds")
#'      #
#'      # Best static3: acc=79.34% auc=0.8589
#'      # smit=2, smws=7, delta=8, nfit=10, npmax=0, maxShift=100, maxCombine=0, maxSnap=30
#'
#'
#'      mdm_grid_dyn2 <- fit_mdm_grid(
#'          spectra, y, pgrid=get_pgrid("dynamic2"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'
#'      # -~-~-~ Full Benchmark -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#'      bm <- benchmark_mdm(spectra, y, sfr, cadir=cadir)
#'
#'      # -~-~-~ Interactive Development -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#'      stub(fit_mdm, spectra=spectra, y=y, sfr=sfr)
#'      stub(fit_mdm_grid, spectra=spectra, y=y, sfr=sfr)
#'      stub(benchmark_mdm, spectra=spectra, y=y, sfr=sfr)
#' }
#'
fit_mdm <- function(
    # Mandatory
    spectra, y,
    # Shared
    sfr = NULL, use_rust = TRUE, nworkers = 1, verbosity = 2,
    seed = 1, cadir = decon_cachedir(), check = TRUE,
    # Deconvolution
    npmax = 1000, nfit = 5, smit = 2, smws = 5, delta = 6.4,
    # Alignment
    maxShift = 200, maxCombine = 0,
    # Snapping
    maxSnap = 50,
    # Glmnet
    nfolds = 10, nfp = 10
) {

    if (check) check_mdm_args(
        spectra = spectra, y = y, sfr = sfr,
        use_rust = use_rust, nworkers = nworkers,
        verbosity = verbosity, seed = seed,
        cadir = cadir, check = check,
        nfolds = nfolds, nfp = nfp,
        npmax = npmax, nfit = nfit, smit = smit,
        smws = smws, delta = delta,
        maxShift = maxShift, maxCombine = maxCombine,
        maxSnap = maxSnap
    )

    # Cache keys
    hash <- attr(spectra, "hash")
    decons_key <- list(
        hash = hash, sfr = sfr, use_rust = use_rust,
        npmax = npmax, nfit = nfit, smit = smit,
        smws = smws, delta = delta
    )
    aligns_key <- c(decons_key, list(
        maxShift = maxShift, maxCombine = maxCombine
    ))
    cache <- getOption("metabodecon.fit_mdm.cache")
    dc <- !is.null(hash) && identical(cache$decons_key, decons_key)
    ac <- !is.null(hash) && identical(cache$aligns_key, aligns_key)

    fmt <- "Deconvoluting (npmax=%d, nfit=%d, smit=%d, smws=%d, delta=%.1f, cached=%s)"
    logv(fmt, npmax, nfit, smit, smws, delta, dc)
    decons <- if (dc) cache$decons else deconvolute(
        x = spectra, sfr = sfr, verbose = as.logical(verbosity - 1),
        use_rust = use_rust, nfit = nfit,
        smopts = c(smit, smws), delta = delta,
        npmax = npmax, nworkers = nworkers, cadir = cadir
    )
    if (!is.null(hash) && !dc) {
        options(metabodecon.fit_mdm.cache = list(
            decons_key = decons_key, decons = decons
        ))
    }
    nps <- unname(unlist(lapply(decons, function(decon) nrow(decon$lcpar))))
    npzero <- sum(nps == 0)
    if (npzero > 0) {
        logv("Aborting: %d/%d spectra produced zero peaks", npzero, length(nps))
        meta <- list(
            sfr = sfr, use_rust = use_rust, npmax = npmax, nfit = nfit,
            smit = smit, smws = smws, delta = delta,
            maxShift = maxShift, maxCombine = maxCombine, maxSnap = maxSnap
        )
        pooled <- data.frame(acc = 0, acc_sd = 0, auc = 0, auc_sd = 0)
        perf <- list(by_seed = NULL, pooled = pooled)
        mdm <- list(model = NULL, ref = NULL, meta = meta,
                    perf = perf, preds = NULL)
        return(structure(mdm, class = "mdm"))
    }
    nprange <- range(nps)
    logv("Peak range: [%d, %d]", nprange[1], nprange[2])

    fmt <- "Aligning (maxShift=%d, maxCombine=%d, cached=%s)"
    logv(fmt, maxShift, maxCombine, ac)
    aligns <- if (ac) cache$aligns else align(
        x = decons, maxShift = maxShift,
        maxCombine = maxCombine, verbose = as.logical(verbosity - 1),
        nworkers = nworkers, method = 2, full = FALSE
    )

    if (!is.null(hash)) {
        options(metabodecon.fit_mdm.cache = list(
            decons_key = decons_key, aligns_key = aligns_key,
            decons = decons, aligns = aligns
        ))
    }

    logv("Constructing feature matrix (maxSnap=%g)", maxSnap)
    pl <- lapply(aligns, get_peak_indices)
    idx <- find_ref(pl)$refInd
    ref <- aligns[[idx]]
    X <- t(get_si_mat(aligns, maxSnap = maxSnap, ref = ref, drop_zero = TRUE))

    logv("Fitting cv.glmnet")
    foldid <- get_foldid(y = y, nfolds = nfolds, seed = seed)
    model <- glmnet::cv.glmnet(
        X, y, family = "binomial", alpha = 1, foldid = foldid
    )

    logv("Estimating performance (roughly)")
    bm <- benchmark_cv_glmnet(
        X = X, y = y, nfolds = nfp, seed = seed, verbosity = verbosity
    )

    logv("Constructing and returning mdm object")
    meta <- list(
        sfr = sfr, use_rust = use_rust,
        npmax = npmax, nfit = nfit, smit = smit,
        smws = smws, delta = delta,
        maxShift = maxShift, maxCombine = maxCombine,
        maxSnap = maxSnap
    )
    mdm <- list(
        model = model, ref = ref, meta = meta,
        perf = bm$perf, preds = bm$preds
    )
    structure(mdm, class = "mdm")
}

#' @export
#' @rdname mdm
fit_mdm_grid <- function(
    # Mandatory
    spectra, y, sfr = NULL,
    # Shared
    use_rust = TRUE, nworkers = 1, verbosity = 2,
    seed = 1, cadir = decon_cachedir(), check = TRUE,
    # Grid
    pgrid = get_pgrid("dynamic2"),
    ignore = NULL, warm_cache = TRUE,
    # Model (passed through to fit_mdm)
    nfolds = 10, nfp = 10
) {

    if (check) check_mdm_args(
        spectra = spectra, y = y, sfr = sfr,
        use_rust = use_rust, nworkers = nworkers,
        verbosity = verbosity, nfolds = nfolds, nfp = nfp,
        seed = seed, cadir = cadir, pgrid = pgrid,
        ignore = ignore, check = check
    )

    if (!is.null(ignore)) {
        logv("Removing %d ignored samples from spectra and y", length(ignore))
        spectra <- spectra[-ignore]
        y <- y[-ignore]
    }

    if (warm_cache) {
        init_cache(
            x = spectra, sfr = sfr,
            verbosity = verbosity - 1,
            nworkers = nworkers, use_rust = use_rust,
            cadir = cadir, pgrid = pgrid
        )
    }

    np <- nrow(pgrid)
    ns <- length(spectra)
    logv("Starting grid search (%d combinations, %s spectra)", np, ns)
    on.exit(options(metabodecon.fit_mdm.cache = NULL), add = TRUE)
    attr(spectra, "hash") <- rlang::hash(spectra)

    pp_names <- c(
        "npmax", "nfit", "smit", "smws", "delta",
        "maxShift", "maxCombine", "maxSnap"
    )
    pp_cols <- intersect(pp_names, names(pgrid))

    best_mdm <- NULL
    best_acc <- -Inf
    best_auc <- NA_real_
    for (i in seq_len(np)) {
        pc <- as.list(pgrid[i, pp_cols, drop = FALSE])
        tag <- paste(
            sprintf("%s=%g", pp_cols, unlist(pc)),
            collapse = ", "
        )
        logv("[%d/%d]: %s", i, np, tag)
        mdm <- do.call(fit_mdm, c(
            list(
                spectra = spectra, y = y, sfr = sfr,
                use_rust = use_rust, nworkers = nworkers,
                verbosity = verbosity - 1, seed = seed,
                cadir = cadir, check = FALSE,
                nfolds = nfolds, nfp = nfp
            ),
            pc
        ))
        pgrid$acc[i] <- mdm$perf$pooled$acc
        pgrid$acc_sd[i] <- mdm$perf$pooled$acc_sd
        pgrid$auc[i] <- mdm$perf$pooled$auc
        pgrid$auc_sd[i] <- mdm$perf$pooled$auc_sd
        if (pgrid$acc[i] > best_acc) {
            best_acc <- pgrid$acc[i]
            best_auc <- pgrid$auc[i]
            best_mdm <- mdm
        }
        logv("[%d/%d] acc=%.2f%% auc=%.4f | best: acc=%.2f%% auc=%.4f",
            i, np, pgrid$acc[i] * 100, pgrid$auc[i],
            best_acc * 100, best_auc)
    }
    ibest <- which.max(pgrid$acc)
    best <- as.list(pgrid[ibest, ])
    logv("Best [%d/%d]: acc=%.2f%% auc=%.4f", ibest, np,
        best$acc * 100, best$auc)
    best_mdm$pgrid <- pgrid
    best_mdm
}

#' @export
#' @rdname mdm
benchmark_mdm <- function(
    # Mandatory
    spectra, y, sfr = NULL,
    # Shared
    use_rust = TRUE, verbosity = 2,
    seed = 1, cadir = decon_cachedir(),
    # Workers
    nwo = 1, nwi = half_cores(),
    # Grid
    pgrid = get_pgrid(),
    # CV folds
    nfo = 5, nfl = 10, nfp = 10
) {
    check_mdm_args(
        spectra=spectra, y=y, sfr=sfr, use_rust=use_rust,
        verbosity=verbosity, nfolds=nfo, nfl=nfl, nfp=nfp,
        nwo=nwo, nwi=nwi, cadir=cadir, pgrid=pgrid, seed=seed
    )

    ns <- length(spectra)
    init_cache(
        x=spectra, sfr=sfr, verbosity=verbosity - 1,
        nworkers=nwi * nwo, use_rust=use_rust, cadir=cadir,
        pgrid=pgrid
    )

    logv("Fitting fit_mdm_grid for %d outer folds", nfo)
    te_list <- get_test_ids(nfolds=nfo, nsamples=ns, seed=seed, y=y)
    fids <- seq_along(te_list)
    nwoa <- min(nwo, nfo)
    mdms <- mcmapply(
        nwoa, fit_mdm_grid,
        ignore=te_list,
        MoreArgs=list(
            spectra=spectra, y=y, sfr=sfr, use_rust=use_rust,
            nworkers=nwi, verbosity=verbosity - 1,
            nfolds=nfl, nfp=nfp, seed=seed, cadir=cadir,
            pgrid=pgrid, warm_cache=FALSE, check=FALSE
        )
    )

    logv("Predicting on held-out test folds")
    te_spectra <- lapply(te_list, function(te) spectra[te])
    pred_metrics <- mcmapply(
        nwoa, predict.mdm,
        object=mdms, newdata=te_spectra,
        MoreArgs=list(type="all", nworkers=nwi),
        SIMPLIFY=FALSE, USE.NAMES=FALSE
    )
    preds <- mapply(
        FUN=function(fid, te, pred) {
            data.frame(fold=fid, true=y[te], link=pred$link, prob=pred$prob, pred=pred$class)
        },
        fids, te_list, pred_metrics,
        SIMPLIFY=FALSE, USE.NAMES=FALSE)

    Y <- do.call(rbind, preds)
    acc <- mean(Y$true == Y$pred, na.rm = TRUE)
    logv("Nested CV accuracy: %.2f%%", acc * 100)
    out <- list(models=mdms, predictions=Y)
    out$mdms <- out$models
    out
}

#' @export
#' @rdname mdm
get_pgrid <- function(
    conf = "dynamic2" # "static1", "dynamic1", "static2", "dynamic2"
) {
    if (conf == "static1") {
        P <- expand.grid(
            smit = c(2, 3),
            smws = c(3, 5, 7, 9),
            delta = 4:8,
            nfit = 5,
            npmax = 0,
            maxShift = 2^(4:8),
            maxCombine = 2^(4:8),
            maxSnap = 0,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else if (conf == "dynamic1") {
        P <- data.frame(
            smit = 0,
            smws = 0,
            delta = 0,
            nfit = 0,
            npmax = seq(800, 1200, 100),
            maxShift = 2^(4:8),
            maxCombine = 2^(4:8),
            maxSnap = 0,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else if (conf == "static2") {
        P <- expand.grid(
            smit = c(2, 3),
            smws = c(3, 5, 7, 9),
            delta = 4:8,
            nfit = c(3, 5, 7),
            npmax = 0,
            maxShift = seq(50, 250, 50),
            maxCombine = 0,
            maxSnap = seq(10, 50, 10),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else if (conf == "static3") {
        P <- expand.grid(
            smit = c(2),
            smws = c(3, 5, 7, 9),
            delta = 4:8,
            nfit = c(10),
            npmax = 0,
            maxShift = seq(50, 250, 50),
            maxCombine = 0,
            maxSnap = seq(10, 50, 10),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else if (conf == "dynamic2") {
        P <- expand.grid(
            nfit = 0,
            smit = 0,
            smws = 0,
            delta = 0,
            npmax = seq(800, 1200, 100),
            maxShift = seq(50, 250, 50),
            maxCombine = 0,
            maxSnap = seq(10, 50, 10),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else {
        stop("Invalid conf")
    }
    ord <- order(
        P$npmax, P$nfit, P$smit, P$smws, P$delta,
        P$maxShift, P$maxCombine, P$maxSnap
    )
    P <- P[ord, ]
    P$acc <- NA_real_
    P$acc_sd <- NA_real_
    P$auc <- NA_real_
    P$auc_sd <- NA_real_
    rownames(P) <- NULL
    P
}

# Helpers #####

init_cache <- function(
    x, sfr, pgrid = get_pgrid(), nworkers = 1,
    verbosity = 1, use_rust = TRUE, cadir = decon_cachedir()
) {

    ns <- length(x)
    nwu <- min(nworkers, ns)

    # npmax-based rows need the grid-decon cache + per-npmax cache
    npmaxs <- if ("npmax" %in% names(pgrid)) {
        unique(pgrid$npmax[pgrid$npmax > 0])
    } else {
        integer(0)
    }

    if (length(npmaxs) > 0) {
        logv("Initializing grid cache for (ns=%d, nworkers=%d)", ns, nwu)
        grid_deconvolute_spectra(
            x = x, sfr = sfr, verbose = as.logical(verbosity - 1),
            nw = nwu, use_rust = use_rust, cadir = cadir
        )
        logv("Initializing npmax cache for %d values", length(npmaxs))
        cache <- disk_cache(cadir)
        for (npmax in npmaxs) {
            arglist <- lapply(x, function(s) {
                list("ds", s, sfr, 0, FALSE, 2, use_rust, npmax, "decon2")
            })
            hashes <- mapply(rlang::hash, arglist, USE.NAMES = FALSE)
            logv("Checking npmax cache for npmax=%d", npmax)
            seen <- sapply(hashes, cache$exists, USE.NAMES = FALSE)
            fmt <- "Reading %d/%d npmax=%d results from cache"
            logv(fmt, sum(seen), length(seen), npmax)
            if (all(seen)) next
            xu <- x[!seen]
            nwu_i <- min(nwu, length(xu))
            fmt <- "Deconvoluting %d unseen spectra (npmax=%d, nworkers=%d)"
            logv(fmt, length(xu), npmax, nwu_i)
            deconvolute(
                x = xu, sfr = sfr, wshw = 0,
                verbose = as.logical(verbosity - 1), nworkers = nwu_i,
                use_rust = use_rust, cadir = cadir,
                npmax = npmax
            )
        }
    }

    logv("Finished cache initialization")
    invisible(NULL)
}

benchmark_cv_glmnet <- function(
    X, y, nfolds, seed, nseeds = 10, verbosity = 1
) {
    n <- length(y)
    ns <- "lambda.min"
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max, nseeds)

    # Pre-allocate predictions (nseeds * n rows)
    nr <- nseeds * n
    sn <- names(y) %||% as.character(seq_len(n))
    preds <- data.frame(
        seed = integer(nr),
        fold = integer(nr),
        name = character(nr),
        true = factor(rep(NA, nr), levels = levels(y)),
        link = numeric(nr),
        prob = numeric(nr),
        pred = factor(rep(NA, nr), levels = levels(y)),
        stringsAsFactors = FALSE
    )

    nf <- nseeds * nfolds
    logv("Running %d CV fits (%d seeds x %d folds)", nf, nseeds, nfolds)
    row <- 0L
    by_seed <- data.frame(seed = seeds, acc = numeric(nseeds), auc = numeric(nseeds))
    for (si in seq_along(seeds)) {
        s <- seeds[si]
        te_list <- get_test_ids(nfolds, n, seed = s, y = y)
        for (fi in seq_len(nfolds)) {
            te <- te_list[[fi]]
            nte <- length(te)
            Xtr <- X[-te, , drop = FALSE]
            Xte <- X[te, , drop = FALSE]
            yte <- y[te]
            fm <- glmnet::cv.glmnet(Xtr, y[-te], family = "binomial", alpha = 1)
            idx <- (row + 1):(row + nte)
            preds$seed[idx] <- s
            preds$fold[idx] <- fi
            preds$name[idx] <- sn[te]
            preds$true[idx] <- yte
            preds$link[idx] <- predict(fm, Xte, s = ns, type = "link")[, 1]
            preds$prob[idx] <- predict(fm, Xte, s = ns, type = "response")[, 1]
            preds$pred[idx] <- predict(fm, Xte, s = ns, type = "class")[, 1]
            row <- row + nte
        }
        d <- preds[preds$seed == s, ]
        by_seed$acc[si] <- mean(d$true == d$pred)
        by_seed$auc[si] <- AUC(d$true, d$prob)
        fmt <- "seed %d/%d: acc=%.2f%%, AUC=%.4f"
        logvv(fmt, si, nseeds, by_seed$acc[si] * 100, by_seed$auc[si])
    }

    # CV summary
    pooled <- data.frame(
        acc = mean(preds$true == preds$pred), acc_sd = stats::sd(by_seed$acc),
        auc = AUC(preds$true, preds$prob), auc_sd = stats::sd(by_seed$auc)
    )
    logv("Per-seed acc: [%.2f%%, %.2f%%]", min(by_seed$acc) * 100, max(by_seed$acc) * 100)
    logv("Per-seed AUC: [%.4f, %.4f]", min(by_seed$auc), max(by_seed$auc))
    logv("Pooled acc: %.2f%% +- %.2f%%", pooled$acc * 100, pooled$acc_sd * 100)
    logv("Pooled AUC: %.4f +- %.4f", pooled$auc, pooled$auc_sd)
    list(preds = preds, perf = list(by_seed = by_seed, pooled = pooled))
}

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

get_foldid <- function(y, nfolds = 5, seed = 1) {
    te_list <- get_test_ids(
        nfolds = nfolds, nsamples = length(y), seed = seed, y = y
    )
    foldid <- integer(length(y))
    for (i in seq_along(te_list)) foldid[te_list[[i]]] <- i
    foldid
}

#' @noRd
#' @title Compute rank-based AUC
#' @description Computes area under the ROC curve using rank statistics.
#' @param y Binary labels coded as 0/1 or coercible to integer.
#' @param yhat Numeric prediction scores.
#' @return Numeric scalar AUC or `NA_real_` if one class is missing.
#' @examples
#' y <- c(0, 0, 1, 1)
#' yhat <- c(0.1, 0.3, 0.6, 0.8)
AUC <- function(y, yhat) {
    y <- as_binary01(y)
    pos <- y == 1
    n1 <- sum(pos)
    n0 <- sum(!pos)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(yhat)
    (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

check_mdm_args <- function(
    spectra, y, sfr,
    use_rust = NULL, nworkers = NULL, verbosity = NULL,
    nfolds = NULL, nfl = NULL, nfp = NULL, nwo = NULL, nwi = NULL,
    cadir = NULL, pgrid = NULL,
    ignore = NULL, seed = NULL, check = NULL,
    npmax = NULL, nfit = NULL, smit = NULL, smws = NULL, delta = NULL,
    maxShift = NULL, maxCombine = NULL, maxSnap = NULL
) {
    pp_names <- c(
        "npmax", "nfit", "smit", "smws", "delta",
        "maxShift", "maxCombine", "maxSnap"
    )
    stopifnot(
        is_spectra(spectra),
        is.factor(y),
        length(y) == length(spectra),
        is_num_or_null(sfr, 2),
        is_int_or_null(npmax, 1),
        is_int_or_null(nfit, 1),
        is_int_or_null(smit, 1),
        is_int_or_null(smws, 1),
        is_num_or_null(delta, 1),
        is_int_or_null(maxShift, 1),
        is_int_or_null(maxCombine, 1),
        is_num_or_null(maxSnap, 1),
        is_bool_or_null(use_rust, 1),
        is_int_or_null(nworkers, 1),
        is_int_or_null(verbosity, 1),
        is.null(nfolds) || (is_int(nfolds, 1) && nfolds >= 2),
        is.null(nfl) || (is_int(nfl, 1) && nfl >= 2),
        is.null(nfp) || (is_int(nfp, 1) && nfp >= 2),
        is_int_or_null(nwo, 1),
        is_int_or_null(nwi, 1),
        is_str_or_null(cadir),
        is.null(pgrid) || (
            is.data.frame(pgrid) &&
            nrow(pgrid) >= 1 &&
            all(pp_names %in% names(pgrid))
        ),
        is_int_or_null(ignore),
        is_int_or_null(seed, 1),
        is_bool_or_null(check, 1)
    )
    if (!is.null(names(y)) && !identical(get_names(spectra), names(y))) {
        stop("Names of `spectra` and `y` must match and be in the same order.", call. = FALSE)
    }
    if (nlevels(y) != 2 || any(table(y) == 0)) {
        stop("`y` must contain exactly 2 non-empty classes.", call. = FALSE)
    }
    if (!is.null(nfolds) && nfolds > length(y)) {
        stop("`nfolds` must not exceed the number of samples.", call. = FALSE)
    }
    if (!is.null(nfl) && nfl > length(y)) {
        stop("`nfl` must not exceed the number of samples.", call. = FALSE)
    }
    if (!is.null(nfp) && nfp > length(y)) {
        stop("`nfp` must not exceed the number of samples.", call. = FALSE)
    }
    if (!is.null(nwo) && nwo < 1) {
        stop("`nwo` must be at least 1.", call. = FALSE)
    }
    if (!is.null(nwi) && nwi < 1) {
        stop("`nwi` must be at least 1.", call. = FALSE)
    }
    if (!is.null(ignore) && any(!ignore %in% seq_along(spectra))) {
        stop("`ignore` contains invalid sample indices.", call. = FALSE)
    }
    invisible(NULL)
}

# S3 methods #####

#' @export
#' @name mdm_methods
#' @rdname mdm_methods
#'
#' @title S3 methods for mdm objects
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' **WARNING: These methods are experimental and must not be used in
#' production. Their API is very likely to change in non-backwards-compatible
#' ways over the next few weeks.**
#'
#' S3 methods for objects of class `mdm` and `summary.mdm`.
#'
#' `predict.mdm()` predicts probabilities, classes, link scores, or all
#' three from an `mdm` object. When `newdata` is a spectra object, the
#' spectra are deconvoluted, aligned and snapped to the reference stored in
#' the model before prediction. When `newdata` is a numeric matrix, it is
#' used directly as the feature matrix.
#'
#' `print.mdm()` prints a compact model summary.
#'
#' `coef.mdm()` returns lasso coefficients.
#'
#' `plot.mdm()` plots the lasso path.
#'
#' `summary.mdm()` builds a compact summary list.
#'
#' `print.summary.mdm()` prints formatted output for `summary.mdm` objects.
#'
#' @param object,x
#' A fitted `mdm` object (for `predict`, `coef`, `summary`, `print` and `plot`)
#' or a `summary.mdm` object (for `print.summary.mdm`).
#' @param newdata Spectra object or numeric feature matrix.
#' @param type Prediction type, one of `"all"`, `"prob"`, `"class"`, `"link"`.
#' @param s Regularization value for lasso predictions.
#' @param nworkers Number of workers to deconvolute and align `newdata`.
#' @param verbosity Integer verbosity level.
#' @param ... Passed to underlying methods where applicable.
#'
#' @return
#' - `predict`: numeric vector of probabilities, classes, and/or link scores.
#' - `print`: invisibly returns `x`.
#' - `coef`: coefficient object from `glmnet`.
#' - `plot`: invisibly returns `NULL`.
#' - `summary`: object of class `summary.mdm`.
#' - `print.summary.mdm`: invisibly returns `x`.
#'
#' @examples
#' m <- structure(
#'   list(
#'     model = NULL,
#'     ref = NULL,
#'     meta = list(npmax = 1000, nfit = 3, smit = 2, smws = 5,
#'                 delta = 6.4, maxShift = 100, maxCombine = 0,
#'                 maxSnap = 20)
#'   ),
#'   class = "mdm"
#' )
#' print(m)
#' summary(m)
#'
#' \dontrun{
#'   m <- fit_mdm_grid(spectra, y, sfr = c(11, -2))
#'   predict(m, test_spectra, type = "prob")
#'   coef(m)
#'   plot(m)
#' }
predict.mdm <- function(object,
                        newdata,
                        type = c("all", "prob", "class", "link"),
                        s = "lambda.min",
                        nworkers = 1,
                        verbosity = 1,
                        ...) {
    stopifnot(
        inherits(object, "mdm"), is_int(nworkers, 1),
        is_num(s, 1) || is_str(s), is_spectra(newdata) || is.matrix(newdata) || is.data.frame(newdata)
    )
    type <- match.arg(type)
    if (is.null(object$model)) {
        n <- if (is_spectra(newdata)) length(newdata)
             else nrow(as.matrix(newdata))
        z <- rep(0, n); h <- rep(0.5, n)
        cl <- factor(rep(NA_character_, n))
        if (type == "all") return(data.frame(link=z, prob=h, class=cl))
        if (type == "class") return(cl)
        if (type == "prob") return(h)
        return(z)
    }
    if (is_spectra(newdata)) {
        m <- object$meta
        logf("Deconvoluting %d spectra with %d nworkers",
             length(newdata), nworkers)
        if (m$npmax > 0) {
            decons <- deconvolute(
                x = newdata, sfr = m$sfr,
                verbose = as.logical(verbosity - 1),
                use_rust = m$use_rust, npmax = m$npmax,
                nworkers = nworkers
            )
        } else {
            decons <- deconvolute(
                x = newdata, sfr = m$sfr,
                verbose = as.logical(verbosity - 1),
                use_rust = m$use_rust, nfit = m$nfit,
                smopts = c(m$smit, m$smws),
                delta = m$delta, npmax = 0,
                nworkers = nworkers
            )
        }
        logf("Aligning spectra with %d nworkers", nworkers)
        als <- align(
            decons, maxShift = m$maxShift,
            maxCombine = m$maxCombine,
            verbose = as.logical(verbosity - 1), nworkers = nworkers,
            ref = object$ref, method = 2, full = FALSE
        )
        Xn <- t(get_si_mat(
            als, maxSnap = m$maxSnap, ref = object$ref,
            drop_zero = TRUE
        ))
    } else {
        Xn <- as.matrix(newdata)
    }
    logf("Predicting with s=%s", as.character(s))
    requireNamespace("glmnet", quietly = TRUE)
    score <- as.numeric(predict(object$model, newx = Xn, s = s, type = "link"))
    prob <- as.numeric(predict(object$model, newx = Xn, s = s, type = "response"))
    pred <- predict(object$model, newx = Xn, s = s, type = "class")[, 1]
    pred <- factor(pred, levels = levels(object$meta$y))
    if (type == "all") return(data.frame(link = score, prob = prob, class = pred))
    if (type == "class") return(pred)
    if (type == "prob") return(prob)
    if (type == "link") return(score)
    prob
}

#' @export
#' @rdname mdm_methods
print.mdm <- function(x, ...) {
    stopifnot(inherits(x, "mdm"), is.list(x$meta))
    pp <- c("npmax", "nfit", "smit", "smws", "delta",
            "maxShift", "maxCombine", "maxSnap")
    cat("metabodecon model (mdm)\n")
    for (nm in pp) {
        lab <- formatC(paste0(nm, ":"), width = -15)
        cat("  ", lab, x$meta[[nm]], "\n", sep = "")
    }
    if (!is.null(x$pgrid)) {
        cat("  grid rows:     ", nrow(x$pgrid), "\n", sep = "")
    }
    invisible(x)
}

#' @export
#' @rdname mdm_methods
coef.mdm <- function(object, ...) {
    stopifnot(inherits(object, "mdm"), !is.null(object$model))
    stats::coef(object$model, s = "lambda.min", ...)
}

#' @export
#' @rdname mdm_methods
plot.mdm <- function(x, ...) {
    stopifnot(inherits(x, "mdm"), !is.null(x$model))
    graphics::plot(x$model, ...)
    invisible(NULL)
}

#' @export
#' @rdname mdm_methods
summary.mdm <- function(object, ...) {
    stopifnot(inherits(object, "mdm"), is.list(object$meta))
    pp <- c("npmax", "nfit", "smit", "smws", "delta",
            "maxShift", "maxCombine", "maxSnap")
    out <- object$meta[pp]
    out$n_peaks <- if (is.null(object$ref)) 0L
                   else length(get_peak_indices(object$ref))
    out$grid_rows <- if (is.null(object$pgrid)) 0L
                     else nrow(object$pgrid)
    class(out) <- "summary.mdm"
    out
}

#' @export
#' @rdname mdm_methods
print.summary.mdm <- function(x, ...) {
    stopifnot(inherits(x, "summary.mdm"))
    cat("Summary of mdm\n")
    for (nm in names(x)) {
        lab <- formatC(paste0(nm, ":"), width = -15)
        cat("  ", lab, x[[nm]], "\n", sep = "")
    }
    invisible(x)
}
