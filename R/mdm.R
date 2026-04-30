
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
#' alignment parameter `maxShift` and peak-combining parameter `maxCombine`
#' are tunable hyperparameters.
#'
#' When `npmax > 0`, deconvolution runs an internal grid search over smoothing
#' and fitting parameters and keeps at most `npmax` peaks. The explicit
#' `nfit`/`smit`/`smws`/`delta` values are ignored in that case. When
#' `npmax == 0`, the explicit parameters are used directly.
#'
#' [metabodecon::fit_mdm()] deconvolutes spectra, aligns detected peaks,
#' combines them via [metabodecon::get_si_mat()] and fits one lasso model via
#' [glmnet::cv.glmnet()]. Lambda is selected internally by cross-validation
#' but no performance metrics are reported. Use [metabodecon::cv_mdm()] to
#' tune preprocessing parameters and [metabodecon::benchmark_mdm()] for
#' unbiased performance estimates.
#'
#' [metabodecon::cv_mdm()] evaluates a grid of preprocessing parameter
#' combinations. For each grid row it builds a feature matrix, runs
#' [glmnet::cv.glmnet()] with a fixed fold assignment, and records the
#' held-out accuracy and AUC at `lambda.min`. Returns the model with the best
#' AUC and the full performance grid.
#'
#' [metabodecon::benchmark_mdm()] wraps [metabodecon::cv_mdm()] in an
#' outer cross-validation loop to estimate end-to-end predictive performance.
#' It returns the per-fold models and held-out predictions.
#'
#' @details
#'
#' ## Disk caching across processes
#'
#' [metabodecon::cv_mdm()] and [metabodecon::benchmark_mdm()] can share
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
#' avoids redundant preprocessing when [metabodecon::cv_mdm()] evaluates
#' several settings on the same spectra.
#'
#' ## Parallelism
#'
#' [metabodecon::benchmark_mdm()] runs outer folds sequentially and
#' delegates all parallel work to [metabodecon::deconvolute()] and
#' [metabodecon::align()] via the `nworkers` argument.
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
#' @param maxCombine
#' Maximum peak-combining distance in datapoints, passed to
#' [metabodecon::get_si_mat()]. During training, no reference is passed so
#' partly-filled neighbouring columns are merged. During prediction, the stored
#' model feature positions are passed as `ref`, so new peaks are snapped to
#' those positions.
#' @param nworkers
#' Number of workers used by [metabodecon::fit_mdm()] and
#' [metabodecon::cv_mdm()] for deconvolution and alignment.
#' @param verbosity Integer. Verbosity level; each nested call
#' decrements by 1. Messages print when `verbosity >= 1`.
#' @param nfolds
#' Number of folds for the `cv.glmnet()` call in [metabodecon::fit_mdm()]
#' and [metabodecon::cv_mdm()]. Default 10.
#' @param cadir
#' Directory used to cache deconvolution results across grid points and
#' processes. Defaults to [metabodecon::decon_cachedir()].
#' @param nfo Number of outer folds in [metabodecon::benchmark_mdm()].
#' @param nfl
#' Number of folds for the `cv.glmnet()` call inside
#' [metabodecon::benchmark_mdm()]. Passed as `nfolds` to
#' [metabodecon::fit_mdm()].
#' @param seed
#' Random seed used for fold assignment in [metabodecon::benchmark_mdm()] and
#' for inner cross-validation splits in [metabodecon::cv_mdm()] and
#' [metabodecon::fit_mdm()].
#' @param check Logical. Whether to validate inputs at function entry.
#' @param pgrid
#' Data frame of preprocessing parameter combinations as returned
#' by [metabodecon::get_pgrid()].
#' @param ignore Optional integer vector of sample indices to exclude.
#' @param warm_cache
#' Logical. Whether to pre-populate the disk cache before
#' starting the grid search.
#' @param conf
#' Character string selecting a predefined parameter grid
#' configuration.
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with elements
#' `model` (a [glmnet::cv.glmnet()] object), `ref` (the reference alignment
#' spectrum) and `meta` (list of preprocessing parameters).
#'
#' [metabodecon::cv_mdm()] returns an `mdm` object with an additional element
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
#'          maxShift=128, maxCombine=256,
#'          nworkers=half_cores(), cadir=cadir
#'      )
#'
#'      mdm <- fit_mdm(
#'          spectra, y, sfr = NULL,
#'          nfit=5, smit=0, smws=0, delta=0, npmax=1000,
#'          maxShift=64, maxCombine=16,
#'          nworkers=half_cores(), cadir=cadir
#'      )
#'
#'
#'
#'      # -~-~-~-~-~-~-~-~-~ Search -~-~-~-~-~-~-~-~-~
#'
#'      mdm_grid_stat1 <- cv_mdm(
#'          spectra, y, pgrid=get_pgrid("static1"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'      saveRDS(mdm_grid_stat1, "tmp/mdm_grid_stat1.rds")
#'
#'
#'      mdm_grid_stat2 <- cv_mdm(
#'          spectra, y, pgrid=get_pgrid("static2"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'      saveRDS(mdm_grid_stat2, "tmp/mdm_grid_stat2.rds")
#'      #
#'      # Best static2: acc=0.82, auc=0.89
#'      # smit=2, smws=3, delta=6, nfit=7, npmax=0, maxShift=100, maxCombine=30
#'
#'
#'      mdm_grid_stat3 <- cv_mdm(
#'          spectra, y, pgrid=get_pgrid("static3"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'      saveRDS(mdm_grid_stat3, "tmp/mdm_grid_stat3.rds")
#'      #
#'      # Best static3: acc=79.34% auc=0.8589
#'      # smit=2, smws=7, delta=8, nfit=10, npmax=0, maxShift=100, maxCombine=30
#'
#'
#'      mdm_grid_dyn2 <- cv_mdm(
#'          spectra, y, pgrid=get_pgrid("dynamic2"),
#'          nworkers=half_cores(), use_rust=TRUE, cadir=cadir
#'      )
#'
#'      # -~-~-~ Full Benchmark -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
#'      bm <- benchmark_mdm(spectra, y, sfr, cadir=cadir)
#'
#'      # -~-~-~ Interactive Development -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
#'      stub(fit_mdm, spectra=spectra, y=y, sfr=sfr)
#'      stub(cv_mdm, spectra=spectra, y=y, sfr=sfr)
#'      stub(benchmark_mdm, spectra=spectra, y=y, sfr=sfr)
#' }
#'
fit_mdm <- function(
    # Mandatory
    spectra, y,
    # Shared
    sfr = NULL, use_rust = 0.5, nworkers = 1, verbosity = 2,
    seed = 1, cadir = decon_cachedir(), check = TRUE,
    # Deconvolution
    npmax = 1000, nfit = 5, smit = 2, smws = 5, delta = 6.4,
    # Alignment
    maxShift = 200,
    # Peak combining
    maxCombine = 50,
    # Glmnet
    nfolds = 10
) {

    if (check) check_mdm_args(
        spectra=spectra, y=y, sfr=sfr, use_rust=use_rust, nworkers=nworkers,
        verbosity=verbosity, seed=seed, cadir=cadir, check=check, nfolds=nfolds,
        npmax=npmax, nfit=nfit, smit=smit, smws=smws, delta=delta,
        maxShift=maxShift, maxCombine=maxCombine
    )

    logv("Deconvoluting (npmax=%d, nfit=%d, smit=%d, smws=%d, delta=%.1f)",
         npmax, nfit, smit, smws, delta)
    shash <- attr(spectra, "hash")
    if (is.null(shash)) stop("spectra must carry a 'hash' attribute.", call.=FALSE)
    decon_hash <- rlang::hash(list(shash, sfr, use_rust, npmax, nfit, smit, smws, delta))
    decons <- deconvolute_spectra(
        x=spectra, sfr=sfr, verbose=verbosity >= 2, use_rust=use_rust,
        nfit=nfit, smit=smit, smws=smws, delta=delta, npmax=npmax,
        nworkers=nworkers, cadir=cadir, hash=decon_hash, full=FALSE
    )
    nps <- unname(unlist(lapply(decons, function(d) nrow(d$lcpar))))
    npzero <- sum(nps == 0)
    if (npzero > 0) {
        logv("Aborting: %d/%d spectra produced zero peaks", npzero, length(nps))
        meta <- list(sfr=sfr, use_rust=use_rust, npmax=npmax, nfit=nfit,
            smit=smit, smws=smws, delta=delta, maxShift=maxShift, maxCombine=maxCombine)
        return(structure(list(model=NULL, ref=NULL, meta=meta), class="mdm"))
    }
    logv("Peak range: [%d, %d]", range(nps)[1], range(nps)[2])
    attr(decons, "hash") <- decon_hash

    logv("Aligning (maxShift=%d, maxCombine=%g)", maxShift, maxCombine)
    aligns <- align_decons(
        xx=decons, maxShift=maxShift, verbose=verbosity >= 2,
        nworkers=nworkers, full=FALSE
    )

    logv("Constructing feature matrix")
    ref <- find_ref(aligns)
    mat <- si_mat(aligns, maxCombine=maxCombine)
    peak_cols <- which(colSums(mat != 0) > 0)
    intervals <- make_intervals(peak_cols, maxCombine, nc=ncol(mat))
    X <- mat[, peak_cols, drop=FALSE]

    logv("Fitting cv.glmnet")
    foldid <- get_foldid(y=y, nfolds=nfolds, seed=seed)
    model <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1, foldid=foldid, keep=TRUE)

    logv("Constructing and returning mdm object")
    meta <- list(sfr=sfr, use_rust=use_rust, npmax=npmax, nfit=nfit, smit=smit,
        smws=smws, delta=delta, maxShift=maxShift, maxCombine=maxCombine, intervals=intervals)
    structure(list(model=model, ref=ref, meta=meta), class="mdm")
}

#' @export
#' @rdname mdm
cv_mdm <- function(
    # Mandatory
    spectra, y, sfr = NULL,
    # Shared
    use_rust = 0.5, nworkers = 1, verbosity = 2,
    seed = 1, cadir = decon_cachedir(), check = TRUE,
    # Grid
    pgrid = get_pgrid("dynamic2"),
    ignore = NULL, warm_cache = TRUE,
    # Glmnet
    nfolds = 10
) {

    if (check) check_mdm_args(
        spectra=spectra, y=y, sfr=sfr, use_rust=use_rust, nworkers=nworkers,
        verbosity=verbosity, nfolds=nfolds, seed=seed, cadir=cadir,
        pgrid=pgrid, ignore=ignore, check=check
    )

    if (!is.null(ignore)) {
        logv("Removing %d ignored samples from spectra and y", length(ignore))
        spectra <- spectra[-ignore]
        y <- y[-ignore]
    }

    if (warm_cache) init_cache(
        x=spectra, sfr=sfr, verbosity=verbosity - 1,
        nworkers=nworkers, use_rust=use_rust, cadir=cadir, pgrid=pgrid
    )

    on.exit(options(metabodecon.align_cache=NULL, metabodecon.decon_cache=NULL), add=TRUE)

    np <- nrow(pgrid)
    ns <- length(spectra)
    logv("Starting grid search (%d combinations, %s spectra)", np, ns)

    foldid <- get_foldid(y=y, nfolds=nfolds, seed=seed)
    pp_cols <- intersect(
        c("npmax", "nfit", "smit", "smws", "delta", "maxShift", "maxCombine"),
        names(pgrid)
    )

    best_mdm <- NULL
    best_auc <- -Inf
    for (i in seq_len(np)) {
        pc <- as.list(pgrid[i, pp_cols, drop=FALSE])
        tag <- paste(sprintf("%s=%g", pp_cols, unlist(pc)), collapse=", ")
        logv("[%d/%d]: %s", i, np, tag)
        mdm <- do.call(fit_mdm, c(
            list(spectra=spectra, y=y, sfr=sfr, use_rust=use_rust, nworkers=nworkers,
                 verbosity=verbosity - 1, seed=seed, cadir=cadir, check=FALSE, nfolds=nfolds),
            pc
        ))
        if (is.null(mdm$model)) {
            pgrid$acc[i] <- 0; pgrid$auc[i] <- 0
            logv("[%d/%d] zero peaks, skipping", i, np)
            next
        }
        cvfit <- mdm$model
        li <- which(cvfit$lambda == cvfit$lambda.min)
        link <- cvfit$fit.preval[, li]
        prob <- 1 / (1 + exp(-link))
        lvs <- levels(y)
        pred <- factor(ifelse(prob > 0.5, lvs[2], lvs[1]), levels=lvs)
        pgrid$acc[i] <- mean(pred == y)
        pgrid$auc[i] <- AUC(y, prob)
        if (pgrid$auc[i] > best_auc) { best_auc <- pgrid$auc[i]; best_mdm <- mdm }
        logv("[%d/%d] acc=%.2f%% auc=%.4f | best auc=%.4f",
            i, np, pgrid$acc[i] * 100, pgrid$auc[i], best_auc)
    }
    ibest <- which.max(pgrid$auc)
    logv("Best [%d/%d]: acc=%.2f%% auc=%.4f", ibest, np,
        pgrid$acc[ibest] * 100, pgrid$auc[ibest])
    best_mdm$pgrid <- pgrid
    best_mdm
}

#' @export
#' @rdname mdm
benchmark_mdm <- function(
    spectra, y, sfr=NULL,
    use_rust=0.5, verbosity=2,  seed=1, cadir=decon_cachedir(),
    nworkers=half_cores(), pgrid=get_pgrid(), nfo=5, nfl=10
) {

    check_mdm_args(
        spectra=spectra, y=y, sfr=sfr, use_rust=use_rust, verbosity=verbosity,
        nfolds=nfo, nfl=nfl, nworkers=nworkers, cadir=cadir, pgrid=pgrid,
        seed=seed
    )

    init_cache(
        x=spectra, sfr=sfr, verbosity=verbosity-1, nworkers=nworkers,
        use_rust=use_rust, cadir=cadir, pgrid=pgrid
    )

    logv("Fitting cv_mdm for %d outer folds", nfo)
    te_list <- get_test_ids(nfolds=nfo, nsamples=length(spectra), seed=seed, y=y)
    fids <- seq_along(te_list)
    mdms <- lapply(te_list, function(ignore) cv_mdm(
        ignore=ignore, spectra=spectra, y=y, sfr=sfr, use_rust=use_rust,
        nworkers=nworkers, verbosity=verbosity - 1, nfolds=nfl, seed=seed,
        cadir=cadir, pgrid=pgrid, warm_cache=FALSE, check=FALSE
    ))

    logv("Predicting on held-out test folds")
    te_spectra <- lapply(te_list, function(te) spectra[te])
    preds_raw <- mapply(
        predict.mdm, object=mdms, newdata=te_spectra,
        MoreArgs=list(type="all", nworkers=nworkers, verbosity=verbosity - 1),
        SIMPLIFY=FALSE, USE.NAMES=FALSE
    )
    preds <- mapply(
        FUN = function(fid, te, p) {
            data.frame(fold=fid, true=y[te], link=p$link, prob=p$prob, pred=p$class)
        },
        fids, te_list, preds_raw,
        SIMPLIFY=FALSE, USE.NAMES=FALSE
    )
    Y <- do.call(rbind, preds)
    acc <- mean(Y$true == Y$pred, na.rm=TRUE)
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
            maxCombine = seq(10, 50, 10),
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
            maxCombine = seq(10, 50, 10),
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
            maxCombine = seq(10, 50, 10),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
    } else {
        stop("Invalid conf")
    }
    ord <- order(
        P$npmax, P$nfit, P$smit, P$smws, P$delta,
        P$maxShift, P$maxCombine
    )
    P <- P[ord, ]
    P$acc <- NA_real_
    P$auc <- NA_real_
    rownames(P) <- NULL
    P
}

# Helpers #####

init_cache <- function(
    x, sfr, pgrid = get_pgrid(), nworkers = 1,
    verbosity = 1, use_rust = 0.5, cadir = decon_cachedir()
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
            x = x, sfr = sfr, verbose = verbosity >= 2,
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
            deconvolute_spectra(x=xu, sfr=sfr, verbose=verbosity>1, nworkers=nwu_i, use_rust=use_rust, npmax=npmax
            )
        }
    }

    logv("Finished cache initialization")
    invisible(NULL)
}

as_binary01 <- function(y) {
    lvs <- sort(unique(y))
    if (length(lvs) != 2) stop("y must have exactly 2 unique levels")
    as.integer(y == lvs[2])
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
    nfolds = NULL, nfl = NULL,
    cadir = NULL, pgrid = NULL,
    ignore = NULL, seed = NULL, check = NULL,
    npmax = NULL, nfit = NULL, smit = NULL, smws = NULL, delta = NULL,
    maxShift = NULL, maxCombine = NULL
) {
    pp_names <- c(
        "npmax", "nfit", "smit", "smws", "delta",
        "maxShift", "maxCombine"
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
        is_bool_or_num(use_rust),
        is_int_or_null(nworkers, 1),
        is_int_or_null(verbosity, 1),
        is.null(nfolds) || (is_int(nfolds, 1) && nfolds >= 2),
        is.null(nfl) || (is_int(nfl, 1) && nfl >= 2),
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
#'                 delta = 6.4, maxShift = 100, maxCombine = 50)
#'   ),
#'   class = "mdm"
#' )
#' print(m)
#' summary(m)
#'
#' \dontrun{
#'   m <- cv_mdm(spectra, y, sfr = c(11, -2))
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
        logv("Deconvoluting %d spectra with %d nworkers", length(newdata), nworkers)
        decons <- if (m$npmax > 0) {
            deconvolute(x=newdata, sfr=m$sfr, verbose=verbosity >= 2,
                use_rust=m$use_rust, npmax=m$npmax, nworkers=nworkers)
        } else {
            deconvolute(x=newdata, sfr=m$sfr, verbose=verbosity >= 2,
                use_rust=m$use_rust, nfit=m$nfit, smit=m$smit, smws=m$smws,
                delta=m$delta, npmax=0, nworkers=nworkers)
        }
        logv("Aligning spectra with %d nworkers", nworkers)
        als <- align_decons(xx=decons, maxShift=m$maxShift, verbose=verbosity >= 2,
            nworkers=nworkers, ref=object$ref, full=FALSE)
        Xn <- si_mat(als, maxCombine=m$maxCombine, intervals=m$intervals)
    } else {
        Xn <- as.matrix(newdata)
    }
    logv("Predicting with s=%s", as.character(s))
    requireNamespace("glmnet", quietly = TRUE)
    score <- as.numeric(predict(object$model, newx = Xn, s = s, type = "link"))
    prob <- as.numeric(predict(object$model, newx = Xn, s = s, type = "response"))
    pred <- predict(object$model, newx = Xn, s = s, type = "class")[, 1]
    lvs <- object$model$glmnet.fit$classnames
    pred <- factor(pred, levels = lvs)
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
            "maxShift", "maxCombine")
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
        "maxShift", "maxCombine")
    out <- object$meta[pp]
    out$n_peaks <- if (is.null(object$meta$intervals)) 0L
                   else nrow(object$meta$intervals)
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
