
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
#' Fit ([metabodecon::fit_mdm()]) or cross-validate
#' ([metabodecon::benchmark()]) a binary classification model built on a
#' set of NMR spectra. Both run the full deconvolute -> align -> snap ->
#' featurize -> fit pipeline with sensible defaults and expose only the
#' parameters a typical user tunes; the classification backend is chosen
#' via `model`. Power users who need to swap individual pipeline stages
#' can call the internal engines `metabodecon:::fit_mdm_internal()` /
#' `metabodecon:::benchmark_internal()`, which take pluggable `decon_fun`
#' / `align_fun` / `snap_fun` / `feat_fun` / `fit_fun` / `predict_fun`
#' arguments.
#'
#' [metabodecon::fit_mdm()] runs the pipeline once, or iterates over the
#' cartesian product of `npmax` / `maxShift` / `maxCombine` when any is a
#' vector and returns the row with the highest `acc` (ties broken by
#' `auc`). [metabodecon::benchmark()] wraps [metabodecon::fit_mdm()] in
#' outer k-fold cross-validation to estimate end-to-end performance on
#' held-out spectra.
#'
#' @param x Spectra object.
#' @param y Factor vector with class labels for each spectrum.
#' @param model Classification backend. One of `"lasso"` (default,
#'   L1-penalised logistic regression via `glmnet`) or `"ranger"`
#'   (probability random forest).
#' @param npmax Max peaks per spectrum. Integer in `{-2, -1, 0, 1, ...}`,
#'   scalar or vector. `-1` (default) selects the median per-spectrum
#'   Kneedle elbow; `-2` selects each spectrum's own elbow.
#' @param maxShift Max CluPA shift in datapoints. Integer >= -1, scalar or
#'   vector. `-1` (default) means auto (sweep to the alignment-correlation
#'   dip).
#' @param maxCombine RefPA snap window in datapoints. Integer, scalar or
#'   vector. Default 10.
#' @param nworkers Number of workers for deconvolution, alignment and the
#'   inner fitter.
#' @param seed Random seed. Forwarded to the fitter; also used for
#'   stratified fold assignment inside [metabodecon::benchmark()]. May be a
#'   vector for repeated CV.
#' @param verbosity Verbosity level.
#' @param k Number of outer folds for [metabodecon::benchmark()].
#' @param ... Further arguments passed on to the internal engine
#'   (`metabodecon:::fit_mdm_internal()` /
#'   `metabodecon:::benchmark_internal()`), e.g. `sfr`, `igrs`, `deg`,
#'   `use_rust`. Rarely needed.
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with elements
#' `model` (trained backend model of the best grid row), `ref` (a list
#' `list(align, snap)` for prediction-time replay), `params` (resolved
#' pipeline parameters of the best row), the scalar performance of the best
#' row (`acc`, `auc`, `acc_se`, `auc_se`), and `mog` (the augmented grid
#' with per-row performance).
#'
#' [metabodecon::benchmark()] returns a list with elements `models` (one
#' fitted model per outer fold), `predictions` (per-spectrum out-of-fold
#' predictions), `performance` (per-fold `acc` / `auc`) and `overall`
#' (pooled `acc` / `auc`).
#'
#' @examples
#' \dontrun{
#'   x <- sim2
#'   y <- attr(sim2, "group")
#'   m  <- fit_mdm(x, y)                     # lasso, full pipeline
#'   mr <- fit_mdm(x, y, model="ranger")     # random forest
#'   bm <- benchmark(x, y, k=5)              # 5-fold CV
#'   fm <- fit_mdm(
#'       x, y, model="ranger",
#'       npmax=25L, maxShift=c(0L, 1L, 2L, 4L, 8L),
#'       maxCombine=c(0L, 1L, 2L, 4L), nworkers=4L
#'   )
#' }
#'
fit_mdm <- function(
    x, y, model=c("lasso", "ranger"),
    npmax=-1L, maxShift=-1L, maxCombine=10L,
    nworkers=1L, seed=1L, verbosity=1L, ...
) {
    model <- match.arg(model)
    fit_fun <- if (model == "ranger") fit_ranger else fit_lasso
    predict_fun <- if (model == "ranger") predict_ranger else predict_lasso
    fit_mdm_internal(
        x=x, y=y, fit_fun=fit_fun, predict_fun=predict_fun,
        npmax=npmax, maxShift=maxShift, maxCombine=maxCombine,
        nworkers=nworkers, seed=seed, verbosity=verbosity, ...
    )
}

#' @export
#' @rdname mdm
benchmark <- function(
    x, y, model=c("lasso", "ranger"),
    npmax=-1L, maxShift=-1L, maxCombine=10L,
    nworkers=1L, seed=1L, verbosity=2L, k=3L, ...
) {
    model <- match.arg(model)
    fit_fun <- if (model == "ranger") fit_ranger else fit_lasso
    predict_fun <- if (model == "ranger") predict_ranger else predict_lasso
    benchmark_internal(
        x=x, y=y, fit_fun=fit_fun, predict_fun=predict_fun,
        npmax=npmax, maxShift=maxShift, maxCombine=maxCombine,
        nworkers=nworkers, seed=seed, verbosity=verbosity, k=k, ...
    )
}

# Engine (private) #####

#' @noRd
#'
#' @title Metabodecon Models (internal pluggable engine)
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' **WARNING: These functions are experimental and must not be used in
#' production. Their API is very likely to change in non-backwards-compatible
#' ways over the next few weeks.**
#'
#' Fit and benchmark binary classification models built on NMR spectra
#' through a pluggable five-stage pipeline (deconvolute, align, snap,
#' featurize, fit). [metabodecon::fit_mdm()] runs the pipeline once, or
#' iterates over the cartesian product of `npmax`/`maxShift`/`maxCombine`
#' when any is a vector and returns the row with the highest `acc`
#' (ties broken by `auc`).
#' [metabodecon::benchmark()] wraps [metabodecon::fit_mdm()] in outer
#' k-fold cross-validation to estimate end-to-end performance on
#' held-out spectra.
#'
#' @param x Spectra object.
#' @param y Factor vector with class labels for each spectrum.
#' @param decon_fun Function
#'   `(x, sfr, igrs, verbose, use_rust, npmax, nworkers) -> decons2`.
#'   Built-ins: [metabodecon::deconvolute()] (default),
#'   [metabodecon::identity2()]. The underlying
#'   `(nfit, smit, smws, delta)` tuple is no longer a `fit_mdm()`
#'   argument — it is picked from each spectrum's `$deg` cache via
#'   `npmax`.
#' @param align_fun Function
#'   `(x, ref, maxShift, verbose, nworkers, full, ...) -> aligns`
#'   (or pass-through). Built-ins: [metabodecon::clupa()] (default),
#'   [metabodecon::identity_align()].
#' @param snap_fun Function `(x, ref=NULL, maxCombine, ...) -> aligns`
#'   with the per-peak `pcisn` / `x0sn` columns populated. Built-ins:
#'   [metabodecon::snap_to_ref()] (default),
#'   [metabodecon::combine_peaks()], [metabodecon::snap_nw_blind()],
#'   [metabodecon::identity_snap()]. Pass `identity2` / `identity_align`
#'   / `identity_snap` to skip a stage and fit baselines (e.g. a binning
#'   model) directly on raw spectra.
#' @param feat_fun Function
#'   `(x, maxCombine, igrs, peakPos=NULL, ...) -> matrix` with one row
#'   per spectrum. With `peakPos=NULL` (training), filters all-zero
#'   columns and attaches the kept indices on `attr(., "peakPos")`;
#'   with `peakPos=<int>` (predict), subsets to those columns.
#'   Built-ins: [metabodecon::peak_mat()] (default),
#'   [metabodecon::si_mat()], [metabodecon::bin()].
#' @param fit_fun Function
#'   `(X, y, seed, nworkers) -> list(model, acc, auc, acc_se, auc_se)`.
#'   The acc/auc scalars are generalization estimates in `[0, 1]`; use
#'   `NA_real_` for SEs when the backend produces a single point estimate.
#'   Must call `requireNamespace("<backend>")` so it works after a fresh
#'   `readRDS()`. Built-ins: [metabodecon::fit_lasso()] (default,
#'   repeated `cv.glmnet` OOF), [metabodecon::fit_ranger()] (OOB).
#' @param predict_fun Function `(model, newx) -> numeric` of
#'   positive-class probabilities. Must call
#'   `requireNamespace("<backend>")`. Built-ins:
#'   [metabodecon::predict_lasso()], [metabodecon::predict_ranger()].
#' @param npmax Max peaks per spectrum. Integer in `{-2, -1, 0, 1,
#'   ...}`, scalar or vector. Drives selection of the underlying
#'   `(nfit, smit, smws, delta)` row from each spectrum's `$deg`
#'   cache. See *Automatic selection sentinels*.
#' @param maxShift Max CluPA shift in datapoints. Integer >= -1, scalar
#'   or vector. `-1` means auto. Default `-1`.
#' @param maxCombine RefPA snap window in datapoints. Integer, scalar
#'   or vector. Default 10.
#' @param deg Deconvolution-parameter grid forwarded to
#'   [metabodecon::grid_deconvolute_spectra()]. When `NULL` (default),
#'   the default grid built into
#'   [metabodecon::grid_deconvolute_spectra()] is used.
#' @param sfr Signal-free region. See [metabodecon::deconvolute()].
#' @param igrs Ignore regions in ppm.
#' @param use_rust Use the Rust backend?
#' @param nworkers Number of workers for deconvolution, alignment and
#'   the inner fitter.
#' @param verbosity Verbosity level.
#' @param seed Random seed. Forwarded to `fit_fun`; also used for
#'   stratified fold assignment inside [metabodecon::benchmark()].
#' @param k Number of outer folds for [metabodecon::benchmark()].
#' @param check Validate inputs at function entry?
#'
#' @details
#'
#' ## Pipeline
#'
#' ```
#' d <- decon_fun(x, npmax, ...)
#' a <- align_fun(d, maxShift, ...)
#' s <- snap_fun(a, maxCombine, ...)
#' X <- feat_fun(s, maxCombine, ...)
#' m <- fit_fun(X, y, ...)
#' ```
#'
#' ## Automatic selection sentinels
#'
#' - `npmax = -1` (default) — median per-spectrum Kneedle elbow on `$deg`.
#' - `npmax = -2` — each spectrum's own Kneedle elbow ("intrinsic");
#'   parameter-free at the cohort level.
#' - `maxShift = -1` (default) — sweep CluPA shifts at powers of 2 and
#'   stop one step before the alignment-correlation dip. Requires a
#'   CluPA-compatible `align_fun` (must populate `sit$supal`).
#'
#' ## Caching and parallelism
#'
#' Grid rows are sorted by `(npmax, maxShift, maxCombine)` so the most
#' recent deconvolution and alignment are reused across rows.
#' [metabodecon::fit_mdm()] runs an idempotent
#' [metabodecon::grid_deconvolute_spectra()] up front to attach `$deg`
#' to each spectrum (pre-enriched spectra skip the slow step).
#' [metabodecon::benchmark()] runs outer folds sequentially; all
#' parallelism is delegated to `fit_fun` via `nworkers`.
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with
#' elements `model` (trained backend model of the best grid row),
#' `ref` (a list `list(align, snap)` for prediction-time replay —
#' both elements are `attr(., "ref")` of the corresponding stage
#' output), `params` (resolved scalar pipeline parameters of the best
#' row plus the pluggable functions, `lvs`, `peakPos`, `sfr`, `igrs`,
#' `use_rust`), the scalar performance of the best row (`acc`, `auc`,
#' `acc_se`, `auc_se`), and `mog` (the augmented grid with per-row
#' `acc`, `auc`, `acc_se`, `auc_se`).
#'
#' [metabodecon::benchmark()] returns a list with elements:
#'
#' - `models`: list of fitted models, one per outer fold.
#' - `predictions`: data frame with columns `fold`, `true`, `link`,
#'   `prob`, `pred`.
#' - `performance`: data frame with per-fold `acc` and `auc`.
#' - `overall`: list with pooled `acc` and `auc`.
#'
#' @examples
#' \dontrun{
#'   # Original examples on a private AKI dataset.
#'   aki <- read_aki_data()
#'   x <- aki$spectra
#'   y <- aki$meta$type; names(y) <- aki$meta$sid
#'   m <- fit_mdm(x, y)
#'   bm <- benchmark(x, y, k=5)
#'   mrf <- fit_mdm(x, y, fit_fun=fit_ranger, predict_fun=predict_ranger)
#'   mnw <- fit_mdm(
#'       x, y, snap_fun=snap_nw_blind,
#'       fit_fun=fit_ranger, predict_fun=predict_ranger
#'   )
#'
#'   # The 3.4-real-perf benchmarks reproduced on the public sim2
#'   # dataset. Still slow (deconvolution + grid search across outer
#'   # folds), but much quicker than the AKI version above.
#'   x <- sim2
#'   y <- attr(sim2, "group")
#'
#'   # (a) Binning baseline: ranger on raw spectra binned to fixed
#'   # ppm columns (no deconvolution, no alignment, no snap). The
#'   # `maxCombine` sweep selects the bin width (in datapoints).
#'   bm_bin <- benchmark(
#'       x, y,
#'       decon_fun=identity2, align_fun=identity_align,
#'       snap_fun=identity_snap, feat_fun=bin,
#'       fit_fun=fit_ranger, predict_fun=predict_ranger,
#'       maxCombine=c(5L, 10L, 20L, 40L), k=5L, nworkers=4L
#'   )
#'
#'   # (b) Full metabodecon pipeline (deconvolute -> CluPA -> RefPA
#'   # -> peak_mat -> ranger), sweeping (maxShift, maxCombine).
#'   bm_mdm <- benchmark(
#'       x, y,
#'       fit_fun=fit_ranger, predict_fun=predict_ranger,
#'       npmax=25L, maxShift=c(0L, 1L, 2L, 4L, 8L, 16L),
#'       maxCombine=c(0L, 1L, 2L, 4L, 8L), k=5L, nworkers=4L
#'   )
#'
#'   # (c) Final model trained on the full sim2 dataset across the
#'   # same grid; `$mog` carries the per-row (acc, auc) heatmap and
#'   # the ACC-best row (ties broken by AUC) is auto-selected.
#'   fm <- fit_mdm(
#'       x, y,
#'       fit_fun=fit_ranger, predict_fun=predict_ranger,
#'       npmax=25L, maxShift=c(0L, 1L, 2L, 4L, 8L, 16L),
#'       maxCombine=c(0L, 1L, 2L, 4L, 8L), nworkers=4L
#'   )
#' }
#'
fit_mdm_internal <- function(
    x, y,
    decon_fun=deconvolute_spectra, align_fun=clupa, snap_fun=snap_to_ref,
    feat_fun=peak_mat, fit_fun=fit_lasso, predict_fun=predict_lasso,
    npmax=-1L, maxShift=-1L, maxCombine=10L,
    deg=NULL, sfr=NULL, igrs=list(), use_rust=0, nworkers=1,
    verbosity=1, seed=1, check=TRUE
) {

    lvs <- levels(y)
    verbose <- verbosity >= 2L
    x <- grid_deconvolute_spectra(x, deg, sfr, igrs, verbose, nworkers, use_rust)
    npmax[npmax == -1L] <- find_npmax_elbow(x)
    g <- get_mog(npmax, maxShift, maxCombine)
    nr <- nrow(g)
    ns <- length(x)
    last_np <- -99L
    last_ms <- -99L;
    d <- NULL
    a <- NULL
    # Grid winner is picked by accuracy with AUC as tiebreaker (matches
    # the rest of the pipeline; ranger's probability output sometimes
    # squeezes toward 0.5 on small samples, making AUC-only picks
    # land on cells with mediocre headline accuracy).
    best_acc <- -Inf
    best_auc <- -Inf
    best_i <- NA_integer_
    best_model <- NULL
    best_pp <- NULL
    best_refs <- NULL

    logv("Grid search (%d rows, %d spectra)", nr, ns)
    perf_cols <- c("acc", "auc", "acc_se", "auc_se")
    for (i in seq_len(nr)) {
        np <- g$npmax[i]
        ms <- g$maxShift[i]
        mc <- g$maxCombine[i]
        if (np != last_np) {
            d <- decon_fun(
                x=x, sfr=sfr, igrs=igrs, verbose=verbose,
                use_rust=use_rust, npmax=np, nworkers=nworkers
            )
            last_np <- np
            last_ms <- -99L
        }
        if (ms == -1L) {
            ms <- find_maxShift_dip(d, align_fun=align_fun, nworkers=nworkers, verbose=verbose)
            g$maxShift[g$npmax == np & g$maxShift == -1L] <- ms
        }
        if (ms != last_ms) {
            a <- align_fun(
                x=d, ref=NULL, maxShift=ms, verbose=verbose,
                nworkers=nworkers, full=FALSE
            )
            last_ms <- ms
        }
        s <- snap_fun(a, ref=NULL, maxCombine=mc, igrs=igrs)
        X <- feat_fun(s, maxCombine=mc, igrs=igrs, peakPos=NULL)
        r <- fit_fun(X, y, seed=seed, nworkers=nworkers)
        g[i, perf_cols] <- r[perf_cols]

        is_best <- !is.na(r$acc) && !is.na(r$auc) && (
            r$acc > best_acc || (r$acc == best_acc && r$auc > best_auc)
        )
        sym <- if (is_best) " <-- BEST" else ""
        fmt <- "[%d/%d] np=%d S=%d C=%d acc=%s auc=%s%s"
        logv(fmt, i, nr, np, ms, mc, r$acc, r$auc, sym)
        if (is_best) {
            best_acc <- r$acc
            best_auc <- r$auc
            best_i <- i
            best_model <- r$model
            best_pp <- attr(X, "peakPos")
            best_refs <- list(align=attr(a, "ref"), snap=attr(s, "ref"))
        }
    }
    if (is.na(best_i)) stop("No grid row produced a model.", call.=FALSE)
    params <- list(
        feat_fun=feat_fun, fit_fun=fit_fun, predict_fun=predict_fun,
        decon_fun=decon_fun, align_fun=align_fun, snap_fun=snap_fun,
        lvs=lvs, sfr=sfr, igrs=igrs, use_rust=use_rust,
        npmax=g$npmax[best_i], maxShift=g$maxShift[best_i],
        maxCombine=g$maxCombine[best_i], peakPos=best_pp
    )
    ret <- list(
        model=best_model, ref=best_refs, params=params, acc=g$acc[best_i],
        auc=g$auc[best_i], acc_se=g$acc_se[best_i], auc_se=g$auc_se[best_i], mog=g
    )
    structure(ret, class="mdm")
}

#' @noRd
benchmark_internal <- function(
    x, y,
    decon_fun=deconvolute,  align_fun=clupa,    snap_fun=snap_to_ref,
    feat_fun=peak_mat,      fit_fun=fit_lasso,  predict_fun=predict_lasso,
    npmax=-1L,              maxShift=-1L,       maxCombine=10L,
    deg=NULL,               sfr=NULL,           igrs=list(),
    use_rust=0,             nworkers=1,         verbosity=2,
    seed=1,                 k=3,                check=TRUE
) {
    stopifnot(
        is.function(decon_fun), is.function(align_fun),
        is.function(snap_fun),  is.function(feat_fun),
        is.function(fit_fun),   is.function(predict_fun),
        is_int(k, 1), k >= 2, k <= length(y),
        is_int(npmax),    all(npmax    >= -2L),
        is_int(maxShift), all(maxShift >= -1L),
        is_int(maxCombine), all(maxCombine >= 0L)
    )
    if (check) check_mdm_args(
        x=x, y=y, sfr=sfr, igrs=igrs, use_rust=use_rust,
        nworkers=nworkers, verbosity=verbosity, seed=seed
    )

    # One-time grid attach so each per-fold fit_mdm() sees pre-enriched
    # spectra and skips this step.
    if (!identical(decon_fun, identity2)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2, nworkers=nworkers, use_rust=use_rust
        )
    }

    # Inner fit_mdm / predict run silently by default — the per-fold result
    # line below is the only log emission per fold. Bumping outer verbosity
    # to 3 lets the user re-enable the fit_mdm grid-search output for
    # debugging.
    inner_v <- max(0L, verbosity - 2L)
    # `te_list` has length k * length(seed). Each test-id vector carries
    # `seed` / `fold` attributes (see get_test_ids). For a scalar `seed`
    # the attributes are unset and we fall back to (seed, fold=i).
    te_list <- get_test_ids(nfolds=k, nsamples=length(x), seed=seed, y=y)
    nf <- length(te_list)
    models <- vector("list", nf)
    fold_preds <- vector("list", nf)
    perf <- data.frame(seed=integer(0), fold=integer(0),
                       acc=numeric(0), auc=numeric(0))
    if (length(seed) > 1L) {
        logv("Running %d-fold outer CV x %d seeds with fit_mdm",
             k, length(seed))
    } else {
        logv("Running %d-fold outer CV with fit_mdm", k)
    }
    for (i in seq_along(te_list)) {
        te <- te_list[[i]]
        s <- attr(te, "seed") %||% seed
        f <- attr(te, "fold") %||% i
        tr <- setdiff(seq_along(x), te)
        m <- fit_mdm_internal(
            x=x[tr], y=y[tr],
            decon_fun=decon_fun, align_fun=align_fun, snap_fun=snap_fun,
            feat_fun=feat_fun, fit_fun=fit_fun, predict_fun=predict_fun,
            npmax=npmax, maxShift=maxShift, maxCombine=maxCombine,
            deg=deg, sfr=sfr, igrs=igrs,
            use_rust=use_rust, nworkers=nworkers, verbosity=inner_v,
            seed=s, check=FALSE
        )
        p <- stats::predict(m, x[te], type="all", nworkers=nworkers, verbosity=inner_v)
        fp <- data.frame(seed=s, fold=f, true=y[te], link=p$link,
                         prob=p$prob, pred=p$class)
        fold_preds[[i]] <- fp
        acc <- mean(fp$pred == fp$true, na.rm=TRUE)
        auc <- AUC(fp$true, fp$prob)
        perf <- rbind(perf, data.frame(seed=s, fold=f, acc=acc, auc=auc))
        if (length(seed) > 1L) {
            fmt <- "[seed %d, fold %d/%d] acc=%.3f auc=%.3f | mean acc=%.3f auc=%.3f"
            logv(fmt, s, f, k, acc, auc, mean(perf$acc), mean(perf$auc))
        } else {
            fmt <- "[fold %d/%d] acc=%.3f auc=%.3f | mean acc=%.3f auc=%.3f"
            logv(fmt, f, k, acc, auc, mean(perf$acc), mean(perf$auc))
        }
        models[[i]] <- m
    }

    preds <- do.call(rbind, fold_preds)
    overall_acc <- mean(preds$true == preds$pred, na.rm=TRUE)
    overall_auc <- AUC(preds$true, preds$prob)
    logv("Overall: acc=%.3f auc=%.3f", overall_acc, overall_auc)
    # When seeds were swept, also expose per-seed mean ± sd so callers
    # can report repeated-CV stability without re-aggregating from `perf`.
    overall <- list(acc=overall_acc, auc=overall_auc)
    if (length(seed) > 1L) {
        per_seed <- stats::aggregate(
            cbind(acc, auc) ~ seed, data=perf, FUN=mean
        )
        overall$acc_seed_mean <- mean(per_seed$acc)
        overall$acc_seed_sd   <- stats::sd(per_seed$acc)
        overall$auc_seed_mean <- mean(per_seed$auc)
        overall$auc_seed_sd   <- stats::sd(per_seed$auc)
    }
    list(
        models=models,
        predictions=preds,
        performance=perf,
        overall=overall
    )
}

#' @noRd
#' @title Identity decon function for fit_mdm
#' @description
#' No-op replacement for the `decon_fun` argument of
#' [metabodecon::fit_mdm()] / [metabodecon::benchmark()]. Returns its
#' first argument unchanged and ignores all other arguments. Use this to
#' skip the deconvolution stage of the pipeline (e.g. for binning
#' baselines).
#' @return `x`, unchanged.
identity2 <- function(x, ...) x

#' @noRd
#' @title Identity snap function for fit_mdm
#' @description
#' No-op replacement for the `snap_fun` argument of
#' [metabodecon::fit_mdm()]. Returns its first argument unchanged so the
#' feature-matrix stage operates on whatever the alignment stage
#' produced (without an extra snap-to-reference step).
#' @return `x`, unchanged.
identity_snap <- function(x, ref=NULL, maxCombine=0L, ...) x

#' @noRd
#' @title Lasso fitter for fit_mdm
#' @description
#' Fits an L1-penalised binomial logistic regression via repeated
#' [glmnet::cv.glmnet()] (default `nreps=5`, internal `nfolds=10`,
#' `keep=TRUE`). All reps share the lambda path discovered by the first
#' rep so per-rep out-of-fold (OOF) predictions are directly comparable
#' lambda-by-lambda. For each lambda the per-rep OOF accuracy and AUC
#' are averaged across reps; the lambda that maximizes the averaged
#' accuracy (with averaged AUC as tiebreaker) is reported as
#' `lambda*` and used at predict time. Averaging the per-lambda
#' performance curve *before* optimizing over lambda keeps `lambda*`
#' stable across reps and the reported acc/AUC reflects model
#' variance rather than the noise of lambda-pick instability across
#' reps. `acc_se` / `auc_se` are the SE across reps at `lambda*`. The
#' model object is the last rep's `cv.glmnet`; its `lambda.min` is
#' overwritten with `lambda*` so [metabodecon::predict_lasso()] picks
#' the right lambda by default.
#' @param X Numeric feature matrix.
#' @param y Factor with two levels.
#' @param seed Random seed for the first rep's inner-CV fold assignment;
#'   subsequent reps use `seed+1, seed+2, …`.
#' @param nworkers When `> 1`, sets up a `doParallel` cluster and asks
#'   `cv.glmnet()` to parallelize its inner-CV folds via foreach. Falls
#'   back to single-threaded when `doParallel` isn't installed.
#' @param nreps Number of `cv.glmnet` repetitions used to estimate the
#'   reported acc/AUC. Default 5.
#' @return A list with `model` (a `cv.glmnet` object whose `lambda.min`
#'   has been overwritten with the ACC-maximizing (AUC-tiebroken)
#'   `lambda*`), `acc`, `auc`, `acc_se`, `auc_se`.
fit_lasso <- function(X, y, seed=1, nworkers=1L, nreps=5L) {
    requireNamespace("glmnet", quietly=TRUE)
    stopifnot(is_int(nreps, 1), nreps >= 1L)
    lvs <- levels(y)
    par_ok <- nworkers > 1L && requireNamespace("doParallel", quietly=TRUE)
    if (par_ok) {
        cl <- parallel::makeCluster(min(nworkers, 10L))
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add=TRUE)
    }
    # Rep 1 discovers the lambda path from the data; subsequent reps
    # reuse it so per-lambda OOF arrays line up across reps.
    cvs <- vector("list", nreps)
    set.seed(seed)
    cvs[[1]] <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1,
                                   nfolds=10, parallel=par_ok, keep=TRUE)
    lambda_path <- cvs[[1]]$lambda
    for (r in seq_len(nreps - 1L) + 1L) {
        set.seed(seed + r - 1L)
        cvs[[r]] <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1,
                                       nfolds=10, parallel=par_ok, keep=TRUE,
                                       lambda=lambda_path)
    }
    nl <- length(lambda_path)
    acc_mat <- matrix(NA_real_, nrow=nreps, ncol=nl)
    auc_mat <- matrix(NA_real_, nrow=nreps, ncol=nl)
    for (r in seq_len(nreps)) {
        pv <- cvs[[r]]$fit.preval[, seq_len(nl), drop=FALSE]
        for (j in seq_len(nl)) {
            link <- pv[, j]
            prob <- 1 / (1 + exp(-link))
            cls  <- factor(ifelse(prob > 0.5, lvs[2], lvs[1]), levels=lvs)
            acc_mat[r, j] <- mean(cls == y, na.rm=TRUE)
            auc_mat[r, j] <- AUC(y, prob)
        }
    }
    mean_acc <- colMeans(acc_mat, na.rm=TRUE)
    mean_auc <- colMeans(auc_mat, na.rm=TRUE)
    # Pick lambda* by averaged accuracy with AUC as tiebreaker. Search
    # from the largest-lambda end (cv.glmnet's path is decreasing, so
    # `rev()` puts the largest lambda first; `which.max` then breaks
    # remaining ties toward the more-regularized end, matching
    # cv.glmnet's lambda.1se sensibility).
    score <- mean_acc * 1e6 + mean_auc
    j_star <- which.max(rev(score))
    j_star <- nl - j_star + 1L
    chosen_lambda <- lambda_path[j_star]
    final_model <- cvs[[nreps]]
    final_model$lambda.min <- chosen_lambda
    list(
        model = final_model,
        acc = mean_acc[j_star],
        auc = mean_auc[j_star],
        acc_se = if (nreps >= 2L) stats::sd(acc_mat[, j_star], na.rm=TRUE) / sqrt(nreps) else NA_real_,
        auc_se = if (nreps >= 2L) stats::sd(auc_mat[, j_star], na.rm=TRUE) / sqrt(nreps) else NA_real_
    )
}

#' @noRd
#' @title Lasso predictor for fit_mdm
#' @description
#' Companion of [metabodecon::fit_lasso()]. Returns the positive-class
#' probability at `lambda.min` for each row of `newx`.
#' @param model Object returned in the `model` slot of
#'   [metabodecon::fit_lasso()].
#' @param newx Numeric feature matrix.
#' @return Numeric vector of length `nrow(newx)`.
predict_lasso <- function(model, newx) {
    requireNamespace("glmnet", quietly=TRUE)
    as.numeric(stats::predict(model, newx=newx, s="lambda.min", type="response"))
}

#' @noRd
#' @title Random-forest fitter for fit_mdm
#' @description
#' Fits a probability random forest with `num.trees` trees via
#' [ranger::ranger()] and reads off the OOB acc/AUC for the positive
#' class. OOB is an unbiased estimate of generalization performance and
#' converges as `num.trees` grows. The default `num.trees=5000` is
#' chosen high enough that OOB acc/AUC has converged for the sample
#' sizes mdm typically targets (50-500); a cheap convergence check is to
#' double `num.trees` and confirm OOB AUC does not move beyond its
#' standard error.
#' @param X Numeric feature matrix.
#' @param y Factor with two levels.
#' @param seed Random seed for ranger.
#' @param nworkers Forwarded to `ranger::ranger(num.threads=...)` for
#'   per-tree parallelism.
#' @param num.trees Number of ranger trees. Default 5000.
#' @param importance Forwarded to `ranger::ranger(importance=...)`.
#'   Default `"none"`. Set to `"permutation"` (or `"impurity"`) to
#'   populate `model$variable.importance` at fit time and avoid a
#'   redundant refit downstream; permutation importance roughly doubles
#'   the per-tree training cost.
#' @return A list with `model` (a `ranger` object with the trained
#'   levels stashed on `model$lvs` for the predict path), `acc`, `auc`.
#'   `acc_se` and `auc_se` are `NA` because OOB produces a single point
#'   estimate.
fit_ranger <- function(
    X, y, seed=1, nworkers=1L, num.trees=5000L, importance="none"
) {
    requireNamespace("ranger", quietly=TRUE)
    stopifnot(is_int(num.trees, 1), num.trees >= 1L)
    lvs <- levels(y)
    rf <- ranger::ranger(
        x=X, y=y, probability=TRUE, num.trees=num.trees,
        importance=importance,
        seed=seed, num.threads=max(1L, nworkers)
    )
    rf$lvs <- lvs
    oob <- rf$predictions[, lvs[2]]
    ok <- is.finite(oob)
    acc <- NA_real_; auc <- NA_real_
    if (any(ok)) {
        cls <- factor(ifelse(oob[ok] > 0.5, lvs[2], lvs[1]), levels=lvs)
        acc <- mean(cls == y[ok])
        auc <- AUC(y[ok], oob[ok])
    }
    list(model=rf, acc=acc, auc=auc, acc_se=NA_real_, auc_se=NA_real_)
}

#' @noRd
#' @title Random-forest predictor for fit_mdm
#' @description
#' Companion of [metabodecon::fit_ranger()]. Returns the positive-class
#' probability for each row of `newx`.
#' @param model Object returned in the `model` slot of [metabodecon::fit_ranger()].
#' @param newx Numeric feature matrix.
#' @return Numeric vector of length `nrow(newx)`.
predict_ranger <- function(model, newx) {
    requireNamespace("ranger", quietly=TRUE)
    colnames(newx) <- model$forest$independent.variable.names
    pm <- stats::predict(model, data=newx)$predictions
    pm[, model$lvs[2]]
}

# Helpers #####

get_mog <- function(npmax, maxShift, maxCombine) {
    g <- expand.grid(
        npmax=as.integer(npmax), maxShift=as.integer(maxShift),
        maxCombine=as.integer(maxCombine),
        KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE
    )
    g <- g[order(g$npmax, g$maxShift, g$maxCombine), , drop=FALSE]
    rownames(g) <- NULL
    g$acc <- NA_real_; g$auc <- NA_real_
    g$acc_se <- NA_real_; g$auc_se <- NA_real_
    g
}

# `find_npmax_elbow` / `find_npmax_elbow_one` live in R/decon.R since
# they back the `npmax = -1` (auto) / `npmax = -2` (intrinsic)
# resolution inside `deconvolute_spectra` / `deconvolute_spectrum`.

# Adaptive maxShift selection by dip detection. Sweeps maxShift through
# {1, 2, 4, 8, ...}, runs CluPA at each step, computes the average
# pairwise Pearson correlation of the aligned superpositions
# (`sit$supal`), and stops the FIRST time the correlation decreases
# compared to the previous step. Returns the maxShift from the step
# *before* the dip (the last one that was still improving). If no dip
# is seen by `max_cap`, returns `max_cap`. Requires a CluPA-compatible
# `align_fun` (writes `sit$supal`); defaults to [metabodecon::clupa()].
find_maxShift_dip <- function(d, align_fun=clupa, max_cap=512L,
                              nworkers=1, verbose=FALSE) {
    avg_pearson <- function(a) {
        M <- do.call(cbind, lapply(a, function(s) s$sit$supal))
        C <- stats::cor(M)
        mean(C[upper.tri(C)])
    }
    p_prev <- NA_real_; ms <- 1L
    repeat {
        a <- align_fun(x=d, ref=NULL, maxShift=ms,
                       verbose=verbose, nworkers=nworkers, full=TRUE)
        p <- avg_pearson(a)
        if (!is.na(p_prev) && p < p_prev) return(as.integer(ms %/% 2L))
        if (ms >= max_cap) return(as.integer(ms))
        p_prev <- p; ms <- ms * 2L
    }
}

# True when any element of a decons2 has zero peaks.
has_zero_peaks <- function(d) {
    if (!inherits(d, "decons2")) return(FALSE)
    any(vapply(d, function(o) nrow(o$lcpar) == 0L, logical(1)))
}

as_binary01 <- function(y) {
    lvs <- sort(unique(y))
    if (length(lvs) != 2) stop("y must have exactly 2 unique levels")
    as.integer(y == lvs[2])
}

get_test_ids <- function(nfolds=5, nsamples, seed=1, y=NULL) {
    # Vectorize over seed: returns a flat length(seed)*nfolds list of
    # test-id vectors, each carrying the originating `seed` as an
    # attribute so callers (e.g. `benchmark`) can stamp it back onto
    # per-fold output rows. `(seed, fold)` ordering is seed-major: all
    # `nfolds` folds of seed[1] first, then seed[2], etc.
    if (length(seed) > 1L) {
        out <- vector("list", length(seed) * nfolds)
        for (si in seq_along(seed)) {
            sub <- get_test_ids(nfolds=nfolds, nsamples=nsamples, seed=seed[si], y=y)
            for (i in seq_len(nfolds)) {
                attr(sub[[i]], "seed") <- seed[si]
                attr(sub[[i]], "fold") <- i
            }
            out[((si - 1L) * nfolds + 1L):(si * nfolds)] <- sub
        }
        return(out)
    }
    set.seed(seed)
    if (is.null(y)) {
        ids <- sample(seq_len(nsamples))
        grp <- split(ids, cut(seq_along(ids), nfolds, labels=FALSE))
        return(lapply(grp, sort))
    }

    y <- as_binary01(y)
    levs <- sort(unique(y))
    out <- vector("list", nfolds)
    for (k in seq_len(nfolds)) out[[k]] <- integer(0)

    for (lev in levs) {
        ids <- sample(which(y == lev))
        grp <- split(ids, cut(seq_along(ids), nfolds, labels=FALSE))
        for (k in seq_len(nfolds)) {
            out[[k]] <- c(out[[k]], grp[[k]])
        }
    }

    lapply(out, sort)
}

get_foldid <- function(y, nfolds=5, seed=1) {
    te_list <- get_test_ids(
        nfolds=nfolds, nsamples=length(y), seed=seed, y=y
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
AUC <- function(y, yhat) {
    y <- as_binary01(y)
    pos <- y == 1
    n1 <- sum(pos); n0 <- sum(!pos)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(yhat)
    (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

check_mdm_args <- function(
    x, y,
    sfr=NULL, igrs=list(), use_rust=NULL, nworkers=NULL,
    verbosity=NULL, seed=NULL
) {
    stopifnot(
        is_spectra(x),
        is.factor(y),
        length(y) == length(x),
        is_num_or_null(sfr, 2),
        is_list_of_nums(igrs, nv=2),
        is_bool_or_num(use_rust),
        is_int_or_null(nworkers, 1),
        is_int_or_null(verbosity, 1),
        is.null(seed) || is_int(seed)
    )
    if (!is.null(names(y)) && !identical(get_names(x), names(y))) {
        stop(
            "Names of `x` and `y` must match and be in the same order.",
            call.=FALSE
        )
    }
    if (nlevels(y) != 2 || any(table(y) == 0)) {
        stop("`y` must contain exactly 2 non-empty classes.", call.=FALSE)
    }
    invisible(NULL)
}

mdm_eval <- function(y, prob, lvs) {
    pred <- factor(ifelse(prob > 0.5, lvs[2], lvs[1]), levels=lvs)
    list(acc=mean(pred == y), auc=AUC(y, prob))
}

# Format a (mean, SE) pair as "nn.n(±x.x)%" with both quantities in
# percentage points. SE is omitted when NA.
fmt_pct_se <- function(m, s) {
    if (is.na(m)) return("NA")
    if (is.na(s)) return(sprintf("%.1f%%", m * 100))
    sprintf("%.1f(\u00b1%.1f)%%", m * 100, s * 100)
}

# Format a fraction as "nn.n%". NA -> "NA".
fmt_pct <- function(m) {
    if (is.na(m)) return("NA")
    sprintf("%.1f%%", m * 100)
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
#' spectra are deconvoluted, aligned and snapped to the references
#' stored in the model before prediction. When `newdata` is a numeric
#' matrix, it is used directly as the feature matrix.
#'
#' `print.mdm()` prints a compact model summary.
#'
#' `coef.mdm()` returns lasso coefficients (or ranger importance).
#'
#' `plot.mdm()` plots the lasso path (or ranger importance bars).
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
#' @param nworkers Number of workers to deconvolute and align `newdata`.
#' @param verbosity Integer verbosity level.
#' @param ... Passed to underlying methods where applicable.
#'
#' @return
#' - `predict`: numeric vector of probabilities, classes, and/or link scores.
#' - `print`: invisibly returns `x`.
#' - `coef`: coefficient object from `glmnet` (or ranger importance).
#' - `plot`: invisibly returns `NULL`.
#' - `summary`: object of class `summary.mdm`.
#' - `print.summary.mdm`: invisibly returns `x`.
#'
predict.mdm <- function(
    object, newdata,
    type=c("all", "prob", "class", "link"),
    nworkers=1, verbosity=1, ...
) {
    stopifnot(
        inherits(object, "mdm"), is_int(nworkers, 1),
        is_spectra(newdata) || is.matrix(newdata) || is.data.frame(newdata)
    )
    type <- match.arg(type)
    p <- object$params
    lvs <- p$lvs

    if (is.null(object$model)) {
        n <- if (is_spectra(newdata)) length(newdata)
             else nrow(as.matrix(newdata))
        z <- rep(0, n); h <- rep(0.5, n)
        cl <- factor(rep(NA_character_, n), levels=lvs)
        if (type == "all") return(data.frame(link=z, prob=h, class=cl))
        if (type == "class") return(cl)
        if (type == "prob") return(h)
        return(z)
    }

    if (is_spectra(newdata)) {
        decon_fun <- p$decon_fun %||% deconvolute
        logv("Deconvoluting %d spectra with %d nworkers",
             length(newdata), nworkers)
        d <- decon_fun(
            x=newdata, sfr=p$sfr, igrs=p$igrs %||% list(),
            verbose=verbosity >= 2, use_rust=p$use_rust,
            npmax=p$npmax, nworkers=nworkers
        )
        align_fun <- p$align_fun %||% clupa
        a_aligned <- align_fun(
            x=d, ref=object$ref$align, maxShift=p$maxShift,
            verbose=verbosity >= 2, nworkers=nworkers, full=FALSE
        )
        a <- p$snap_fun(
            a_aligned, ref=object$ref$snap,
            maxCombine=p$maxCombine, igrs=p$igrs %||% list()
        )
        Xn <- p$feat_fun(
            a, maxCombine=p$maxCombine,
            igrs=p$igrs %||% list(), peakPos=p$peakPos
        )
    } else {
        Xn <- as.matrix(newdata)
    }

    logv("Predicting with stored predict_fun")
    prob <- p$predict_fun(object$model, Xn)
    if (type == "prob") return(prob)
    cls <- factor(ifelse(prob > 0.5, lvs[2], lvs[1]), levels=lvs)
    if (type == "class") return(cls)
    eps <- .Machine$double.eps
    link <- log(pmin(1 - eps, pmax(eps, prob)) / pmin(1 - eps, pmax(eps, 1 - prob)))
    if (type == "link") return(link)
    out <- data.frame(link=link, prob=prob, class=cls)
    if (type == "all") return(out)
    prob
}

#' @export
#' @rdname mdm_methods
print.mdm <- function(x, ...) {
    stopifnot(inherits(x, "mdm"), is.list(x$params))
    pp <- c("npmax", "maxShift", "maxCombine")
    cat("metabodecon model (mdm)\n")
    cat("  ", formatC("model:", width=-15), paste(class(x$model), collapse=", "),
        "\n", sep="")
    for (nm in pp) {
        v <- x$params[[nm]]
        if (is.null(v)) next
        lab <- formatC(paste0(nm, ":"), width=-15)
        cat("  ", lab, v, "\n", sep="")
    }
    if (!is.null(x$auc)) {
        cat("  ", formatC("acc:", width=-15),
            fmt_pct_se(x$acc, x$acc_se), "\n", sep="")
        cat("  ", formatC("auc:", width=-15),
            fmt_pct_se(x$auc, x$auc_se), "\n", sep="")
    }
    invisible(x)
}

#' @export
#' @rdname mdm_methods
coef.mdm <- function(object, ...) {
    stopifnot(inherits(object, "mdm"), !is.null(object$model))
    if (inherits(object$model, "ranger")) {
        return(object$model$variable.importance)
    }
    stats::coef(object$model, s="lambda.min", ...)
}

#' @export
#' @rdname mdm_methods
plot.mdm <- function(x, ...) {
    stopifnot(inherits(x, "mdm"), !is.null(x$model))
    if (inherits(x$model, "ranger")) {
        vi <- sort(x$model$variable.importance, decreasing=TRUE)
        graphics::barplot(vi, las=2, ...)
        return(invisible(NULL))
    }
    graphics::plot(x$model, ...)
    invisible(NULL)
}

#' @export
#' @rdname mdm_methods
summary.mdm <- function(object, ...) {
    stopifnot(inherits(object, "mdm"), is.list(object$params))
    pp <- c("npmax", "maxShift", "maxCombine")
    out <- object$params[pp]
    out$model <- paste(class(object$model), collapse=", ")
    out$n_peaks <- length(object$params$peakPos %||% integer(0))
    out$acc <- fmt_pct_se(object$acc, object$acc_se)
    out$auc <- fmt_pct_se(object$auc, object$auc_se)
    class(out) <- "summary.mdm"
    out
}

#' @export
#' @rdname mdm_methods
print.summary.mdm <- function(x, ...) {
    stopifnot(inherits(x, "summary.mdm"))
    cat("Summary of mdm\n")
    for (nm in names(x)) {
        lab <- formatC(paste0(nm, ":"), width=-15)
        cat("  ", lab, x[[nm]], "\n", sep="")
    }
    invisible(x)
}
