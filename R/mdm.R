
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
#' An `mdm` is a binary classification model fitted on a feature matrix
#' built from NMR spectra by the following pipeline:
#'
#' \preformatted{
#'   x  --decon_fun-->  decons  --align_fun-->  aligns  --feat_fun-->  X
#'   X, y  --fit_fun-->  model  --predict_fun-->  prob
#' }
#'
#' All five stages are pluggable: pass `decon_fun=identity2` and/or
#' `align_fun=identity_align` to skip them and fit baselines that
#' operate directly on raw spectra (e.g. a binning model). The
#' `feat_fun`, `fit_fun` and `predict_fun` stages are pluggable as well.
#'
#' [metabodecon::fit_mdm()] iterates over the rows of a 'model fitting
#' grid' (`mog`) — see [metabodecon::get_mog()] — and for each row applies
#' the pipeline and records held-out accuracy and AUC. The best model and
#' the augmented grid are returned. A single-row `mog` degenerates to
#' fitting one model.
#'
#' [metabodecon::benchmark()] runs an outer k-fold cross-validation over
#' [metabodecon::fit_mdm()] to estimate end-to-end predictive performance.
#'
#' @details
#'
#' ## Pluggable interfaces
#'
#' The `feat_fun`, `fit_fun` and `predict_fun` functions speak a fixed
#' parameter vocabulary. Replacements must accept the listed arguments
#' (extras via `...`).
#'
#' \describe{
#'   \item{`feat_fun(x, maxCombine, peakPos, igrs, ...)`}{Returns a
#'         numeric matrix with one row per spectrum.}
#'   \item{`fit_fun(X, y, seed, nworkers)`}{Trains and returns a model.
#'         `y` is a factor with the two class levels (positive class is
#'         `levels(y)[2]`). `nworkers` requests intra-fit parallelism
#'         (e.g. `ranger`'s `num.threads`, `cv.glmnet`'s parallel CV);
#'         implementations may ignore it. Must call
#'         `requireNamespace("<backend>")` so it works after a fresh
#'         `readRDS()`.}
#'   \item{`predict_fun(model, newx)`}{Returns a numeric vector of
#'         positive-class probabilities, one per row of `newx`. Must
#'         also call `requireNamespace("<backend>")`.}
#' }
#'
#' Built-in implementations: [metabodecon::peak_mat()],
#' [metabodecon::si_mat()], [metabodecon::bin()],
#' [metabodecon::fit_lasso()] / [metabodecon::predict_lasso()] and
#' [metabodecon::fit_ranger500()] / [metabodecon::predict_ranger500()].
#'
#' ## Caching within `fit_mdm`
#'
#' Rows of `mog` are sorted by `(npmax, nfit, smit, smws, delta, maxShift,
#' maxCombine)` so identical decon-tuples and align-tuples cluster.
#' Inside the loop, the most recent deconvolution and alignment are kept
#' and reused whenever the relevant subset of parameters is unchanged.
#'
#' ## Grid-search results carried by spectra
#'
#' When any row of `mog` has `npmax > 0` and `decon_fun` is not
#' [metabodecon::identity2()], [metabodecon::fit_mdm()] calls
#' [metabodecon::grid_deconvolute_spectra()] once up front to attach a
#' `$deg` element to each spectrum. The enriched spectra are reused
#' across rows and outer folds. [metabodecon::benchmark()] does the same
#' up-front attachment.
#'
#' ## Parallelism
#'
#' [metabodecon::benchmark()] runs outer folds sequentially and delegates
#' all parallelism to the inner fitter via `nworkers`.
#'
#' @param x Spectra object.
#' @param y Factor vector with class labels for each spectrum.
#' @param mog Model-fitting grid as returned by [metabodecon::get_mog()].
#' @param deg Deconvolution-parameter grid forwarded to
#'   [metabodecon::grid_deconvolute_spectra()] when any row of `mog` has
#'   `npmax > 0`. When `NULL` (default), the default grid built into
#'   [metabodecon::grid_deconvolute_spectra()] is used.
#' @param sfr Signal-free region. See [metabodecon::deconvolute()].
#' @param use_rust Use the Rust backend?
#' @param nworkers Number of workers for deconvolution and alignment.
#' @param verbosity Verbosity level.
#' @param seed Random seed for fold assignments and `fit`.
#' @param nfolds Number of folds for the inner CV used by `fit`.
#' @param check Validate inputs at function entry?
#' @param decon_fun Deconvolution function. Must accept the parameter
#'   vocabulary used by [metabodecon::deconvolute()] (i.e. `x`, `sfr`,
#'   `igrs`, `use_rust`, `nfit`, `smit`, `smws`, `delta`, `npmax`,
#'   `nworkers`, `verbose`). Pass [metabodecon::identity2()] to skip
#'   deconvolution and feed the raw spectra into the next stage.
#' @param align_fun Alignment function called after deconvolution. Must
#'   have signature `f(x, ref, maxShift, verbose, nworkers, full, ...)`.
#'   Built-in choices: [metabodecon::clupa()] (default, CluPA
#'   hierarchical-clustering peak alignment), [metabodecon::vopa()]
#'   (vote-based peak alignment, faster), [metabodecon::glopa()] (single
#'   global integer shift via FFT cross-correlation),
#'   [metabodecon::identity_align()] (no alignment).
#' @param feat_fun Feature-matrix function. See *Pluggable interfaces*.
#' @param fit_fun Inner-model fitter. See *Pluggable interfaces*.
#' @param predict_fun Inner-model predictor paired with `fit_fun`. See
#'   *Pluggable interfaces*.
#' @param igrs Ignore regions in ppm.
#' @param maxCombine Bin width / peak-merging tolerance for
#'   `peak_mat()` / `si_mat()` / `bin()`.
#' @param k Number of outer folds for [metabodecon::benchmark()].
#' @param conf Character string selecting a predefined `mog` configuration.
#' @param x_tr Training spectra for [metabodecon::bootstrap_mdm_one()].
#' @param y_tr Training labels (factor) for [metabodecon::bootstrap_mdm_one()].
#' @param x_te Test spectra for [metabodecon::bootstrap_mdm_one()].
#' @param y_te Test labels (factor) for [metabodecon::bootstrap_mdm_one()].
#' @param row Single-row data frame (one row of `mog`) passed to
#'   [metabodecon::bootstrap_mdm_one()].
#' @param test_frac Fraction of spectra to hold out as the test set per round.
#'   Actual count is `max(round(n * test_frac), 5)` so small datasets still
#'   have a usable test set.
#' @param max_rounds Maximum number of bootstrap rounds.
#' @param se_thresh Stop early when the SE of the best AUC (across completed
#'   rounds) falls below this. Set to `0` to always run `max_rounds` rounds.
#' @param nrounds For [metabodecon::fit_mdm()]: number of outer-CV rounds
#'   (each round is a fresh stratified `nfolds`-fold split). Also used by
#'   [metabodecon::cv_mdm()] as the number of repeated `kfold`-fold CVs.
#' @param kfold Number of folds per CV for [metabodecon::cv_mdm()].
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with elements
#' `model` (best fitted backend model, refit on the full data), `ref`
#' (reference alignment spectrum, `NULL` when `align_fun=identity_align`),
#' `params` (everything needed to reproduce predictions: chosen grid row,
#' `feat_fun`, `predict_fun`, `lvs`, `peakPos`, `sfr`, `igrs`, `use_rust`,
#' `decon_fun`, `align_fun`) and `mog` (input grid augmented with `acc`,
#' `acc_se`, `auc` and `auc_se` columns). Per-row scoring uses `nrounds`
#' repeats of `nfolds`-fold outer CV, with deconvolution and alignment
#' reused across folds (mild over-optimism vs. a full `cv_mdm()`
#' re-decon/re-align per fold, but much faster); `acc`/`auc` are means
#' of the per-fold acc/AUC over all `nrounds * nfolds` folds, and the
#' `_se` columns are the across-folds SE of that mean (`NA` when
#' `nrounds * nfolds < 2`).
#'
#' [metabodecon::benchmark()] returns a list with elements:
#' - `models`: list of fitted models, one per outer fold.
#' - `predictions`: data frame with columns `fold`, `true`, `link`, `prob`,
#'   `pred`.
#' - `performance`: data frame with per-fold `acc` and `auc`.
#' - `overall`: list with pooled `acc` and `auc`.
#'
#' [metabodecon::bootstrap_mdm_one()] returns a list with `prob` (named
#' numeric vector of test-set class probabilities), `acc` and `auc`.
#'
#' [metabodecon::bootstrap_mdm()] returns a list with:
#' - `best_row`: 1-row data frame, the `mog` row with the highest mean AUC.
#' - `mog`: grid augmented with mean `acc`, `auc` and `auc_se` columns.
#' - `probs`: list of length `nrow(mog)`; each element is a list of
#'   per-round probability vectors (named by test-sample name).
#' - `acc_mat`, `auc_mat`: numeric matrices (rounds × `nrow(mog)`).
#'
#' [metabodecon::cv_mdm()] returns the same shape as
#' [metabodecon::bootstrap_mdm()], with rows of `acc_mat`/`auc_mat`
#' corresponding to the `nrounds` repeated k-fold CVs and each per-row
#' AUC/acc computed from the pooled out-of-fold predictions of one
#' k-fold partition.
#'
#' @examples
#' \dontrun{
#'   m <- fit_mdm(spectra, y, mog=get_mog("default"))
#'   bm <- benchmark(spectra, y, k=5, mog=get_mog("default"))
#'   bs <- bootstrap_mdm(spectra, y, mog=get_mog("default"), max_rounds=5)
#'   cv <- cv_mdm(spectra, y, mog=get_mog("default"), nrounds=1, kfold=10)
#' }
#'
fit_mdm <- function(x, y,
    decon_fun=deconvolute, align_fun=clupa,           feat_fun=peak_mat,
    fit_fun=fit_lasso,     predict_fun=predict_lasso, mog=get_mog("default"),
    deg=NULL,              sfr=NULL,                  igrs=list(),
    use_rust=0,            nworkers=1,                verbosity=1,
    seed=1,                nfolds=5,                  check=TRUE,
    nrounds=5
) {
    if (check) check_mdm_args(
        x=x, y=y, mog=mog, sfr=sfr, igrs=igrs,
        use_rust=use_rust, nworkers=nworkers, verbosity=verbosity,
        seed=seed, nfolds=nfolds
    )
    stopifnot(
        is.function(feat_fun), is.function(fit_fun),
        is.function(predict_fun), is.function(decon_fun),
        is.function(align_fun),
        is_int(nrounds, 1), nrounds >= 1
    )
    skip_decon <- identical(decon_fun, identity2)
    decon_fn <- if (identical(decon_fun, deconvolute)) deconvolute_spectra else decon_fun
    lvs <- levels(y)

    # Sort rows so identical decon/align tuples cluster.
    ord <- with(mog, order(npmax, nfit, smit, smws, delta, maxShift, maxCombine))
    mog <- mog[ord, , drop=FALSE]
    rownames(mog) <- NULL
    mog$acc <- NA_real_
    mog$auc <- NA_real_
    mog$acc_se <- NA_real_
    mog$auc_se <- NA_real_

    # Pre-attach per-spectrum `$deg` tables when any row uses npmax > 0.
    if (!skip_decon && any(mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2,
            nworkers=min(nworkers, length(x)), use_rust=use_rust
        )
    }

    nr <- nrow(mog); ns <- length(x)
    logv("Starting grid search (%d combinations, %d spectra, %d rounds x %d-fold CV)",
         nr, ns, nrounds, nfolds)
    last_dkey <- NULL; last_akey <- NULL
    d <- NULL; a <- NULL; a_clupa <- NULL
    best_mdm <- NULL; best_auc <- -Inf
    for (i in seq_len(nr)) {
        r <- mog[i, , drop=FALSE]
        dkey <- list(r$npmax, r$nfit, r$smit, r$smws, r$delta)
        if (!identical(dkey, last_dkey)) {
            d <- decon_fn(
                x=x, sfr=sfr, igrs=igrs, verbose=verbosity >= 2,
                use_rust=use_rust, nfit=r$nfit, smit=r$smit, smws=r$smws,
                delta=r$delta, npmax=r$npmax, nworkers=nworkers
            )
            if (inherits(d, "decons2")) {
                nps <- vapply(d, function(o) nrow(o$lcpar), integer(1))
                if (any(nps == 0)) {
                    logv("[%d/%d] %d spectra produced zero peaks; skipping",
                         i, nr, sum(nps == 0))
                    mog$acc[i] <- 0; mog$auc[i] <- 0
                    last_dkey <- dkey; last_akey <- NULL
                    next
                }
            }
            last_dkey <- dkey
            last_akey <- NULL
        }

        akey <- list(dkey, r$maxShift)
        if (!identical(akey, last_akey)) {
            a_clupa <- align_fun(x=d, ref=NULL, maxShift=r$maxShift, verbose=verbosity >= 2, nworkers=nworkers, full=FALSE)
            last_akey <- akey
        }

        # RefPA snap (cheap peak-list op; not cached). Skipped when
        # align_fun is identity_align (then a_clupa is decons2/spectra).
        if (inherits(a_clupa, "aligns")) {
            a <- snap_to_ref(a_clupa, ref=NULL, maxCombine=r$maxCombine)
        } else {
            a <- a_clupa
        }

        # Outer CV: flatten nrounds * nfolds test folds into one iteration.
        # Score each fold's model on its own test fold (per-fold AUC) and
        # average across all folds for a single mean +- SE. Pooling is
        # avoided so that score-scale differences between fold models do
        # not depress AUC. A progress line is emitted at every round
        # boundary so the caller can watch the mean +- SE converge.
        all_folds <- do.call(c, lapply(seq_len(nrounds), function(rnd) {
            get_test_ids(nfolds=nfolds, nsamples=ns, seed=seed+rnd, y=y)
        }))
        n_total <- length(all_folds)
        fold_accs <- rep(NA_real_, n_total)
        fold_aucs <- rep(NA_real_, n_total)
        m_acc <- NA_real_; m_auc <- NA_real_
        se_acc <- NA_real_; se_auc <- NA_real_
        params_str <- sprintf("p=%d S=%d C=%g", r$npmax, r$maxShift, r$maxCombine)
        for (idx in seq_len(n_total)) {
            te <- all_folds[[idx]]
            rnd <- (idx - 1L) %/% nfolds + 1L
            fld <- (idx - 1L) %% nfolds + 1L
            tr <- setdiff(seq_len(ns), te)
            a_tr <- a[tr]; class(a_tr) <- class(a)
            a_te <- a[te]; class(a_te) <- class(a)
            X_tr_full <- feat_fun(a_tr, maxCombine=r$maxCombine, igrs=igrs)
            pp <- which(colSums(X_tr_full != 0) > 0)
            if (length(pp) == 0L) next
            X_tr <- X_tr_full[, pp, drop=FALSE]
            X_te <- feat_fun(a_te, maxCombine=r$maxCombine, igrs=igrs)[, pp, drop=FALSE]
            model <- fit_fun(X_tr, y[tr], seed=seed+rnd, nworkers=nworkers)
            prob_te <- predict_fun(model, X_te)
            y_te <- y[te]
            pred_te <- factor(ifelse(prob_te > 0.5, lvs[2], lvs[1]), levels=lvs)
            fold_accs[idx] <- mean(pred_te == y_te)
            fold_aucs[idx] <- AUC(y_te, prob_te)
            if (fld != nfolds) next
            done <- seq_len(idx)
            n_ok <- sum(!is.na(fold_aucs[done]))
            if (n_ok > 0) {
                m_acc <- mean(fold_accs[done], na.rm=TRUE)
                m_auc <- mean(fold_aucs[done], na.rm=TRUE)
            }
            if (n_ok >= 2) {
                se_acc <- stats::sd(fold_accs[done], na.rm=TRUE) / sqrt(n_ok)
                se_auc <- stats::sd(fold_aucs[done], na.rm=TRUE) / sqrt(n_ok)
            }
            is_last_round <- rnd == nrounds
            is_best <- is_last_round && !is.na(m_auc) && m_auc > best_auc
            sym <- if (is_best) " <-- BEST" else ""
            header <- if (nrounds == 1L) sprintf("c=%d/%d", i, nr) else sprintf("c=%d/%d r=%d/%d", i, nr, rnd, nrounds)
            logv("%s %s acc=%s auc=%s%s", header, params_str,
                 fmt_pct(m_acc), fmt_pct(m_auc), sym)
            if (!is_last_round) next
            mog$acc[i] <- m_acc; mog$auc[i] <- m_auc
            mog$acc_se[i] <- se_acc; mog$auc_se[i] <- se_auc
            if (!is_best) next
            best_auc <- m_auc
            # Full-data refit: needed only for the final returned mdm. We do
            # it inline at each new best so the cheaper alternative (defer to
            # after the grid loop) doesn't require re-running decon+align for
            # the winning row.
            ref <- if (inherits(a, "aligns")) find_ref(a) else NULL
            mat <- feat_fun(a, maxCombine=r$maxCombine, igrs=igrs)
            peakPos <- which(colSums(mat != 0) > 0)
            full_model <- fit_fun(mat[, peakPos, drop=FALSE], y, seed=seed, nworkers=nworkers)
            params <- list(
                feat_fun=feat_fun, predict_fun=predict_fun, lvs=lvs,
                decon_fun=decon_fun, align_fun=align_fun,
                sfr=sfr, igrs=igrs, use_rust=use_rust,
                npmax=r$npmax, nfit=r$nfit, smit=r$smit, smws=r$smws,
                delta=r$delta, maxShift=r$maxShift, maxCombine=r$maxCombine,
                peakPos=peakPos
            )
            best_mdm <- structure(
                list(model=full_model, ref=ref, params=params), class="mdm"
            )
        }
    }

    best_mdm$mog <- mog
    ibest <- which.max(mog$auc)
    logv(
        "Best c=%d/%d: acc=%s auc=%s", ibest, nr,
        fmt_pct_se(mog$acc[ibest], mog$acc_se[ibest]),
        fmt_pct_se(mog$auc[ibest], mog$auc_se[ibest])
    )
    best_mdm
}

#' @export
#' @rdname mdm
#'
#' @title Fit an mdm with ranger using out-of-bag error
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Streamlined variant of [metabodecon::fit_mdm()] for ranger
#' probability forests. Replaces the `nrounds * nfolds` outer CV with
#' ranger's out-of-bag (OOB) error, which is an unbiased estimate of
#' generalization error and converges as `num.trees` grows.
#' Deconvolution and alignment are cached across grid rows in the same
#' way as [metabodecon::fit_mdm()]. Selects the best row by OOB
#' accuracy at the 0.5 threshold (AUC breaks ties).
#'
#' The grid is built internally as the full cross product of `npmax`,
#' `maxShift` and `maxCombine`. Deconvolution uses the defaults of
#' [metabodecon::deconvolute()] (`nfit=3, smit=2, smws=5, delta=6.4`);
#' alignment is always [metabodecon::clupa()].
#'
#' Any **negative** value in `maxCombine` acts as an on/off switch: it is
#' replaced row-wise by that row's `maxShift`. Pass `maxCombine=-1` to
#' tie peak-snapping width to the alignment tolerance and collapse the
#' grid from a 2D to a 1D sweep over `maxShift`.
#'
#' @param num.trees Number of ranger trees. Defaults to `5000`, which
#'   gives a well-converged OOB error for typical mdm sample sizes.
#'   Doubling `num.trees` and confirming the OOB AUC does not move
#'   beyond its standard error is a cheap convergence check.
#'
#' @return An object of class `mdm` with the same shape as
#'   [metabodecon::fit_mdm()]. `mog$acc_se` and `mog$auc_se` are `NA`
#'   because OOB produces a single point estimate per grid row.
#'
fit_mdm2 <- function(x, y,
    npmax=1000,       maxShift=2^(1:8), maxCombine=0,
    deg=NULL,         sfr=NULL,         igrs=list(),  use_rust=FALSE,
    use_speaq=FALSE,  num.trees=5000,   nworkers=1,   verbosity=1,
    seed=1,           check=TRUE,       supsh="lorentz",
    post_align=c("snap_to_ref", "combine_peaks")
) {
    supsh <- match.arg(supsh, c("lorentz", "sparse"))
    post_align <- match.arg(post_align)
    post_align_fun <- switch(post_align,
        snap_to_ref=snap_to_ref, combine_peaks=combine_peaks
    )
    if (check) {
        stopifnot(
            is_spectra(x), is.factor(y), length(y) == length(x),
            is_num(npmax), all(npmax >= 0),
            is_num(maxShift), all(maxShift >= 0),
            is_num(maxCombine),
            is_num_or_null(sfr, 2), is_list_of_nums(igrs, nv=2),
            is_bool_or_num(use_rust), is_bool(use_speaq),
            is_int(num.trees, 1), num.trees >= 1,
            is_int(nworkers, 1), is_int(verbosity, 1), is_int(seed, 1)
        )
        if (nlevels(y) != 2 || any(table(y) == 0)) {
            stop("`y` must contain exactly 2 non-empty classes.", call.=FALSE)
        }
        if (!is.null(names(y)) && !identical(get_names(x), names(y))) {
            stop("Names of `x` and `y` must match and be in the same order.",
                 call.=FALSE)
        }
    }
    requireNamespace("ranger", quietly=TRUE)
    lvs <- levels(y)

    mog <- expand.grid(npmax=npmax, maxShift=maxShift, maxCombine=maxCombine, KEEP.OUT.ATTRS=FALSE)
    # Negative maxCombine acts as a switch: use the row's maxShift instead.
    neg <- mog$maxCombine < 0
    mog$maxCombine[neg] <- mog$maxShift[neg]
    ord <- with(mog, order(npmax, maxShift, maxCombine))
    mog <- mog[ord, , drop=FALSE]
    rownames(mog) <- NULL
    mog$acc <- NA_real_; mog$auc <- NA_real_
    mog$acc_se <- NA_real_; mog$auc_se <- NA_real_

    if (any(mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2,
            nworkers=min(nworkers, length(x)), use_rust=use_rust
        )
    }

    nr <- nrow(mog); ns <- length(x)
    logv("fit_mdm2: %d combos, %d spectra, ranger OOB (num.trees=%d)", nr, ns, num.trees)
    last_dkey <- NULL; last_akey <- NULL
    d <- NULL; a <- NULL; a_clupa <- NULL
    best_mdm <- NULL; best_acc <- -Inf; best_tie_auc <- -Inf
    for (i in seq_len(nr)) {
        r <- mog[i, , drop=FALSE]
        dkey <- r$npmax
        if (!identical(dkey, last_dkey)) {
            d <- deconvolute_spectra(
                x=x, sfr=sfr, igrs=igrs, verbose=verbosity >= 2,
                use_rust=use_rust, npmax=r$npmax, nworkers=nworkers,
                full=FALSE, supsh=supsh
            )
            nps <- vapply(d, function(o) nrow(o$lcpar), integer(1))
            if (any(nps == 0)) {
                logv("[c=%d/%d] %d spectra produced zero peaks; skipping",
                     i, nr, sum(nps == 0))
                mog$acc[i] <- 0; mog$auc[i] <- 0
                last_dkey <- dkey; last_akey <- NULL
                next
            }
            last_dkey <- dkey; last_akey <- NULL
        }

        akey <- list(dkey, r$maxShift)
        if (!identical(akey, last_akey)) {
            a_clupa <- clupa(x=d, ref=NULL, maxShift=r$maxShift,
                             verbose=verbosity >= 2, nworkers=nworkers,
                             full=FALSE, use_speaq=use_speaq,
                             supsh=supsh)
            last_akey <- akey
        }
        a <- post_align_fun(a_clupa, ref=NULL, maxCombine=r$maxCombine)

        X <- peak_mat(a, igrs=igrs)
        pp <- which(colSums(X != 0) > 0)
        if (length(pp) == 0L) next
        Xp <- X[, pp, drop=FALSE]

        rf <- ranger::ranger(
            x=Xp, y=y, probability=TRUE, num.trees=num.trees,
            seed=seed, num.threads=max(1L, nworkers)
        )
        rf$lvs <- lvs
        class(rf) <- c("ranger500", class(rf))

        # OOB probability for the positive class. Samples that were
        # in-bag for every tree get NaN; with num.trees=5000 this is
        # vanishingly rare but we still filter to be safe.
        oob <- rf$predictions[, lvs[2]]
        ok <- is.finite(oob)
        m_acc <- NA_real_; m_auc <- NA_real_
        if (any(ok)) {
            cls <- factor(ifelse(oob[ok] > 0.5, lvs[2], lvs[1]), levels=lvs)
            m_acc <- mean(cls == y[ok])
            m_auc <- AUC(y[ok], oob[ok])
        }
        mog$acc[i] <- m_acc; mog$auc[i] <- m_auc

        # Optimize accuracy; break ties on AUC.
        is_best <- !is.na(m_acc) && (
            m_acc > best_acc ||
            (m_acc == best_acc && !is.na(m_auc) && m_auc > best_tie_auc)
        )
        sym <- if (is_best) " <-- BEST" else ""
        logv("c=%d/%d p=%d S=%d C=%g acc=%s auc=%s%s",
             i, nr, r$npmax, r$maxShift, r$maxCombine,
             fmt_pct(m_acc), fmt_pct(m_auc), sym)
        if (!is_best) next
        best_acc <- m_acc; best_tie_auc <- if (is.na(m_auc)) -Inf else m_auc
        ref <- if (inherits(a, "aligns")) find_ref(a) else NULL
        params <- list(
            feat_fun=peak_mat, predict_fun=predict_ranger500, lvs=lvs,
            decon_fun=deconvolute, align_fun=clupa,
            post_align=post_align,
            sfr=sfr, igrs=igrs, use_rust=use_rust,
            npmax=r$npmax, nfit=3, smit=2, smws=5, delta=6.4,
            maxShift=r$maxShift, maxCombine=r$maxCombine,
            peakPos=pp
        )
        best_mdm <- structure(
            list(model=rf, ref=ref, params=params), class="mdm"
        )
    }

    best_mdm$mog <- mog
    # Pick the final best by acc, AUC as tiebreaker.
    ord_final <- order(-mog$acc, -mog$auc, na.last=TRUE)
    ibest <- ord_final[1]
    logv("Best c=%d/%d: acc=%s auc=%s", ibest, nr,
         fmt_pct(mog$acc[ibest]), fmt_pct(mog$auc[ibest]))
    best_mdm
}

#' @export
#' @rdname mdm
#'
#' @title Repeated k-fold CV benchmark of fit_mdm2
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Slim outer-CV harness hardwired to [metabodecon::fit_mdm2()] (ranger
#' + clupa + peak_mat snapping). Runs one round of stratified
#' `nfolds`-fold CV per entry of `seeds`: in every fold a fresh
#' `fit_mdm2()` is trained on the train ids and scored on the held-out
#' ids via [stats::predict()] on the returned `mdm`. Spectra are
#' deconvoluted once up front so every fold reuses the same `$deg`
#' cache.
#'
#' Per-round mean and SD are taken across the `nfolds` per-fold values.
#' Overall mean and SE are taken across the `length(seeds)` per-round
#' means (`SE = sd(round_means) / sqrt(length(seeds))`).
#'
#' @return A list with elements:
#' - `acc_mean`, `acc_se`, `auc_mean`, `auc_se`: overall point
#'   estimates and standard errors.
#' - `round_acc`, `round_auc`: length-`length(seeds)` per-round means.
#' - `round_acc_sd`, `round_auc_sd`: length-`length(seeds)` within-round SDs.
#' - `fold_acc`, `fold_auc`: `length(seeds) × nfolds` matrices.
#' - `predictions`: data frame with `round`, `seed`, `fold`, `true`,
#'   `prob`, `pred`.
#' - `seeds`, `nfolds`: echo of the inputs.
#'
benchmark2 <- function(x, y,
    npmax=1000,       maxShift=2^(1:8), maxCombine=0,
    deg=NULL,         sfr=NULL,         igrs=list(),  use_rust=FALSE,
    use_speaq=FALSE,  num.trees=5000,   nworkers=1,   verbosity=1,
    seeds=1:10,       nfolds=10,        check=TRUE,    supsh="lorentz"
) {
    supsh <- match.arg(supsh, c("lorentz", "sparse"))
    if (check) {
        stopifnot(
            is_spectra(x), is.factor(y), length(y) == length(x),
            is_num(npmax), all(npmax >= 0),
            is_num(maxShift), all(maxShift >= 0),
            is_num(maxCombine),
            is_num_or_null(sfr, 2), is_list_of_nums(igrs, nv=2),
            is_bool_or_num(use_rust), is_bool(use_speaq),
            is_int(num.trees, 1), num.trees >= 1,
            is_int(nworkers, 1), is_int(verbosity, 1),
            is_num(seeds), length(seeds) >= 1L, !anyDuplicated(seeds),
            is_int(nfolds, 1), nfolds >= 2, nfolds <= length(y)
        )
        if (nlevels(y) != 2 || any(table(y) == 0)) {
            stop("`y` must contain exactly 2 non-empty classes.", call.=FALSE)
        }
    }
    requireNamespace("ranger", quietly=TRUE)

    # Pre-attach $deg once so every fold's fit_mdm2 reuses the cache.
    if (any(npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs, verbose=verbosity>=2,
            nworkers=min(nworkers, length(x)), use_rust=use_rust
        )
    }

    ns <- length(x); inner_v <- max(0L, verbosity - 1L)
    nr <- length(seeds)
    fold_acc <- matrix(NA_real_, nrow=nr, ncol=nfolds)
    fold_auc <- matrix(NA_real_, nrow=nr, ncol=nfolds)
    preds <- vector("list", nr * nfolds)
    logv("benchmark2: %d rounds x %d folds (num.trees=%d)", nr, nfolds, num.trees)
    for (r in seq_along(seeds)) {
        s <- seeds[r]
        te_list <- get_test_ids(nfolds=nfolds, nsamples=ns, seed=s, y=y)
        for (f in seq_along(te_list)) {
            te <- te_list[[f]]; tr <- setdiff(seq_len(ns), te)
            m <- fit_mdm2(
                x=x[tr], y=y[tr],
                npmax=npmax, maxShift=maxShift, maxCombine=maxCombine,
                deg=deg, sfr=sfr, igrs=igrs, use_rust=use_rust,
                use_speaq=use_speaq, num.trees=num.trees,
                nworkers=nworkers, verbosity=inner_v, seed=s,
                check=FALSE, supsh=supsh
            )
            p <- stats::predict(m, x[te], type="all", nworkers=nworkers, verbosity=inner_v)
            fold_acc[r, f] <- mean(p$class == y[te])
            fold_auc[r, f] <- AUC(y[te], p$prob)
            row <- data.frame(round=r, seed=s, fold=f, true=y[te], prob=p$prob, pred=p$class)
            preds[[((r-1) * nfolds) + f]] <- row
            acc <- fmt_pct(fold_acc[r, f])
            auc <- fmt_pct(fold_auc[r, f])
            logv("[r=%d/%d f=%d/%d] acc=%s auc=%s", r, nr, f, nfolds, acc, auc)
        }
        acc <- fmt_pct(mean(fold_acc[r, ], na.rm=TRUE))
        auc <- fmt_pct(mean(fold_auc[r, ], na.rm=TRUE))
        logv("[r=%d/%d s=%d] mean acc=%s auc=%s", r, nr, s, acc, auc)
    }
    round_acc <- rowMeans(fold_acc); round_auc <- rowMeans(fold_auc)
    round_acc_sd <- apply(fold_acc, 1, stats::sd)
    round_auc_sd <- apply(fold_auc, 1, stats::sd)
    acc_mean <- mean(round_acc); auc_mean <- mean(round_auc)
    acc_se <- if (nr >= 2) stats::sd(round_acc) / sqrt(nr) else NA_real_
    auc_se <- if (nr >= 2) stats::sd(round_auc) / sqrt(nr) else NA_real_
    logv("Overall: acc=%s auc=%s",
         fmt_pct_se(acc_mean, acc_se), fmt_pct_se(auc_mean, auc_se))
    list(
        acc_mean=acc_mean, acc_se=acc_se,
        auc_mean=auc_mean, auc_se=auc_se,
        round_acc=round_acc, round_acc_sd=round_acc_sd,
        round_auc=round_auc, round_auc_sd=round_auc_sd,
        fold_acc=fold_acc, fold_auc=fold_auc,
        predictions=do.call(rbind, preds),
        seeds=seeds, nfolds=nfolds
    )
}

#' @export
#' @rdname mdm
get_mog <- function(conf="default") {
    g <- expand.grid2(
        nfit = switch(conf, dynamic=0, 5),
        smit = switch(conf, dynamic=0, 2),
        smws = switch(conf, dynamic=0, static=c(3,5,7,9), 5),
        delta = switch(conf, dynamic=0, static=(1:5)*1.6, 6.4),
        npmax = switch(conf, dynamic=2^(6:11), 0),
        maxShift = switch(conf, default=50, 2^(1:8)),
        maxCombine = switch(conf, default=5, 2^(1:6))
    )
    ord <- order(g$npmax, g$nfit, g$smit, g$smws, g$delta, g$maxShift, g$maxCombine)
    g <- g[ord, , drop=FALSE]
    rownames(g) <- NULL
    g
}

#' @export
#' @rdname mdm
benchmark <- function(x, y,
    feat_fun=peak_mat,      fit_fun=fit_lasso,    predict_fun=predict_lasso,
    decon_fun=deconvolute,  align_fun=clupa,      mog=get_mog("default"),
    deg=NULL,               sfr=NULL,             igrs=list(),
    use_rust=0,             nworkers=1,           verbosity=2,
    seed=1,                 nfolds=5,             check=TRUE,
    k=3,                    nrounds=5
) {
    if (check) check_mdm_args(
        x=x, y=y, mog=mog, sfr=sfr, igrs=igrs,
        use_rust=use_rust, nworkers=nworkers, verbosity=verbosity,
        seed=seed, nfolds=nfolds
    )
    stopifnot(
        is.function(feat_fun), is.function(fit_fun),
        is.function(predict_fun), is.function(decon_fun),
        is.function(align_fun),
        is_int(k, 1), k >= 2, k <= length(y)
    )

    # One-time grid attach when any row of mog uses npmax > 0.
    if (!identical(decon_fun, identity2) && any(mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2, nworkers=nworkers, use_rust=use_rust
        )
    }

    # Forward verbosity-1 to the inner fitter.
    inner_verb <- max(0L, verbosity - 1L)
    te_list <- get_test_ids(nfolds=k, nsamples=length(x), seed=seed, y=y)
    models <- vector("list", k)
    fold_preds <- vector("list", k)
    perf <- data.frame(fold=integer(0), acc=numeric(0), auc=numeric(0))
    logv("Running %d-fold outer CV with fit_mdm", k)
    for (i in seq_along(te_list)) {
        te <- te_list[[i]]
        tr <- setdiff(seq_along(x), te)
        logv("[fold %d/%d] fitting", i, k)
        m <- fit_mdm(
            x=x[tr], y=y[tr],
            feat_fun=feat_fun, fit_fun=fit_fun, predict_fun=predict_fun,
            decon_fun=decon_fun, align_fun=align_fun, mog=mog,
            deg=deg, sfr=sfr, igrs=igrs,
            use_rust=use_rust, nworkers=nworkers, verbosity=inner_verb,
            seed=seed, nfolds=nfolds, nrounds=nrounds, check=FALSE
        )
        p <- stats::predict(m, x[te], type="all", nworkers=nworkers, verbosity=inner_verb)
        fp <- data.frame(fold=i, true=y[te], link=p$link, prob=p$prob, pred=p$class)
        fold_preds[[i]] <- fp
        acc <- mean(fp$pred == fp$true, na.rm=TRUE)
        auc <- AUC(fp$true, fp$prob)
        perf <- rbind(perf, data.frame(fold=i, acc=acc, auc=auc))
        logv("[fold %d/%d] acc=%.2f%% auc=%.4f", i, k, acc * 100, auc)
        models[[i]] <- m
    }

    preds <- do.call(rbind, fold_preds)
    overall_acc <- mean(preds$true == preds$pred, na.rm=TRUE)
    overall_auc <- AUC(preds$true, preds$prob)
    logv("Overall: acc=%.2f%% auc=%.4f", overall_acc * 100, overall_auc)
    list(
        models=models,
        predictions=preds,
        performance=perf,
        overall=list(acc=overall_acc, auc=overall_auc)
    )
}

# Experimental #####

#' @export
#' @rdname mdm
bootstrap_mdm <- function(
    x, y,
    decon_fun=deconvolute, align_fun=clupa,            feat_fun=peak_mat,
    fit_fun=fit_lasso,     predict_fun=predict_lasso,  mog=get_mog("default"),
    sfr=NULL,              igrs=list(),                use_rust=0,
    nworkers=1,            verbosity=1,                seed=1,
    test_frac=0.2,         max_rounds=30,              se_thresh=0.01
) {

    ns <- length(x); nr <- nrow(mog)
    nfolds_split <- ceiling(1 / test_frac)
    skip_decon <- identical(decon_fun, identity2)
    decon_fn <- if (identical(decon_fun, deconvolute)) deconvolute_spectra else decon_fun
    lvs <- levels(y); vb <- verbosity >= 2
    has_zero <- function(d) {
        inherits(d, "decons2") &&
            any(vapply(d, function(o) nrow(o$lcpar) == 0L, logical(1)))
    }
    ord <- with(mog, order(npmax, nfit, smit, smws, delta, maxShift, maxCombine))
    mog <- mog[ord, , drop=FALSE]; rownames(mog) <- NULL
    if (!skip_decon && any(mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=NULL, sfr=sfr, igrs=igrs,
            verbose=vb, nworkers=nworkers, use_rust=use_rust
        )
    }

    acc_mat <- matrix(NA_real_, nrow=0, ncol=nr)
    auc_mat <- matrix(NA_real_, nrow=0, ncol=nr)
    probs <- vector("list", nr)
    rnd <- 0L; se <- Inf

    while (rnd < max_rounds && se > se_thresh) {
        rnd <- rnd + 1L
        te <- get_test_ids(nfolds=nfolds_split, nsamples=ns, seed=seed + rnd, y=y)[[1]]
        tr <- setdiff(seq_len(ns), te)
        x_tr <- x[tr]; y_tr <- y[tr]; x_te <- x[te]; y_te <- y[te]
        logv("[bootstrap %d] %d train / %d test", rnd, length(tr), length(te))

        last_dkey <- NULL; last_akey <- NULL
        d_tr <- NULL; d_te <- NULL; a_tr <- NULL; a_te <- NULL
        a_tr_clupa <- NULL; a_te_clupa <- NULL; ref <- NULL
        cur_means_auc <- if (nrow(auc_mat) == 0) rep(-Inf, nr) else {
            m <- colMeans(auc_mat, na.rm=TRUE); ifelse(is.nan(m), -Inf, m)
        }
        rnd_acc <- rep(NA_real_, nr); rnd_auc <- rep(NA_real_, nr)

        for (i in seq_len(nr)) {
            r <- mog[i, , drop=FALSE]

            dkey <- list(r$npmax, r$nfit, r$smit, r$smws, r$delta)
            if (!identical(dkey, last_dkey)) {
                if (skip_decon) {
                    d_tr <- x_tr; d_te <- x_te
                } else {
                    d_tr <- decon_fn(
                        x=x_tr, sfr=sfr, igrs=igrs, verbose=vb, use_rust=use_rust,
                        nfit=r$nfit, smit=r$smit, smws=r$smws,
                        delta=r$delta, npmax=r$npmax, nworkers=nworkers
                    )
                    d_te <- decon_fn(
                        x=x_te, sfr=sfr, igrs=igrs, verbose=vb, use_rust=use_rust,
                        nfit=r$nfit, smit=r$smit, smws=r$smws,
                        delta=r$delta, npmax=r$npmax, nworkers=nworkers
                    )
                }
                last_dkey <- dkey; last_akey <- NULL
                if (has_zero(d_tr) || has_zero(d_te)) {
                    logv("[%d|%d/%d] zero peaks; skipping dkey group", rnd, i, nr)
                    d_tr <- NULL
                }
            }
            if (is.null(d_tr)) next

            akey <- list(dkey, r$maxShift)
            if (!identical(akey, last_akey)) {
                a_tr_clupa <- align_fun(x=d_tr, ref=NULL, maxShift=r$maxShift, verbose=vb, nworkers=nworkers, full=FALSE)
                ref <- if (inherits(a_tr_clupa, "aligns")) find_ref(a_tr_clupa) else NULL
                a_te_clupa <- align_fun(x=d_te, ref=ref, maxShift=r$maxShift, verbose=vb, nworkers=nworkers, full=FALSE)
                last_akey <- akey
            }
            # RefPA snap (cheap; not cached).
            if (inherits(a_tr_clupa, "aligns")) {
                a_tr <- snap_to_ref(a_tr_clupa, ref=ref, maxCombine=r$maxCombine)
                a_te <- snap_to_ref(a_te_clupa, ref=ref, maxCombine=r$maxCombine)
            } else {
                a_tr <- a_tr_clupa
                a_te <- a_te_clupa
            }

            X_tr_full <- feat_fun(a_tr, maxCombine=r$maxCombine, igrs=igrs)
            peakPos <- which(colSums(X_tr_full != 0) > 0)
            X_tr <- X_tr_full[, peakPos, drop=FALSE]
            X_te <- feat_fun(a_te, maxCombine=r$maxCombine, igrs=igrs)[, peakPos, drop=FALSE]
            model <- fit_fun(X_tr, y_tr, seed=seed + rnd, nworkers=nworkers)
            prob_te <- predict_fun(model, X_te)
            prob <- stats::setNames(prob_te, get_names(x_te))
            perf <- mdm_eval(y_te, prob, lvs)
            rnd_acc[i] <- perf$acc; rnd_auc[i] <- perf$auc
            probs[[i]] <- c(probs[[i]], list(prob))

            all_acc <- c(if (nrow(acc_mat) > 0) acc_mat[, i], perf$acc)
            all_auc <- c(if (nrow(auc_mat) > 0) auc_mat[, i], perf$auc)
            n <- sum(!is.na(all_auc))
            m_acc <- mean(all_acc, na.rm=TRUE); m_auc <- mean(all_auc, na.rm=TRUE)
            se_acc <- if (n >= 2) stats::sd(all_acc, na.rm=TRUE) / sqrt(n) else NA_real_
            se_auc <- if (n >= 2) stats::sd(all_auc, na.rm=TRUE) / sqrt(n) else NA_real_
            cur_means_auc[i] <- m_auc
            is_best <- !is.na(m_auc) && m_auc >= max(cur_means_auc, na.rm=TRUE)
            acc_s <- if (is.na(se_acc)) sprintf("%.2f%%", m_acc * 100) else sprintf("%.2f%% +- %.2f%%", m_acc * 100, se_acc * 100)
            auc_s <- if (is.na(se_auc)) sprintf("%.4f", m_auc) else sprintf("%.4f +- %.4f", m_auc, se_auc)
            logv("[%d|%d/%d] p=%d f=%d i=%d w=%d d=%g S=%d C=%g acc=%s auc=%s%s",
                 rnd, i, nr, r$npmax, r$nfit, r$smit, r$smws, r$delta,
                 r$maxShift, r$maxCombine, acc_s, auc_s,
                 if (is_best) " <-- BEST" else "")
        }
        acc_mat <- rbind(acc_mat, rnd_acc)
        auc_mat <- rbind(auc_mat, rnd_auc)

        if (rnd >= 3) {
            mean_aucs <- colMeans(auc_mat, na.rm=TRUE)
            bi <- which.max(mean_aucs)
            valid <- sum(!is.na(auc_mat[, bi]))
            if (valid >= 2) se <- stats::sd(auc_mat[, bi], na.rm=TRUE) / sqrt(valid)
        }
    }

    mean_accs <- colMeans(acc_mat, na.rm=TRUE)
    mean_aucs <- colMeans(auc_mat, na.rm=TRUE)
    auc_ses <- apply(auc_mat, 2, function(v) {
        n <- sum(!is.na(v))
        if (n < 2) NA_real_ else stats::sd(v, na.rm=TRUE) / sqrt(n)
    })
    bi <- which.max(mean_aucs)
    mog$acc <- mean_accs; mog$auc <- mean_aucs; mog$auc_se <- auc_ses
    logv("Best [%d/%d]: acc=%.2f%% auc=%.4f se=%.4f",
         bi, nr, mean_accs[bi] * 100, mean_aucs[bi], auc_ses[bi])
    list(best=mog[bi, ], mog=mog, probs=probs, acc_mat=acc_mat, auc_mat=auc_mat)
}

#' @export
#' @rdname mdm
cv_mdm <- function(
    x, y,
    decon_fun=deconvolute, align_fun=clupa,            feat_fun=peak_mat,
    fit_fun=fit_lasso,     predict_fun=predict_lasso,  mog=get_mog("default"),
    sfr=NULL,              igrs=list(),                use_rust=0,
    nworkers=1,            verbosity=1,                seed=1,
    nrounds=1,             kfold=10
) {

    ns <- length(x); nr <- nrow(mog)
    stopifnot(is_int(nrounds, 1), nrounds >= 1, is_int(kfold, 1), kfold >= 2, kfold <= ns)
    skip_decon <- identical(decon_fun, identity2)
    decon_fn <- if (identical(decon_fun, deconvolute)) deconvolute_spectra else decon_fun
    lvs <- levels(y); vb <- verbosity >= 2
    has_zero <- function(d) {
        inherits(d, "decons2") &&
            any(vapply(d, function(o) nrow(o$lcpar) == 0L, logical(1)))
    }
    ord <- with(mog, order(npmax, nfit, smit, smws, delta, maxShift, maxCombine))
    mog <- mog[ord, , drop=FALSE]; rownames(mog) <- NULL
    if (!skip_decon && any(mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=NULL, sfr=sfr, igrs=igrs,
            verbose=vb, nworkers=nworkers, use_rust=use_rust
        )
    }

    acc_mat <- matrix(NA_real_, nrow=nrounds, ncol=nr)
    auc_mat <- matrix(NA_real_, nrow=nrounds, ncol=nr)
    probs <- vector("list", nr)
    nm <- get_names(x)

    for (rnd in seq_len(nrounds)) {
        fold_ids <- get_test_ids(nfolds=kfold, nsamples=ns, seed=seed + rnd, y=y)
        pred_round <- matrix(NA_real_, nrow=ns, ncol=nr)
        logv("[cv %d/%d] %d-fold CV (ns=%d)", rnd, nrounds, kfold, ns)
        for (f in seq_along(fold_ids)) {
            te <- fold_ids[[f]]
            tr <- setdiff(seq_len(ns), te)
            x_tr <- x[tr]; y_tr <- y[tr]; x_te <- x[te]; y_te <- y[te]
            logvv("[cv %d/%d fold %d/%d] %d train / %d test", rnd, nrounds, f, kfold, length(tr), length(te))
            last_dkey <- NULL; last_akey <- NULL
            d_tr <- NULL; d_te <- NULL; a_tr <- NULL; a_te <- NULL
            a_tr_clupa <- NULL; a_te_clupa <- NULL; ref <- NULL
            for (i in seq_len(nr)) {
                r <- mog[i, , drop=FALSE]
                dkey <- list(r$npmax, r$nfit, r$smit, r$smws, r$delta)
                if (!identical(dkey, last_dkey)) {
                    if (skip_decon) {
                        d_tr <- x_tr; d_te <- x_te
                    } else {
                        d_tr <- decon_fn(
                            x=x_tr, sfr=sfr, igrs=igrs, verbose=vb, use_rust=use_rust,
                            nfit=r$nfit, smit=r$smit, smws=r$smws,
                            delta=r$delta, npmax=r$npmax, nworkers=nworkers
                        )
                        d_te <- decon_fn(
                            x=x_te, sfr=sfr, igrs=igrs, verbose=vb, use_rust=use_rust,
                            nfit=r$nfit, smit=r$smit, smws=r$smws,
                            delta=r$delta, npmax=r$npmax, nworkers=nworkers
                        )
                    }
                    last_dkey <- dkey; last_akey <- NULL
                    if (has_zero(d_tr) || has_zero(d_te)) {
                        logvv("[%d|%d/%d fold %d] zero peaks; skipping dkey group",
                              rnd, i, nr, f)
                        d_tr <- NULL
                    }
                }
                if (is.null(d_tr)) next

                akey <- list(dkey, r$maxShift)
                if (!identical(akey, last_akey)) {
                    a_tr_clupa <- align_fun(x=d_tr, ref=NULL, maxShift=r$maxShift, verbose=vb, nworkers=nworkers, full=FALSE)
                    ref <- if (inherits(a_tr_clupa, "aligns")) find_ref(a_tr_clupa) else NULL
                    a_te_clupa <- align_fun(x=d_te, ref=ref, maxShift=r$maxShift, verbose=vb, nworkers=nworkers, full=FALSE)
                    last_akey <- akey
                }
                # RefPA snap (cheap; not cached).
                if (inherits(a_tr_clupa, "aligns")) {
                    a_tr <- snap_to_ref(a_tr_clupa, ref=ref, maxCombine=r$maxCombine)
                    a_te <- snap_to_ref(a_te_clupa, ref=ref, maxCombine=r$maxCombine)
                } else {
                    a_tr <- a_tr_clupa
                    a_te <- a_te_clupa
                }

                X_tr_full <- feat_fun(a_tr, maxCombine=r$maxCombine, igrs=igrs)
                peakPos <- which(colSums(X_tr_full != 0) > 0)
                X_tr <- X_tr_full[, peakPos, drop=FALSE]
                X_te <- feat_fun(a_te, maxCombine=r$maxCombine, igrs=igrs)[, peakPos, drop=FALSE]
                model <- fit_fun(X_tr, y_tr, seed=seed + rnd, nworkers=nworkers)
                prob_te <- predict_fun(model, X_te)
                pred_round[te, i] <- prob_te
            }
        }

        for (i in seq_len(nr)) {
            if (anyNA(pred_round[, i])) next
            prob <- stats::setNames(pred_round[, i], nm)
            perf <- mdm_eval(y, prob, lvs)
            acc_mat[rnd, i] <- perf$acc; auc_mat[rnd, i] <- perf$auc
            probs[[i]] <- c(probs[[i]], list(prob))
        }
        running_means <- colMeans(auc_mat[seq_len(rnd), , drop=FALSE], na.rm=TRUE)
        running_means[is.nan(running_means)] <- -Inf
        bi <- which.max(running_means)
        for (i in seq_len(nr)) {
            if (is.na(auc_mat[rnd, i])) next
            all_acc <- acc_mat[seq_len(rnd), i]
            all_auc <- auc_mat[seq_len(rnd), i]
            n <- sum(!is.na(all_auc))
            m_acc <- mean(all_acc, na.rm=TRUE); m_auc <- mean(all_auc, na.rm=TRUE)
            se_acc <- if (n >= 2) stats::sd(all_acc, na.rm=TRUE) / sqrt(n) else NA_real_
            se_auc <- if (n >= 2) stats::sd(all_auc, na.rm=TRUE) / sqrt(n) else NA_real_
            acc_s <- if (is.na(se_acc)) sprintf("%.2f%%", m_acc * 100) else sprintf("%.2f%% +- %.2f%%", m_acc * 100, se_acc * 100)
            auc_s <- if (is.na(se_auc)) sprintf("%.4f", m_auc) else sprintf("%.4f +- %.4f", m_auc, se_auc)
            r <- mog[i, , drop=FALSE]
            logv("[%d|%d/%d] p=%d f=%d i=%d w=%d d=%g S=%d C=%g acc=%s auc=%s%s",
                 rnd, i, nr, r$npmax, r$nfit, r$smit, r$smws, r$delta,
                 r$maxShift, r$maxCombine, acc_s, auc_s,
                 if (i == bi) " <-- BEST" else "")
        }
    }

    mean_accs <- colMeans(acc_mat, na.rm=TRUE)
    mean_aucs <- colMeans(auc_mat, na.rm=TRUE)
    auc_ses <- apply(auc_mat, 2, function(v) {
        n <- sum(!is.na(v))
        if (n < 2) NA_real_ else stats::sd(v, na.rm=TRUE) / sqrt(n)
    })
    bi <- which.max(mean_aucs)
    mog$acc <- mean_accs; mog$auc <- mean_aucs; mog$auc_se <- auc_ses
    logv("Best [%d/%d]: acc=%.2f%% auc=%.4f se=%.4f",
         bi, nr, mean_accs[bi] * 100, mean_aucs[bi], auc_ses[bi])
    list(best=mog[bi, ], mog=mog, probs=probs, acc_mat=acc_mat, auc_mat=auc_mat)
}

# Helpers #####

#' @export
#' @rdname mdm
#' @title Identity decon function for fit_mdm
#' @description
#' No-op replacement for the `decon_fun` argument of
#' [metabodecon::fit_mdm()] / [metabodecon::benchmark()]. Returns its first
#' argument unchanged and ignores all other arguments. Use this to skip
#' the deconvolution stage of the pipeline (e.g. for binning baselines).
#' @return `x`, unchanged.
identity2 <- function(x, ...) x

as_binary01 <- function(y) {
    lvs <- sort(unique(y))
    if (length(lvs) != 2) stop("y must have exactly 2 unique levels")
    as.integer(y == lvs[2])
}

get_test_ids <- function(nfolds=5, nsamples, seed=1, y=NULL) {
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
    x, y, mog,
    sfr=NULL, igrs=list(), use_rust=NULL, nworkers=NULL,
    verbosity=NULL, seed=NULL, nfolds=NULL
) {
    cols <- c(
        "nfit", "smit", "smws", "delta", "npmax",
        "maxShift", "maxCombine"
    )
    stopifnot(
        is_spectra(x),
        is.factor(y),
        length(y) == length(x),
        is.data.frame(mog),
        nrow(mog) >= 1,
        all(cols %in% names(mog)),
        is_num_or_null(sfr, 2),
        is_list_of_nums(igrs, nv=2),
        is_bool_or_num(use_rust),
        is_int_or_null(nworkers, 1),
        is_int_or_null(verbosity, 1),
        is_int_or_null(seed, 1),
        is.null(nfolds) || (is_int(nfolds, 1) && nfolds >= 2)
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
    if (!is.null(nfolds) && nfolds > length(y)) {
        stop("`nfolds` must not exceed the number of samples.", call.=FALSE)
    }
    invisible(NULL)
}

check_bm_args <- function(x, y, igrs, nbin, seed, nfolds, verbosity) {
    stopifnot(
        is_spectra(x), is.factor(y), length(y) == length(x),
        is_list_of_nums(igrs, nv=2),
        is_int(nbin, 1), nbin >= 1,
        is_int_or_null(seed, 1),
        is_int(nfolds, 1), nfolds >= 2, nfolds <= length(y),
        is_int_or_null(verbosity, 1)
    )
    if (nlevels(y) != 2 || any(table(y) == 0)) {
        stop("`y` must contain exactly 2 non-empty classes.", call.=FALSE)
    }
    invisible(NULL)
}

mdm_eval <- function(y, prob, lvs) {
    pred <- factor(ifelse(prob > 0.5, lvs[2], lvs[1]), levels=lvs)
    list(acc=mean(pred == y), auc=AUC(y, prob))
}

# Format a (mean, SE) pair as "nn.n(±x.x)%" with both quantities
# expressed in percentage points. SE is omitted when NA.
fmt_pct_se <- function(m, s) {
    if (is.na(m)) return("NA")
    if (is.na(s)) return(sprintf("%.1f%%", m * 100))
    sprintf("%.1f(±%.1f)%%", m * 100, s * 100)
}

# Format a fraction as "nn.n%". NA -> "NA".
fmt_pct <- function(m) {
    if (is.na(m)) return("NA")
    sprintf("%.1f%%", m * 100)
}

#' @export
#' @rdname mdm
#' @title Lasso fitter for fit_mdm
#' @description
#' Fits an L1-penalised binomial logistic regression via
#' [glmnet::cv.glmnet()] with hardcoded internal `nfolds=10`. Returns the
#' fitted `cv.glmnet` object directly.
#' @param X Numeric feature matrix.
#' @param y Factor with two levels.
#' @param seed Random seed for inner-CV fold assignment.
#' @param nworkers When `> 1`, sets up a `doParallel` cluster and asks
#'   `cv.glmnet()` to parallelize its inner-CV folds via foreach. Falls
#'   back to single-threaded when `doParallel` isn't installed.
#' @return A `cv.glmnet` object (sub-classed `lasso`).
fit_lasso <- function(X, y, seed=1, nworkers=1L) {
    requireNamespace("glmnet", quietly=TRUE)
    set.seed(seed)
    par_ok <- nworkers > 1L && requireNamespace("doParallel", quietly=TRUE)
    if (par_ok) {
        cl <- parallel::makeCluster(min(nworkers, 10L))
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add=TRUE)
    }
    cvfit <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1, nfolds=10, parallel=par_ok)
    structure(cvfit, class=c("lasso", class(cvfit)))
}

#' @export
#' @rdname mdm
#' @title Lasso predictor for fit_mdm
#' @description
#' Companion of [metabodecon::fit_lasso()]. Returns the positive-class
#' probability at `lambda.min` for each row of `newx`.
#' @param model Object returned by [metabodecon::fit_lasso()].
#' @param newx Numeric feature matrix.
#' @return Numeric vector of length `nrow(newx)`.
predict_lasso <- function(model, newx) {
    requireNamespace("glmnet", quietly=TRUE)
    as.numeric(stats::predict(model, newx=newx, s="lambda.min", type="response"))
}

#' @export
#' @rdname mdm
#' @title Random-forest fitter for fit_mdm
#' @description
#' Fits a probability random forest with 1000 trees via
#' [ranger::ranger()]. Returns the fitted ranger object.
#' @param X Numeric feature matrix.
#' @param y Factor with two levels.
#' @param seed Random seed for ranger.
#' @param nworkers Forwarded to `ranger::ranger(num.threads=...)` for
#'   per-tree parallelism.
#' @return A `ranger` object (sub-classed `ranger500`). The trained
#'   levels of `y` are stored on the object for prediction-time use.
fit_ranger500 <- function(X, y, seed=1, nworkers=1L) {
    requireNamespace("ranger", quietly=TRUE)
    rf <- ranger::ranger(x=X, y=y, probability=TRUE, num.trees=1000, seed=seed,
                         num.threads=max(1L, nworkers))
    rf$lvs <- levels(y)
    structure(rf, class=c("ranger500", class(rf)))
}

#' @export
#' @rdname mdm
#' @title Random-forest predictor for fit_mdm
#' @description
#' Companion of [metabodecon::fit_ranger500()]. Returns the positive-class
#' probability for each row of `newx`.
#' @param model Object returned by [metabodecon::fit_ranger500()].
#' @param newx Numeric feature matrix.
#' @return Numeric vector of length `nrow(newx)`.
predict_ranger500 <- function(model, newx) {
    requireNamespace("ranger", quietly=TRUE)
    colnames(newx) <- model$forest$independent.variable.names
    pm <- stats::predict(model, data=newx)$predictions
    pm[, model$lvs[2]]
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
#'     model=NULL,
#'     ref=NULL,
#'     params=list(npmax=1000, nfit=3, smit=2, smws=5,
#'                 delta=6.4, maxShift=100, maxCombine=50)
#'   ),
#'   class="mdm"
#' )
#' print(m)
#' summary(m)
#'
#' \dontrun{
#'   m <- fit_mdm(spectra, y, mog=get_mog("default"), sfr=c(11, -2))
#'   predict(m, test_spectra, type="prob")
#'   coef(m)
#'   plot(m)
#' }
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
            nfit=p$nfit, smit=p$smit, smws=p$smws, delta=p$delta,
            npmax=p$npmax, nworkers=nworkers
        )
        logv("Aligning spectra with %d nworkers", nworkers)
        align_fun <- p$align_fun %||% clupa
        a <- align_fun(x=d, ref=object$ref, maxShift=p$maxShift,
                verbose=verbosity >= 2, nworkers=nworkers, full=FALSE)
        # Prediction always uses snap_to_ref: the training peakPos was
        # locked in once the model was fit, so new spectra must be
        # snapped to those exact reference columns regardless of which
        # post-align method was used during discovery.
        if (inherits(a, "aligns") && (p$maxCombine %||% 0L) > 0L) {
            a <- snap_to_ref(a, ref=object$ref, maxCombine=p$maxCombine)
        }
        Xn <- p$feat_fun(a, maxCombine=p$maxCombine, igrs=p$igrs %||% list())
        Xn <- Xn[, p$peakPos, drop=FALSE]
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
    pp <- c(
        "npmax", "nfit", "smit", "smws", "delta",
        "maxShift", "maxCombine"
    )
    cat("metabodecon model (mdm)\n")
    cat("  ", formatC("model:", width=-15), paste(class(x$model), collapse=", "),
        "\n", sep="")
    for (nm in pp) {
        v <- x$params[[nm]]
        if (is.null(v)) next
        lab <- formatC(paste0(nm, ":"), width=-15)
        cat("  ", lab, v, "\n", sep="")
    }
    if (!is.null(x$mog)) {
        cat("  grid rows:     ", nrow(x$mog), "\n", sep="")
    }
    invisible(x)
}

#' @export
#' @rdname mdm_methods
coef.mdm <- function(object, ...) {
    stopifnot(inherits(object, "mdm"), !is.null(object$model))
    if (inherits(object$model, "ranger500")) {
        return(object$model$variable.importance)
    }
    stats::coef(object$model, s="lambda.min", ...)
}

#' @export
#' @rdname mdm_methods
plot.mdm <- function(x, ...) {
    stopifnot(inherits(x, "mdm"), !is.null(x$model))
    if (inherits(x$model, "ranger500")) {
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
    pp <- c(
        "npmax", "nfit", "smit", "smws", "delta",
        "maxShift", "maxCombine"
    )
    out <- object$params[pp]
    out$model <- paste(class(object$model), collapse=", ")
    out$n_peaks <- length(object$params$peakPos %||% integer(0))
    out$grid_rows <- if (is.null(object$mog)) 0L else nrow(object$mog)
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
