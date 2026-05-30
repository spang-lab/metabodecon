
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
#' Utilities for fitting and benchmarking 'metabodecon models' (mdm).
#'
#' An `mdm` is a binary classification model fitted on a feature matrix
#' built from NMR spectra by the following pipeline:
#'
#' \preformatted{
#'   x  --decon_fun-->  d  --align_fun-->  a  --snap_fun-->  s
#'   s  --feat_fun-->  X  --fit_fun-->  (model, acc, auc)
#' }
#'
#' All five stages are pluggable. The `fit_fun` is responsible for
#' returning **both** a trained model and a stable estimate of the
#' model's generalization performance — OOB for [metabodecon::fit_ranger()],
#' repeated `cv.glmnet` OOF for [metabodecon::fit_lasso()]. That single
#' contract removes the outer-CV machinery from `fit_mdm` itself:
#' `fit_mdm` only sweeps preprocessing rows in `mog` and keeps the row
#' with the highest reported AUC.
#'
#' Pass `decon_fun=identity2`, `align_fun=identity_align` and/or
#' `snap_fun=identity_snap` to skip those stages and fit baselines (e.g.
#' a binning model) directly on raw spectra.
#'
#' [metabodecon::benchmark()] wraps [metabodecon::fit_mdm()] in an outer
#' k-fold cross-validation to estimate end-to-end predictive performance
#' on held-out spectra.
#'
#' @details
#'
#' ## Pluggable interfaces
#'
#' `decon_fun`, `align_fun`, `snap_fun`, `feat_fun`, `fit_fun` and
#' `predict_fun` speak a fixed parameter vocabulary. Replacements must
#' accept the listed arguments (extras via `...`).
#'
#' \describe{
#'   \item{`decon_fun(x, sfr, igrs, verbose, use_rust, nfit, smit, smws,
#'         delta, npmax, nworkers)`}{Returns a `decons2` object.
#'         [metabodecon::deconvolute()] (default), [metabodecon::identity2()].}
#'   \item{`align_fun(x, ref, maxShift, verbose, nworkers, full, ...)`}{Returns
#'         an `aligns` (or pass-through). [metabodecon::clupa()] (default),
#'         [metabodecon::identity_align()].}
#'   \item{`snap_fun(x, ref=NULL, maxCombine, ...)`}{Returns an `aligns`
#'         with the per-peak `pcisn` / `x0sn` columns populated.
#'         [metabodecon::snap_to_ref()] (default),
#'         [metabodecon::combine_peaks()], [metabodecon::snap_nw_blind()],
#'         [metabodecon::identity_snap()].}
#'   \item{`feat_fun(x, maxCombine, igrs, ...)`}{Returns a numeric matrix
#'         with one row per spectrum. [metabodecon::peak_mat()] (default),
#'         [metabodecon::si_mat()], [metabodecon::bin()].}
#'   \item{`fit_fun(X, y, seed, nworkers)`}{Returns a list with elements
#'         `model` (trained backend object), `acc` and `auc` (scalar
#'         generalization estimates in `[0, 1]`), and optionally `acc_se`
#'         / `auc_se`. Must call `requireNamespace("<backend>")` so it
#'         works after a fresh `readRDS()`. Built-ins:
#'         [metabodecon::fit_lasso()], [metabodecon::fit_ranger()].}
#'   \item{`predict_fun(model, newx)`}{Returns a numeric vector of
#'         positive-class probabilities. Must call
#'         `requireNamespace("<backend>")`. Built-ins:
#'         [metabodecon::predict_lasso()], [metabodecon::predict_ranger()].}
#' }
#'
#' ## Caching within `fit_mdm`
#'
#' Rows of `mog` are sorted by `(npmax, nfit, smit, smws, delta, maxShift,
#' maxCombine)` so identical decon-tuples and align-tuples cluster.
#' Inside the loop, the most recent deconvolution and alignment are kept
#' and reused whenever the relevant subset of parameters is unchanged.
#'
#' ## "auto" sentinels
#'
#' `npmax`, `maxShift` and `maxCombine` cells may be `NA` to request
#' automatic selection:
#' \itemize{
#'   \item `npmax=NA` ("auto") — resolved once via
#'         [metabodecon::find_npmax_elbow()] on the pre-attached `$deg`
#'         cache (median per-spectrum Kneedle elbow).
#'   \item `maxShift=NA` ("auto") — resolved per-`npmax` via
#'         `find_maxShift_dip()` (sweep CluPA shifts at powers of 2 and
#'         stop one step before the alignment-correlation dip). Requires
#'         `align_fun=clupa`.
#'   \item `maxCombine=NA` ("auto") — set to the row's resolved `maxShift`.
#' }
#' A negative `maxCombine` is also treated as a switch and replaced by
#' the row's `maxShift`.
#'
#' ## Grid-search results carried by spectra
#'
#' When any row of `mog` has `npmax > 0` (or `NA`), [metabodecon::fit_mdm()]
#' calls [metabodecon::grid_deconvolute_spectra()] once up front to
#' attach a `$deg` element to each spectrum. The enriched spectra are
#' reused across rows and outer folds. [metabodecon::benchmark()] does
#' the same up-front attachment.
#'
#' ## Parallelism
#'
#' [metabodecon::benchmark()] runs outer folds sequentially and delegates
#' all parallelism to the inner fitter via `nworkers`.
#'
#' @param x Spectra object.
#' @param y Factor vector with class labels for each spectrum.
#' @param mog Model-fitting grid as returned by [metabodecon::get_mog()].
#'   Missing columns are filled in with [metabodecon::deconvolute()] and
#'   [metabodecon::snap_to_ref()] defaults.
#' @param deg Deconvolution-parameter grid forwarded to
#'   [metabodecon::grid_deconvolute_spectra()] when any row of `mog` has
#'   `npmax > 0` (or `NA`). When `NULL` (default), the default grid built
#'   into [metabodecon::grid_deconvolute_spectra()] is used.
#' @param sfr Signal-free region. See [metabodecon::deconvolute()].
#' @param use_rust Use the Rust backend?
#' @param nworkers Number of workers for deconvolution, alignment and the
#'   inner fitter.
#' @param verbosity Verbosity level.
#' @param seed Random seed. Forwarded to `fit_fun`; also used for
#'   stratified fold assignment inside [metabodecon::benchmark()].
#' @param check Validate inputs at function entry?
#' @param decon_fun,align_fun,snap_fun,feat_fun,fit_fun,predict_fun
#'   See *Pluggable interfaces*.
#' @param igrs Ignore regions in ppm.
#' @param k Number of outer folds for [metabodecon::benchmark()].
#' @param conf Character string selecting a predefined `mog` configuration.
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with elements
#' `model` (best fitted backend model, refit on the full data by the
#' chosen `fit_fun`), `ref` (a list `list(align, snap)` carrying the
#' references needed to replay the pipeline at prediction time; `NULL`
#' when `align_fun=identity_align`), `params` (everything needed to
#' reproduce predictions: chosen grid row, the pluggable functions,
#' `lvs`, `peakPos`, `sfr`, `igrs`, `use_rust`, `snap_kind`, `snap_arg`),
#' and `mog` (input grid augmented with `acc`, `acc_se`, `auc` and
#' `auc_se` columns reported by `fit_fun`).
#'
#' [metabodecon::benchmark()] returns a list with elements:
#' - `models`: list of fitted models, one per outer fold.
#' - `predictions`: data frame with columns `fold`, `true`, `link`, `prob`,
#'   `pred`.
#' - `performance`: data frame with per-fold `acc` and `auc`.
#' - `overall`: list with pooled `acc` and `auc`.
#'
#' @examples
#' \dontrun{
#'   m <- fit_mdm(spectra, y, mog=get_mog("default"))
#'   bm <- benchmark(spectra, y, k=5, mog=get_mog("default"))
#'   mrf <- fit_mdm(spectra, y, fit_fun=fit_ranger,
#'                  predict_fun=predict_ranger)
#'   mnw <- fit_mdm(spectra, y, snap_fun=snap_nw_blind,
#'                  fit_fun=fit_ranger, predict_fun=predict_ranger)
#' }
#'
fit_mdm <- function(x, y,
    decon_fun=deconvolute, align_fun=clupa,           snap_fun=snap_to_ref,
    feat_fun=peak_mat,     fit_fun=fit_lasso,         predict_fun=predict_lasso,
    mog=get_mog("default"),
    deg=NULL,              sfr=NULL,                  igrs=list(),
    use_rust=0,            nworkers=1,                verbosity=1,
    seed=1,                check=TRUE
) {
    stopifnot(
        is.function(decon_fun), is.function(align_fun),
        is.function(snap_fun),  is.function(feat_fun),
        is.function(fit_fun),   is.function(predict_fun)
    )
    mog <- normalize_mog(mog)
    if (check) check_mdm_args(
        x=x, y=y, mog=mog, sfr=sfr, igrs=igrs,
        use_rust=use_rust, nworkers=nworkers, verbosity=verbosity,
        seed=seed
    )
    skip_decon <- identical(decon_fun, identity2)
    decon_fn <- if (identical(decon_fun, deconvolute)) deconvolute_spectra else decon_fun
    lvs <- levels(y)

    # Sort rows so identical decon/align tuples cluster.
    ord <- with(mog, order(npmax, nfit, smit, smws, delta, maxShift, maxCombine,
                            na.last=TRUE))
    mog <- mog[ord, , drop=FALSE]
    rownames(mog) <- NULL
    mog$acc <- NA_real_; mog$auc <- NA_real_
    mog$acc_se <- NA_real_; mog$auc_se <- NA_real_

    # Negative maxCombine acts as a switch (NA-safe).
    neg <- !is.na(mog$maxCombine) & mog$maxCombine < 0
    mog$maxCombine[neg] <- mog$maxShift[neg]

    # Pre-attach per-spectrum `$deg` tables when any row uses npmax > 0
    # or "auto" (NA).
    if (!skip_decon && any(is.na(mog$npmax) | mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2,
            nworkers=min(nworkers, length(x)), use_rust=use_rust
        )
    }
    # Resolve npmax="auto" -> median per-spectrum elbow.
    if (any(is.na(mog$npmax))) {
        auto_np <- find_npmax_elbow(x)
        logv("auto-pick npmax=%d (median elbow over %d spectra)",
             auto_np, length(x))
        mog$npmax[is.na(mog$npmax)] <- auto_np
    }

    nr <- nrow(mog); ns <- length(x)
    logv("Starting grid search (%d combinations, %d spectra)", nr, ns)
    last_dkey <- NULL; last_akey <- NULL
    d <- NULL; a_aligned <- NULL
    best_mdm <- NULL; best_auc <- -Inf
    auto_picks <- list()  # cache: npmax -> auto-selected maxShift
    for (i in seq_len(nr)) {
        r <- mog[i, , drop=FALSE]
        dkey <- list(r$npmax, r$nfit, r$smit, r$smws, r$delta)
        if (!identical(dkey, last_dkey)) {
            if (skip_decon) {
                d <- x
            } else {
                d <- decon_fn(
                    x=x, sfr=sfr, igrs=igrs, verbose=verbosity >= 2,
                    use_rust=use_rust, nfit=r$nfit, smit=r$smit, smws=r$smws,
                    delta=r$delta, npmax=r$npmax, nworkers=nworkers
                )
            }
            if (has_zero_peaks(d)) {
                logv("[c=%d/%d] zero peaks; skipping dkey group", i, nr)
                mog$acc[i] <- 0; mog$auc[i] <- 0
                last_dkey <- dkey; last_akey <- NULL
                next
            }
            last_dkey <- dkey; last_akey <- NULL
        }

        # Resolve maxShift="auto" rows. Cached per npmax.
        if (is.na(r$maxShift)) {
            if (!identical(align_fun, clupa)) stop(
                "maxShift='auto' is only supported with align_fun=clupa.",
                call.=FALSE
            )
            ck <- as.character(r$npmax)
            if (is.null(auto_picks[[ck]])) {
                pr <- find_maxShift_dip(d, nworkers=nworkers,
                                        verbose=verbosity >= 2)
                auto_picks[[ck]] <- pr$pick
                logv("auto-pick maxShift=%d for npmax=%d (ks=[%s] stopped=%s)",
                     pr$pick, r$npmax,
                     paste(pr$ks, collapse=","), pr$stopped)
            }
            r$maxShift <- auto_picks[[ck]]
            mog$maxShift[i] <- r$maxShift
            last_akey <- NULL
        }
        if (is.na(r$maxCombine)) {
            r$maxCombine <- r$maxShift
            mog$maxCombine[i] <- r$maxCombine
        }

        akey <- list(dkey, r$maxShift)
        if (!identical(akey, last_akey)) {
            a_aligned <- align_fun(x=d, ref=NULL, maxShift=r$maxShift,
                                    verbose=verbosity >= 2,
                                    nworkers=nworkers, full=FALSE)
            last_akey <- akey
        }

        a_snap <- snap_fun(a_aligned, ref=NULL, maxCombine=r$maxCombine,
                            igrs=igrs)

        # Extract refs to store for predict-time replay.
        snap_kind <- attr(a_snap, "snap_kind") %||% "ref"
        snap_arg  <- attr(a_snap, "snap_arg")  %||% r$maxCombine
        snap_ref  <- attr(a_snap, "ref")
        align_ref <- if (inherits(a_aligned, "aligns")) find_ref(a_aligned) else NULL
        snap_ref  <- snap_ref %||% align_ref
        refs <- if (is.null(align_ref) && is.null(snap_ref)) NULL
                else list(align=align_ref, snap=snap_ref)

        X_full <- feat_fun(a_snap, maxCombine=r$maxCombine, igrs=igrs)
        pp <- which(colSums(X_full != 0) > 0)
        if (length(pp) == 0L) next
        X <- X_full[, pp, drop=FALSE]
        res <- fit_fun(X, y, seed=seed, nworkers=nworkers)
        mog$acc[i] <- res$acc; mog$auc[i] <- res$auc
        mog$acc_se[i] <- res$acc_se %||% NA_real_
        mog$auc_se[i] <- res$auc_se %||% NA_real_

        is_best <- !is.na(res$auc) && res$auc > best_auc
        sym <- if (is_best) " <-- BEST" else ""
        logv("c=%d/%d p=%d S=%d C=%g acc=%s auc=%s%s",
             i, nr, r$npmax, r$maxShift, r$maxCombine,
             fmt_pct_se(res$acc, mog$acc_se[i]),
             fmt_pct_se(res$auc, mog$auc_se[i]), sym)
        if (!is_best) next
        best_auc <- res$auc
        params <- list(
            feat_fun=feat_fun, fit_fun=fit_fun, predict_fun=predict_fun,
            decon_fun=decon_fun, align_fun=align_fun, snap_fun=snap_fun,
            lvs=lvs, snap_kind=snap_kind, snap_arg=snap_arg,
            sfr=sfr, igrs=igrs, use_rust=use_rust,
            npmax=r$npmax, nfit=r$nfit, smit=r$smit, smws=r$smws,
            delta=r$delta, maxShift=r$maxShift, maxCombine=r$maxCombine,
            peakPos=pp
        )
        best_mdm <- structure(
            list(model=res$model, ref=refs, params=params), class="mdm"
        )
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
    decon_fun=deconvolute, align_fun=clupa,           snap_fun=snap_to_ref,
    feat_fun=peak_mat,     fit_fun=fit_lasso,         predict_fun=predict_lasso,
    mog=get_mog("default"),
    deg=NULL,              sfr=NULL,                  igrs=list(),
    use_rust=0,            nworkers=1,                verbosity=2,
    seed=1,                k=3,                       check=TRUE
) {
    stopifnot(
        is.function(decon_fun), is.function(align_fun),
        is.function(snap_fun),  is.function(feat_fun),
        is.function(fit_fun),   is.function(predict_fun),
        is_int(k, 1), k >= 2, k <= length(y)
    )
    mog <- normalize_mog(mog)
    if (check) check_mdm_args(
        x=x, y=y, mog=mog, sfr=sfr, igrs=igrs,
        use_rust=use_rust, nworkers=nworkers, verbosity=verbosity, seed=seed
    )

    # One-time grid attach when any row of mog uses npmax > 0 or NA.
    if (!identical(decon_fun, identity2) &&
        any(is.na(mog$npmax) | mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=deg, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2, nworkers=nworkers, use_rust=use_rust
        )
    }

    inner_v <- max(0L, verbosity - 1L)
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
            decon_fun=decon_fun, align_fun=align_fun, snap_fun=snap_fun,
            feat_fun=feat_fun, fit_fun=fit_fun, predict_fun=predict_fun,
            mog=mog, deg=deg, sfr=sfr, igrs=igrs,
            use_rust=use_rust, nworkers=nworkers, verbosity=inner_v,
            seed=seed, check=FALSE
        )
        p <- stats::predict(m, x[te], type="all", nworkers=nworkers,
                            verbosity=inner_v)
        fp <- data.frame(fold=i, true=y[te], link=p$link, prob=p$prob,
                         pred=p$class)
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

#' @export
#' @rdname mdm
#' @title Identity decon function for fit_mdm
#' @description
#' No-op replacement for the `decon_fun` argument of
#' [metabodecon::fit_mdm()] / [metabodecon::benchmark()]. Returns its
#' first argument unchanged and ignores all other arguments. Use this to
#' skip the deconvolution stage of the pipeline (e.g. for binning
#' baselines).
#' @return `x`, unchanged.
identity2 <- function(x, ...) x

#' @export
#' @rdname mdm
#' @title Identity snap function for fit_mdm
#' @description
#' No-op replacement for the `snap_fun` argument of
#' [metabodecon::fit_mdm()]. Returns its first argument unchanged so the
#' feature-matrix stage operates on whatever the alignment stage
#' produced (without an extra snap-to-reference step).
#' @return `x`, unchanged.
identity_snap <- function(x, ref=NULL, maxCombine=0L, ...) x

#' @export
#' @rdname mdm
#'
#' @title Label-blind Needleman-Wunsch snap for fit_mdm
#'
#' @description
#' Thin wrapper around [metabodecon::build_consensus()] +
#' [metabodecon::snap_nw()] that drops into the `snap_fun` slot of
#' [metabodecon::fit_mdm()]. The consensus is built label-blind
#' (`y = NULL`) so cross-method comparisons (NW vs.
#' [metabodecon::snap_to_ref()]) and the OOB / repeated-CV scores
#' reported by `fit_fun` stay honest.
#'
#' The `maxCombine` window (in shared-grid columns) is translated into a
#' ppm `gap_tol` via the median spacing of `x[[1]]$cssh`. The translated
#' `gap_tol` is stashed on `attr(out, "snap_arg")` so [predict.mdm] can
#' replay the same NW snap at prediction time, and the consensus is
#' stashed on `attr(out, "ref")` for the same purpose. `attr(out,
#' "snap_kind") = "nw"` tells [predict.mdm] which branch to take.
#'
#' @param x A `decons2` / `aligns` object.
#' @param ref Optional pre-built consensus (used at predict time). When
#'   `NULL`, [metabodecon::build_consensus()] is called on `x` with
#'   `y=NULL`.
#' @param maxCombine Snap window in shared-grid columns. Translated to
#'   `gap_tol = maxCombine * median(diff(cssh))` (in ppm).
#' @param ... Ignored (signature compatibility).
#'
#' @return An `aligns` object with `pcisn` / `x0sn` populated, plus
#'   the attributes `snap_kind="nw"`, `snap_arg=<gap_tol>` and
#'   `ref=<consensus>`.
#'
snap_nw_blind <- function(x, ref=NULL, maxCombine=20, w_A=0, ...) {
    stopifnot(inherits(x, "decons2"))
    x <- ensure_cssh(x)
    cssh <- x[[1]]$cssh
    spacing <- abs(stats::median(diff(cssh)))
    gap_tol <- max(spacing, as.numeric(maxCombine) * spacing)
    pos_field <- if (!is.null(x[[1]]$lcpar$x0al)) "x0al" else "x0"
    if (is.null(ref)) {
        ref <- build_consensus(x, y=NULL, gap_tol=gap_tol, pos_field=pos_field)
    }
    out <- snap_nw(x, ref=ref, gap_tol=gap_tol, pos_field=pos_field, w_A=w_A)
    attr(out, "snap_kind") <- "nw"
    attr(out, "snap_arg")  <- gap_tol
    attr(out, "ref")       <- ref
    attr(out, "w_A")       <- w_A
    out
}

#' @export
#' @rdname mdm
#' @title Lasso fitter for fit_mdm
#' @description
#' Fits an L1-penalised binomial logistic regression via repeated
#' [glmnet::cv.glmnet()] (default `nreps=5`, internal `nfolds=10`,
#' `keep=TRUE`). All reps share the lambda path discovered by the first
#' rep so per-rep out-of-fold (OOF) predictions are directly comparable
#' lambda-by-lambda. For each lambda the per-rep OOF accuracy and AUC
#' are averaged across reps; the lambda that maximizes the averaged AUC
#' (`lambda*`) is the one reported (and the one used at predict time).
#' This averages the per-lambda performance curve *before* optimizing
#' over lambda, so the chosen lambda is stable across reps and the
#' reported acc/AUC reflects model variance rather than the noise of
#' lambda-pick instability across reps. `acc_se` / `auc_se` are the SE
#' across reps at `lambda*`. The model object is the last rep's
#' `cv.glmnet`; its `lambda.min` is overwritten with `lambda*` so
#' [metabodecon::predict_lasso()] picks the right lambda by default.
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
#'   has been overwritten with the AUC-maximizing `lambda*`), `acc`,
#'   `auc`, `acc_se`, `auc_se`.
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
    # Pick lambda* by averaged AUC; ties broken by larger lambda
    # (which.max returns the first index → smallest lambda in cv.glmnet's
    # decreasing path; flipping the search keeps the more regularized
    # solution among ties, matching cv.glmnet's lambda.1se sensibility).
    j_star <- which.max(rev(mean_auc))
    j_star <- nl - j_star + 1L
    chosen_lambda <- lambda_path[j_star]
    final_model <- cvs[[nreps]]
    final_model$lambda.min <- chosen_lambda
    list(
        model = final_model,
        acc = mean_acc[j_star],
        auc = mean_auc[j_star],
        acc_se = if (nreps >= 2L)
            stats::sd(acc_mat[, j_star], na.rm=TRUE) / sqrt(nreps) else NA_real_,
        auc_se = if (nreps >= 2L)
            stats::sd(auc_mat[, j_star], na.rm=TRUE) / sqrt(nreps) else NA_real_
    )
}

#' @export
#' @rdname mdm
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

#' @export
#' @rdname mdm
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
#' @return A list with `model` (a `ranger` object with the trained
#'   levels stashed on `model$lvs` for the predict path), `acc`, `auc`.
#'   `acc_se` and `auc_se` are `NA` because OOB produces a single point
#'   estimate.
fit_ranger <- function(X, y, seed=1, nworkers=1L, num.trees=5000L) {
    requireNamespace("ranger", quietly=TRUE)
    stopifnot(is_int(num.trees, 1), num.trees >= 1L)
    lvs <- levels(y)
    rf <- ranger::ranger(x=X, y=y, probability=TRUE, num.trees=num.trees,
                          seed=seed, num.threads=max(1L, nworkers))
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

#' @export
#' @rdname mdm
#' @title Random-forest predictor for fit_mdm
#' @description
#' Companion of [metabodecon::fit_ranger()]. Returns the positive-class
#' probability for each row of `newx`.
#' @param model Object returned in the `model` slot of
#'   [metabodecon::fit_ranger()].
#' @param newx Numeric feature matrix.
#' @return Numeric vector of length `nrow(newx)`.
predict_ranger <- function(model, newx) {
    requireNamespace("ranger", quietly=TRUE)
    colnames(newx) <- model$forest$independent.variable.names
    pm <- stats::predict(model, data=newx)$predictions
    pm[, model$lvs[2]]
}

# Helpers #####

# Ensure `mog` has every column fit_mdm expects, filling in
# deconvolute() / snap_to_ref() defaults for missing ones. Also coerces
# "auto" strings in npmax/maxShift/maxCombine to NA so the downstream
# sentinel handling is uniform.
normalize_mog <- function(mog) {
    stopifnot(is.data.frame(mog), nrow(mog) >= 1L)
    fill <- list(
        nfit=3L, smit=2L, smws=5L, delta=6.4,
        npmax=0L, maxShift=50L, maxCombine=20L
    )
    for (nm in names(fill)) {
        if (is.null(mog[[nm]])) mog[[nm]] <- fill[[nm]]
    }
    for (nm in c("npmax", "maxShift", "maxCombine")) {
        v <- mog[[nm]]
        if (is.character(v)) mog[[nm]] <- parse_int_with_auto(v, nm)
    }
    mog
}

# Convert an integer-valued argument that may contain the string "auto"
# into an integer vector with NA marking the "auto" entries. Numeric
# input is coerced to integer. Anything else triggers an error.
parse_int_with_auto <- function(x, name) {
    if (is.character(x)) {
        ok <- x == "auto" | !is.na(suppressWarnings(as.integer(x)))
        if (!all(ok)) {
            stop(sprintf(
                "%s entries must be non-negative integers or 'auto'.", name
            ), call.=FALSE)
        }
        out <- suppressWarnings(as.integer(x))
        out[x == "auto"] <- NA_integer_
        return(out)
    }
    as.integer(x)
}

# Per-spectrum npmax elbow from the (np, cum-min ar) frontier of `s$deg`.
# Same Kneedle-on-cum-min-frontier idea as `mdp::find_ellbow` but returns
# only the np at the knee. `s$deg` must be populated upstream (typically
# by grid_deconvolute_spectra()).
find_npmax_elbow_one <- function(s) {
    d <- s$deg
    if (is.null(d) || nrow(d) == 0L) return(NA_integer_)
    by_np <- split(seq_len(nrow(d)), d$np)
    idx <- vapply(by_np, function(ii) ii[which.min(d$ar[ii])], integer(1))
    f <- d[idx, , drop=FALSE]
    f <- f[order(f$np), , drop=FALSE]
    f$cum_ar <- cummin(f$ar)
    np_rng <- diff(range(f$np))
    ar_rng <- diff(range(f$cum_ar))
    if (np_rng == 0 || ar_rng == 0) return(as.integer(f$np[1]))
    nn <- (f$np - min(f$np)) / np_rng
    yn <- (f$cum_ar - min(f$cum_ar)) / ar_rng
    k <- which.max((1 - nn) - yn)
    as.integer(f$np[k])
}

# Aggregate per-spectrum elbows into a single npmax via the median.
# Requires every `x[[i]]` to carry a non-empty `$deg` grid.
find_npmax_elbow <- function(x) {
    picks <- vapply(x, find_npmax_elbow_one, integer(1))
    picks <- picks[!is.na(picks)]
    if (length(picks) == 0L) {
        stop("find_npmax_elbow: no spectra have a $deg grid; ",
             "call grid_deconvolute_spectra() first.", call.=FALSE)
    }
    as.integer(stats::median(picks))
}

# Adaptive maxShift selection by dip detection. Sweeps maxShift through
# {1, 2, 4, 8, ...}, runs CluPA at each step, computes the average
# pairwise Pearson correlation of the aligned superpositions
# (`sit$supal`), and stops the FIRST time the correlation decreases
# compared to the previous step. Returns the maxShift from the step
# *before* the dip (the last one that was still improving). If no dip is
# seen by `max_cap`, returns `max_cap`. Always uses [metabodecon::clupa()]
# because the dip metric reads `sit$supal`.
find_maxShift_dip <- function(d, max_cap=512L, nworkers=1, verbose=FALSE) {
    avg_pearson <- function(a) {
        M <- do.call(cbind, lapply(a, function(s) s$sit$supal))
        C <- stats::cor(M)
        mean(C[upper.tri(C)])
    }
    ks <- integer(0); ps <- numeric(0); ms <- 1L
    repeat {
        a <- clupa(x=d, ref=NULL, maxShift=as.integer(ms),
                   verbose=verbose, nworkers=nworkers, full=TRUE)
        p <- avg_pearson(a)
        ks <- c(ks, ms); ps <- c(ps, p)
        if (length(ps) >= 2L && p < ps[length(ps) - 1L]) {
            return(list(pick=ks[length(ks) - 1L], ks=ks, ps=ps,
                        stopped="dip"))
        }
        if (ms >= max_cap)
            return(list(pick=ms, ks=ks, ps=ps, stopped="cap"))
        ms <- ms * 2L
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
    x, y, mog,
    sfr=NULL, igrs=list(), use_rust=NULL, nworkers=NULL,
    verbosity=NULL, seed=NULL
) {
    cols <- c("nfit", "smit", "smws", "delta", "npmax",
              "maxShift", "maxCombine")
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
        is_int_or_null(seed, 1)
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
    sprintf("%.1f(±%.1f)%%", m * 100, s * 100)
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
            nfit=p$nfit, smit=p$smit, smws=p$smws, delta=p$delta,
            npmax=p$npmax, nworkers=nworkers
        )
        align_fun <- p$align_fun %||% clupa
        a_aligned <- align_fun(
            x=d, ref=object$ref$align, maxShift=p$maxShift,
            verbose=verbosity >= 2, nworkers=nworkers, full=FALSE
        )
        if (identical(p$snap_kind, "nw")) {
            pos_field <- if (!is.null(a_aligned[[1]]$lcpar$x0al)) "x0al" else "x0"
            a <- snap_nw(a_aligned, ref=object$ref$snap,
                          gap_tol=p$snap_arg, pos_field=pos_field)
        } else if (inherits(a_aligned, "aligns") && (p$maxCombine %||% 0L) > 0L) {
            a <- snap_to_ref(a_aligned, ref=object$ref$snap,
                              maxCombine=p$maxCombine)
        } else {
            a <- a_aligned
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
    pp <- c("npmax", "nfit", "smit", "smws", "delta",
            "maxShift", "maxCombine", "snap_kind")
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
    pp <- c("npmax", "nfit", "smit", "smws", "delta",
            "maxShift", "maxCombine", "snap_kind")
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
