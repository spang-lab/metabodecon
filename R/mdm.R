
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
#' Utilities for fitting, tuning and benchmarking 'metabodecon models' (mdm)
#' and 'binning models' (bm).
#'
#' A mdm is essentially a [glmnet::cv.glmnet()] lasso model, fitted on a
#' feature matrix obtained by deconvoluting and aligning spectra and snapping
#' their peaks to a shared reference. Deconvolution parameters
#' (`npmax`/`nfit`/`smit`/`smws`/`delta`), alignment parameter `maxShift` and
#' peak-combining parameter `maxCombine` are tunable hyperparameters.
#'
#' A bm is a lasso model fitted on a binned-intensity feature matrix and
#' serves as a simple baseline for comparison against `mdm` models.
#'
#' [metabodecon::fit_mdm()] iterates over the rows of a 'model fitting grid'
#' (`mog`) — see [metabodecon::get_mog()] — and for each row deconvolutes,
#' aligns, builds the feature matrix and fits a [glmnet::cv.glmnet()] model.
#' Held-out accuracy and AUC at `lambda.min` are recorded for each row. The
#' best model and the augmented grid are returned. A single-row `mog`
#' degenerates to fitting one model.
#'
#' [metabodecon::fit_bm()] fits a single binning-based lasso model.
#'
#' [metabodecon::benchmark()] runs an outer k-fold cross-validation over any
#' fitter function (`fit_mdm`, `fit_bm`, ...) to estimate end-to-end
#' predictive performance, and returns per-fold models, predictions and
#' performance metrics.
#'
#' @details
#'
#' ## Caching within `fit_mdm`
#'
#' Rows of `mog` are sorted by `(npmax, nfit, smit, smws, delta, maxShift,
#' maxCombine)` so identical decon-tuples and align-tuples cluster. Inside
#' the loop, the most recent deconvolution and alignment are kept and reused
#' whenever the relevant subset of parameters is unchanged.
#'
#' ## Grid-search results carried by spectra
#'
#' When any row of `mog` has `npmax > 0`, [metabodecon::fit_mdm()] calls
#' [metabodecon::grid_deconvolute_spectra()] once up front to attach a
#' `$deg` element to each spectrum. The enriched spectra are reused across
#' rows and outer folds. [metabodecon::benchmark()] does the same up-front
#' attachment for `fun = "fit_mdm"`.
#'
#' ## Parallelism
#'
#' [metabodecon::benchmark()] runs outer folds sequentially and delegates
#' all parallelism to the inner fitter via `nworkers` (`fit_mdm`).
#'
#' @param x Spectra object. May already carry per-spectrum `$deg` tables.
#' @param y Factor vector with class labels for each spectrum.
#' @param mog Model-fitting grid as returned by [metabodecon::get_mog()].
#' @param sfr Signal-free region. See [metabodecon::deconvolute()].
#' @param use_rust Use the Rust backend?
#' @param nworkers Number of workers for deconvolution and alignment.
#' @param verbosity Verbosity level.
#' @param seed Random seed for fold assignments.
#' @param nfolds Number of folds for the inner [glmnet::cv.glmnet()] call.
#' @param check Validate inputs at function entry?
#' @param igrs Ignore regions passed to `fun`.
#' @param nbin Number of bins in the non-ignored part of the ppm range.
#' @param fun Name of fitter function. Either "fit_mdm" or "fit_bm"`.
#' @param k Number of outer folds for [metabodecon::benchmark()].
#' @param ... Extra arguments forwarded to `fun`.
#' @param conf Character string selecting a predefined `mog` configuration.
#'
#' @return
#' [metabodecon::fit_mdm()] returns an object of class `mdm` with elements
#' `model` (best [glmnet::cv.glmnet()]), `ref` (reference alignment spectrum),
#' `params` (all settings needed to reproduce predictions: chosen grid row
#' plus non-grid arguments such as `sfr`, `igrs`, `use_rust` and `peakPos`)
#' and `mog` (input grid augmented with `acc`/`auc` columns).
#'
#' [metabodecon::fit_bm()] returns an object of class `bm` with elements
#' `model` and `params`.
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
#'   bm <- benchmark(spectra, y, fun="fit_mdm", k=5, mog=get_mog("default"))
#' }
fit_mdm <- function(
    x, y, mog=get_mog("default"),
    sfr=NULL, igrs=list(), use_rust=0.5, nworkers=1, verbosity=2,
    seed=1, nfolds=10, check=TRUE
) {
    if (check) check_mdm_args(
        x=x, y=y, mog=mog, sfr=sfr, igrs=igrs,
        use_rust=use_rust, nworkers=nworkers, verbosity=verbosity,
        seed=seed, nfolds=nfolds
    )

    # Sort rows so identical decon/align tuples cluster.
    ord <- with(
        mog, order(npmax, nfit, smit, smws, delta, maxShift, maxCombine)
    )
    mog <- mog[ord, , drop=FALSE]
    rownames(mog) <- NULL
    mog$acc <- NA_real_
    mog$auc <- NA_real_

    # Pre-attach per-spectrum `$deg` tables when any row uses npmax > 0.
    if (any(mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=mog, sfr=sfr, igrs=igrs,
            verbose=verbosity >= 2,
            nworkers=min(nworkers, length(x)), use_rust=use_rust
        )
    }

    nr <- nrow(mog)
    ns <- length(x)
    logv("Starting grid search (%d combinations, %d spectra)", nr, ns)
    foldid <- get_foldid(y=y, nfolds=nfolds, seed=seed)
    last_dkey <- NULL; last_akey <- NULL
    best_mdm <- NULL;  best_auc <- -Inf
    row_fmt <- "[%d/%d] npmax=%d nfit=%d smit=%d smws=%d delta=%g maxShift=%d maxCombine=%g"
    for (i in seq_len(nr)) {
        r <- mog[i, , drop=FALSE]
        logv(row_fmt, i, nr, r$npmax, r$nfit, r$smit, r$smws, r$delta, r$maxShift, r$maxCombine)
        dkey <- list(r$npmax, r$nfit, r$smit, r$smws, r$delta)
        if (!identical(dkey, last_dkey)) {
            d <- deconvolute_spectra(
                x=x, sfr=sfr, igrs=igrs, verbose=verbosity>=2,
                use_rust=use_rust, nfit=r$nfit, smit=r$smit, smws=r$smws,
                delta=r$delta, npmax=r$npmax, nworkers=nworkers, full=FALSE
            )
            last_dkey <- dkey
            last_akey <- NULL
        }
        nps <- vapply(d, function(o) nrow(o$lcpar), integer(1))
        if (any(nps == 0)) {
            logv("[%d/%d] %d spectra produced zero peaks; skipping", i, nr, sum(nps == 0))
            mog$acc[i] <- 0; mog$auc[i] <- 0
            next
        }

        akey <- list(dkey, r$maxShift)
        if (!identical(akey, last_akey)) {
            a <- align_decons(
                x=d, maxShift=r$maxShift, verbose=verbosity >= 2,
                nworkers=nworkers, full=FALSE
            )
            last_akey <- akey
        }

        ref <- find_ref(a)
        mat <- si_mat(a, maxCombine=r$maxCombine)
        peakPos <- which(colSums(mat != 0) > 0)
        X <- mat[, peakPos, drop=FALSE]

        cvfit <- glmnet::cv.glmnet(X, y, family="binomial", alpha=1, foldid=foldid, keep=TRUE)
        li <- which(cvfit$lambda == cvfit$lambda.min)
        link <- cvfit$fit.preval[, li]
        prob <- 1 / (1 + exp(-link))
        lvs <- levels(y)
        pred <- factor(ifelse(prob > 0.5, lvs[2], lvs[1]), levels=lvs)
        mog$acc[i] <- mean(pred == y)
        mog$auc[i] <- AUC(y, prob)
        logv("[%d/%d] acc=%.2f%% auc=%.4f", i, nr, mog$acc[i] * 100, mog$auc[i])

        if (!is.na(mog$auc[i]) && mog$auc[i] > best_auc) {
            best_auc <- mog$auc[i]
            params <- list(
                sfr=sfr, igrs=igrs, use_rust=use_rust, npmax=r$npmax,
                nfit=r$nfit, smit=r$smit, smws=r$smws, delta=r$delta,
                maxShift=r$maxShift, maxCombine=r$maxCombine, peakPos=peakPos
            )
            best_mdm <- structure(list(model=cvfit, ref=ref, params=params), class="mdm")
        }
    }

    best_mdm$mog <- mog
    ibest <- which.max(mog$auc)
    fmt <- "Best [%d/%d]: acc=%.2f%% auc=%.4f"
    logv(fmt, ibest, nr, mog$acc[ibest] * 100, mog$auc[ibest])
    best_mdm
}

#' @export
#' @rdname mdm
get_mog <- function(conf="default") {
    g <- expand.grid2(
        nfit = switch(conf, dynamic=0, 5),
        smit = switch(conf, dynamic=0, 2),
        smws = switch(conf, dynamic=0, static=c(3, 5, 7, 9), 5),
        delta = switch(conf, dynamic=0, static=c(3.2, 4.8, 6.4, 8.0), 6.4),
        npmax = switch(conf, dynamic=seq(400,1600,200), 0),
        maxShift = switch(conf, default=50, c(50,100,150,200,250)),
        maxCombine = switch(conf, default=5, c(5,10,20,30,40,50))
    )
    ord <- order(g$npmax, g$nfit, g$smit, g$smws, g$delta, g$maxShift, g$maxCombine)
    g <- g[ord, , drop=FALSE]
    rownames(g) <- NULL
    g
}

#' @export
#' @rdname mdm
fit_bm <- function(
    x, y, igrs=list(), nbin=1000,
    seed=1, nfolds=10, verbosity=2, check=TRUE
) {
    if (check) check_bm_args(
        x=x, y=y, igrs=igrs, nbin=nbin,
        seed=seed, nfolds=nfolds, verbosity=verbosity
    )
    logv(
        "Binning %d spectra into %d bins (igrs excluded)",
        length(x), nbin
    )
    X <- bin_spectra(x, igrs=igrs, nbin=nbin)
    foldid <- get_foldid(y=y, nfolds=nfolds, seed=seed)
    logv("Fitting cv.glmnet on %d features", ncol(X))
    model <- glmnet::cv.glmnet(
        X, y, family="binomial", alpha=1, foldid=foldid, keep=TRUE
    )
    params <- list(igrs=igrs, nbin=nbin, feat_names=colnames(X))
    structure(list(model=model, params=params), class="bm")
}

#' @export
#' @rdname mdm
benchmark <- function(
    x, y, ..., fun="fit_mdm",
    k=5, seed=1, verbosity=2
) {
    stopifnot(
        is_spectra(x), is.factor(y), length(y) == length(x),
        is_str(fun), is_int(k, 1), k >= 2,
        is_int(seed, 1), is_int(verbosity, 1)
    )
    if (nlevels(y) != 2 || any(table(y) == 0)) {
        stop("`y` must contain exactly 2 non-empty classes.", call.=FALSE)
    }
    if (k > length(y)) {
        stop("`k` must not exceed the number of samples.", call.=FALSE)
    }
    if (!fun %in% c("fit_mdm", "fit_bm")) {
        stop("Unsupported fun=", fun, call.=FALSE)
    }

    dots <- list(...)
    fitter <- match.fun(fun)
    pred_fn <- if (fun == "fit_mdm") predict.mdm else predict.bm

    # One-time grid attach when fitting mdm with npmax > 0.
    if (fun == "fit_mdm" && !is.null(dots$mog) && any(dots$mog$npmax > 0)) {
        x <- grid_deconvolute_spectra(
            x=x, deg=dots$mog, sfr=dots$sfr,
            igrs=dots$igrs %||% list(),
            verbose=verbosity >= 2,
            nworkers=dots$nworkers %||% 1L,
            use_rust=dots$use_rust %||% FALSE
        )
    }

    # Forward verbosity-1 to the fitter (overrides any user-passed value).
    dots$verbosity <- max(0L, verbosity - 1L)

    te_list <- get_test_ids(
        nfolds=k, nsamples=length(x), seed=seed, y=y
    )
    models <- vector("list", k)
    fold_preds <- vector("list", k)
    perf <- data.frame(
        fold=integer(0), acc=numeric(0), auc=numeric(0)
    )

    logv("Running %d-fold outer CV with fun=%s", k, fun)
    for (i in seq_along(te_list)) {
        te <- te_list[[i]]
        tr <- setdiff(seq_along(x), te)
        logv("[fold %d/%d] fitting", i, k)
        m <- do.call(fitter, c(list(x=x[tr], y=y[tr]), dots))
        p <- pred_fn(m, x[te], type="all")
        fp <- data.frame(
            fold=i, true=y[te], link=p$link,
            prob=p$prob, pred=p$class
        )
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
    logv(
        "Overall: acc=%.2f%% auc=%.4f",
        overall_acc * 100, overall_auc
    )
    list(
        models=models,
        predictions=preds,
        performance=perf,
        overall=list(acc=overall_acc, auc=overall_auc)
    )
}

# Helpers #####

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

# Bin a spectra object into `nbin` equal-width bins distributed across the
# non-ignored part of the ppm range. The kept domain is `[min(cs), max(cs)]`
# minus all intervals in `igrs`. Its total length is split into `nbin`
# bins of equal width; each kept sub-interval gets a proportional share.
# Returns a numeric matrix with one row per spectrum and named columns
# (rounded bin centers).
bin_spectra <- function(spectra, igrs=list(), nbin=1000) {
    cs0 <- spectra[[1]]$cs
    lo <- min(cs0); hi <- max(cs0)
    # Build kept intervals as [lo, hi] minus the union of igrs.
    if (length(igrs) == 0) {
        kept <- list(c(lo, hi))
    } else {
        ig <- t(vapply(igrs, function(r) c(min(r), max(r)), numeric(2)))
        ig <- ig[order(ig[, 1]), , drop=FALSE]
        # Merge overlapping ignore regions.
        merged <- ig[1, , drop=FALSE]
        for (i in seq_len(nrow(ig))[-1]) {
            if (ig[i, 1] <= merged[nrow(merged), 2]) {
                merged[nrow(merged), 2] <- max(
                    merged[nrow(merged), 2], ig[i, 2]
                )
            } else {
                merged <- rbind(merged, ig[i, , drop=FALSE])
            }
        }
        # Complement within [lo, hi].
        kept <- list()
        cur <- lo
        for (i in seq_len(nrow(merged))) {
            a <- max(lo, merged[i, 1]); b <- min(hi, merged[i, 2])
            if (cur < a) kept[[length(kept) + 1]] <- c(cur, a)
            cur <- max(cur, b)
        }
        if (cur < hi) kept[[length(kept) + 1]] <- c(cur, hi)
    }
    if (length(kept) == 0) stop("All ppm range is ignored.", call.=FALSE)

    # Distribute nbin bins across kept intervals proportional to length.
    lens <- vapply(kept, function(r) r[2] - r[1], numeric(1))
    total <- sum(lens)
    bw <- total / nbin
    # Per-interval bin counts, rounded but summing to nbin.
    raw <- lens / bw
    nb_each <- floor(raw)
    rem <- nbin - sum(nb_each)
    if (rem > 0) {
        ord <- order(raw - nb_each, decreasing=TRUE)
        nb_each[ord[seq_len(rem)]] <- nb_each[ord[seq_len(rem)]] + 1L
    }

    elo <- numeric(0); ehi <- numeric(0); centers <- numeric(0)
    for (i in seq_along(kept)) {
        if (nb_each[i] == 0) next
        e <- seq(kept[[i]][1], kept[[i]][2], length.out=nb_each[i] + 1)
        elo <- c(elo, e[-length(e)])
        ehi <- c(ehi, e[-1])
        centers <- c(centers, (e[-length(e)] + e[-1]) / 2)
    }

    n <- length(spectra)
    X <- matrix(0, nrow=n, ncol=length(centers))
    for (i in seq_len(n)) {
        cs <- spectra[[i]]$cs
        si <- spectra[[i]]$si
        ord <- order(cs)
        cs <- cs[ord]; si <- si[ord]
        cum <- c(0, cumsum(si))
        idx_hi <- findInterval(ehi, cs)
        idx_lo <- findInterval(elo, cs)
        X[i, ] <- cum[idx_hi + 1] - cum[idx_lo + 1]
    }
    colnames(X) <- sprintf("%.4f", centers)
    rownames(X) <- get_names(spectra)
    X
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
    s="lambda.min", nworkers=1, verbosity=1, ...
) {
    stopifnot(
        inherits(object, "mdm"), is_int(nworkers, 1),
        is_num(s, 1) || is_str(s),
        is_spectra(newdata) || is.matrix(newdata) || is.data.frame(newdata)
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
        m <- object$params
        logv(
            "Deconvoluting %d spectra with %d nworkers",
            length(newdata), nworkers
        )
        decons <- if (m$npmax > 0) {
            deconvolute(
                x=newdata, sfr=m$sfr, igrs=m$igrs %||% list(),
                verbose=verbosity >= 2,
                use_rust=m$use_rust, npmax=m$npmax, nworkers=nworkers
            )
        } else {
            deconvolute(
                x=newdata, sfr=m$sfr, igrs=m$igrs %||% list(),
                verbose=verbosity >= 2,
                use_rust=m$use_rust, nfit=m$nfit, smit=m$smit,
                smws=m$smws, delta=m$delta, npmax=0, nworkers=nworkers
            )
        }
        logv("Aligning spectra with %d nworkers", nworkers)
        als <- align_decons(
            x=decons, maxShift=m$maxShift, verbose=verbosity >= 2,
            nworkers=nworkers, ref=object$ref, full=FALSE
        )
        Xn <- si_mat(als, maxCombine=m$maxCombine, peakPos=m$peakPos)
        Xn <- Xn[, m$peakPos, drop=FALSE]
    } else {
        Xn <- as.matrix(newdata)
    }
    logv("Predicting with s=%s", as.character(s))
    requireNamespace("glmnet", quietly=TRUE)
    score <- as.numeric(predict(object$model, newx=Xn, s=s, type="link"))
    prob <- as.numeric(predict(object$model, newx=Xn, s=s, type="response"))
    pred <- predict(object$model, newx=Xn, s=s, type="class")[, 1]
    lvs <- object$model$glmnet.fit$classnames
    pred <- factor(pred, levels=lvs)
    if (type == "all") {
        return(data.frame(link=score, prob=prob, class=pred))
    }
    if (type == "class") return(pred)
    if (type == "prob") return(prob)
    if (type == "link") return(score)
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
    stats::coef(object$model, s="lambda.min", ...)
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
    stopifnot(inherits(object, "mdm"), is.list(object$params))
    pp <- c(
        "npmax", "nfit", "smit", "smws", "delta",
        "maxShift", "maxCombine"
    )
    out <- object$params[pp]
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

#' @export
#' @rdname mdm_methods
predict.bm <- function(
    object, newdata,
    type=c("all", "prob", "class", "link"),
    s="lambda.min", ...
) {
    stopifnot(inherits(object, "bm"))
    type <- match.arg(type)
    Xn <- if (is_spectra(newdata)) {
        Xb <- bin_spectra(
            newdata, igrs=object$params$igrs, nbin=object$params$nbin
        )
        Xb[, object$params$feat_names, drop=FALSE]
    } else {
        as.matrix(newdata)
    }
    requireNamespace("glmnet", quietly=TRUE)
    score <- as.numeric(
        stats::predict(object$model, newx=Xn, s=s, type="link")
    )
    prob <- as.numeric(
        stats::predict(object$model, newx=Xn, s=s, type="response")
    )
    pred <- stats::predict(object$model, newx=Xn, s=s, type="class")[, 1]
    pred <- factor(pred, levels=object$model$glmnet.fit$classnames)
    if (type == "all") {
        return(data.frame(link=score, prob=prob, class=pred))
    }
    if (type == "class") return(pred)
    if (type == "prob") return(prob)
    score
}
