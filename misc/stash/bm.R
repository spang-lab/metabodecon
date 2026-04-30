fit_bm <- function(
    spectra, y,
    regions = rbind(c(9.5, 6.5), c(4.5, 0.5)),
    binwidth = 0.01,
    seed = 1,
    nfolds = 10,
    check = TRUE
) {
    if (check) check_bm_args(spectra, y, regions, binwidth, nfolds)
    X <- bin_spectra(spectra, regions = regions, binwidth = binwidth)
    foldid <- get_foldid(y = y, nfolds = nfolds, seed = seed)
    model <- glmnet::cv.glmnet(
        X, y, family = "binomial", alpha = 1,
        foldid = foldid, keep = TRUE
    )
    structure(list(model = model, meta = list(
        regions = regions,
        binwidth = binwidth,
        feat_names = colnames(X)
    )), class = "bm")
}

predict.bm <- function(object,
                       newdata,
                       type = c("all", "prob", "class", "link"),
                       s = "lambda.min",
                       ...) {
    stopifnot(inherits(object, "bm"))
    type <- match.arg(type)
    Xn <- if (is_spectra(newdata)) {
        X <- bin_spectra(
            newdata,
            regions = object$meta$regions,
            binwidth = object$meta$binwidth
        )
        X[, object$meta$feat_names, drop = FALSE]
    } else {
        as.matrix(newdata)
    }
    score <- as.numeric(stats::predict(object$model, newx = Xn, s = s,
                                       type = "link"))
    prob <- as.numeric(stats::predict(object$model, newx = Xn, s = s,
                                      type = "response"))
    pred <- stats::predict(object$model, newx = Xn, s = s, type = "class")[, 1]
    pred <- factor(pred, levels = object$model$glmnet.fit$classnames)
    if (type == "all") return(data.frame(link = score, prob = prob, class = pred))
    if (type == "class") return(pred)
    if (type == "prob") return(prob)
    score
}

benchmark_bm <- function(
    spectra, y,
    regions = rbind(c(9.5, 6.5), c(4.5, 0.5)),
    binwidth = 0.01,
    verbosity = 2,
    seed = 1,
    nwo = 1,
    nfo = 5,
    nfl = 10
) {
    check_bm_args(spectra, y, regions, binwidth, nfo)
    check_bm_args(spectra, y, regions, binwidth, nfl)
    te_list <- get_test_ids(nfolds = nfo, nsamples = length(spectra), seed = seed, y = y)
    fids <- seq_along(te_list)
    nwoa <- min(nwo, nfo)
    models <- mcmapply(
        nwoa, fit_bm_fold,
        ignore = te_list,
        MoreArgs = list(
            spectra = spectra, y = y,
            regions = regions, binwidth = binwidth,
            seed = seed, nfolds = nfl
        )
    )
    te_spectra <- lapply(te_list, function(te) spectra[te])
    pred_metrics <- mcmapply(
        nwoa, predict.bm,
        object = models, newdata = te_spectra,
        MoreArgs = list(type = "all"),
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
    preds <- mapply(function(fid, te, pred) {
        data.frame(fold = fid, true = y[te], link = pred$link,
                   prob = pred$prob, pred = pred$class)
    }, fids, te_list, pred_metrics, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    logv("Nested CV accuracy: %.2f%%", mean(do.call(rbind, preds)$true == do.call(rbind, preds)$pred) * 100)
    list(models = models, predictions = do.call(rbind, preds))
}

fit_bm_fold <- function(ignore, spectra, y, regions, binwidth, seed, nfolds) {
    fit_bm(
        spectra = spectra[-ignore], y = y[-ignore],
        regions = regions, binwidth = binwidth,
        seed = seed, nfolds = nfolds, check = FALSE
    )
}

check_bm_args <- function(spectra, y, regions, binwidth, nfolds) {
    stopifnot(
        is_spectra(spectra), is.factor(y), length(y) == length(spectra),
        is.matrix(regions), is.numeric(regions), ncol(regions) == 2,
        is_num(binwidth, 1), binwidth > 0,
        is_int(nfolds, 1), nfolds >= 2, nfolds <= length(y)
    )
    if (!is.null(names(y)) && !identical(get_names(spectra), names(y))) {
        stop("Names of `spectra` and `y` must match and be in the same order.", call. = FALSE)
    }
    if (nlevels(y) != 2 || any(table(y) == 0)) {
        stop("`y` must contain exactly 2 non-empty classes.", call. = FALSE)
    }
    if (any(regions[, 1] <= regions[, 2])) {
        stop("Each row of `regions` must satisfy upper > lower.", call. = FALSE)
    }
    invisible(NULL)
}