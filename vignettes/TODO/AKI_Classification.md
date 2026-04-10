---
title: "AKI Classification"
output: rmarkdown::html_vignette
css: styles.css
vignette: >
  %\VignetteIndexEntry{AKI Classification}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r set-defaults, echo=FALSE, results=FALSE, message=FALSE}
knitr::opts_chunk$set(
    fig.dim = c(5, 5),
    fig.show = "hold",
    out.height = "auto",
    eval = FALSE,
    echo = FALSE,
    results = FALSE,
    message = FALSE
)
knitr::opts_chunk$set(
  fig.dim = c(7, 5),
  out.width = "49%",
  eval = FALSE
)
```

This article benchmarks three preprocessing strategies for supervised AKI
classification on 1D urine NMR spectra:

1. metabodecon deconvolution + alignment,
2. speaq-style wavelet alignment,
3. equidistant binning.

The AKI data comprise spectra measured 24 h after surgery in Erlangen.
Phenodata and file paths are expected in an AKI directory containing `pheno.csv`
and Bruker spectrum subfolders.

All expensive analysis steps are cached to:

`R_user_dir("metabodecon", "cache")/<functionName>/<argsHash>/result.rds`.

# Setup

```{r setup, echo=TRUE, message=FALSE, warning=FALSE}
library(metabodecon)

aki_dir_candidates <- c(
  normalizePath(file.path("..", "misc", "example_datasets", "bruker", "aki"), mustWork = FALSE),
  normalizePath(file.path("misc", "example_datasets", "bruker", "aki"), mustWork = FALSE),
  tryCatch(metabodecon::metabodecon_file("example_datasets/bruker/aki"), error = function(e) "")
)
aki_dir_candidates <- unique(aki_dir_candidates[nzchar(aki_dir_candidates)])
aki_dir <- aki_dir_candidates[file.exists(file.path(aki_dir_candidates, "pheno.csv"))][1]
if (is.na(aki_dir) || !nzchar(aki_dir)) {
  stop("Could not locate AKI example dataset folder containing pheno.csv")
}
pheno <- read.csv(file.path(aki_dir, "pheno.csv"), stringsAsFactors = FALSE)

paths <- file.path(aki_dir, pheno$File)
labels <- as.integer(pheno$HasAKI)
ok <- file.exists(paths)
paths <- paths[ok]
labels <- labels[ok]
pheno <- pheno[ok, , drop = FALSE]

cache_root <- tools::R_user_dir("metabodecon", "cache")
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

cache_key <- function(args) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(args, algo = "xxhash64"))
  }
  txt <- paste(capture.output(str(args)), collapse = "|")
  as.character(abs(sum(utf8ToInt(txt))))
}

cached_eval <- function(functionName, args, expr, overwrite = FALSE) {
  key <- cache_key(args)
  cache_dir <- file.path(cache_root, functionName, key)
  cache_file <- file.path(cache_dir, "result.rds")
  if (!overwrite && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  value <- eval(substitute(expr), envir = parent.frame())
  saveRDS(value, cache_file)
  value
}

split <- metabodecon:::aki_train_test_split(labels, ratio = c(2, 1), seed = 1)
train_paths <- paths[split$train_idx]
test_paths <- paths[split$test_idx]
train_labels <- labels[split$train_idx]
test_labels <- labels[split$test_idx]

pick_rep_indices <- function(mse) {
  mse <- as.numeric(mse)
  stopifnot(length(mse) > 4)
  q <- stats::quantile(mse, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  nearest <- function(target) which.min(abs(mse - target))
  idx <- c(
    which.max(mse),
    nearest(q[1]),
    nearest(q[2]),
    nearest(q[3]),
    which.min(mse)
  )
  idx <- unique(idx)
  if (length(idx) < 5) {
    ord <- order(mse, decreasing = TRUE)
    idx <- unique(c(idx, ord))[1:5]
  }
  idx
}

plot_rep_grid <- function(specs, mse, idx, title_prefix = "") {
  stopifnot(length(idx) == 5)
  labs <- c("worst", "q25", "median", "q75", "best")
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(3, 2), mar = c(1.8, 2.2, 2.2, 0.8))
  for (i in seq_along(idx)) {
    j <- idx[i]
    metabodecon::plot_spectrum(specs[[j]], sub2 = FALSE, sub3 = FALSE, frame = TRUE)
    graphics::title(main = sprintf("%s (%s)\nMSE = %.4g", labs[i], basename(names(specs)[j]), mse[j]), cex.main = 0.8)
  }
  graphics::plot(
    seq_along(mse), mse,
    pch = 16, cex = 0.6,
    xlab = "Spectrum index", ylab = "MSE",
    main = sprintf("%sMSE overview", ifelse(nchar(title_prefix) > 0, paste0(title_prefix, ": "), ""))
  )
  graphics::abline(h = stats::quantile(mse, probs = c(0.25, 0.5, 0.75)), lty = 2)
  graphics::points(idx, mse[idx], col = "red", pch = 19)
}

plot_alignment_heatmap <- function(mat, main = "Aligned spectra heatmap") {
  z <- log1p(mat)
  z <- t(scale(t(z)))
  z[!is.finite(z)] <- 0
  graphics::image(
    x = seq_len(ncol(z)),
    y = seq_len(nrow(z)),
    z = t(z[, ncol(z):1, drop = FALSE]),
    col = grDevices::hcl.colors(64, "YlOrRd", rev = TRUE),
    xlab = "Aligned feature index",
    ylab = "Spectrum index",
    main = main
  )
}

as_plot_spectrum <- function(template, si_new) {
  x <- template
  si_new <- as.numeric(si_new)
  if (length(si_new) == length(x$cs)) {
    cs_new <- x$cs
  } else {
    cs_new <- seq(max(x$cs), min(x$cs), length.out = length(si_new))
  }
  x$cs <- cs_new
  x$si <- si_new
  x$meta$fq <- stats::approx(template$cs, template$meta$fq, xout = cs_new, rule = 2)$y
  class(x) <- "spectrum"
  x
}

method_pipeline <- function(paths,
              method = c("metabodecon", "speaq", "binning"),
              install_deps = FALSE,
              verbose = FALSE) {
  method <- match.arg(method)
  raw_specs <- metabodecon:::aki_read_spectra(paths)
  sfr <- metabodecon:::aki_default_sfr(raw_specs[[1]], NULL)
  best <- metabodecon:::find_best_params(raw_specs[[1]], sfr = sfr, verbose = verbose)

  decons <- metabodecon::deconvolute(
    raw_specs,
    sfr = sfr,
    smopts = best$smopts,
    delta = best$delta,
    nfit = best$nfit,
    ask = FALSE,
    verbose = verbose,
    use_rust = TRUE
  )
  decons <- metabodecon::as_decons2(decons)

  decon_mse <- vapply(decons, function(x) mean((x$si - x$sit$sup)^2), numeric(1))

  if (method == "metabodecon") {
    aligns <- metabodecon::align(decons, maxShift = 50, maxCombine = 5, install_deps = install_deps, verbose = verbose)
    align_mse <- vapply(aligns, function(x) mean((x$si - x$sit$supal)^2), numeric(1))
    align_specs <- aligns
    align_mat <- t(metabodecon::get_si_mat(aligns))
  } else if (method == "speaq") {
    decons1 <- metabodecon::as_decons1(decons)
    feat <- metabodecon::gen_feat_mat(decons1)
    shifted <- metabodecon::speaq_align(
      feat,
      maxShift = 50,
      spectrum_data = decons1,
      show = FALSE,
      verbose = verbose
    )
    comb <- metabodecon::combine_peaks(shifted, range = 5)
    align_mat <- comb$long
    align_specs <- structure(
      lapply(seq_len(nrow(align_mat)), function(i) as_plot_spectrum(raw_specs[[i]], align_mat[i, ])),
      class = "spectra"
    )
    align_mse <- vapply(seq_along(raw_specs), function(i) {
      y <- stats::approx(align_specs[[i]]$cs, align_specs[[i]]$si, xout = raw_specs[[i]]$cs, rule = 2)$y
      mean((raw_specs[[i]]$si - y)^2)
    }, numeric(1))
  } else {
    bins <- metabodecon:::aki_bin_matrix(raw_specs, bin_width = 0.01, bin_range = NULL, bin_edges = NULL)
    align_mat <- bins$mat
    align_specs <- structure(
      lapply(seq_len(nrow(align_mat)), function(i) as_plot_spectrum(raw_specs[[i]], align_mat[i, ])),
      class = "spectra"
    )
    align_mse <- vapply(seq_along(raw_specs), function(i) {
      y <- stats::approx(align_specs[[i]]$cs, align_specs[[i]]$si, xout = raw_specs[[i]]$cs, rule = 2)$y
      mean((raw_specs[[i]]$si - y)^2)
    }, numeric(1))
  }

  list(
    raw_specs = raw_specs,
    decon_specs = decons,
    align_specs = align_specs,
    decon_mse = decon_mse,
    align_mse = align_mse,
    align_mat = align_mat
  )
}

train_metabodecon <- cached_eval(
  functionName = "method_pipeline",
  args = list(
    method = "metabodecon",
    train_paths = train_paths,
    install_deps = TRUE,
    verbose = FALSE
  ),
  expr = method_pipeline(train_paths, "metabodecon", install_deps = TRUE)
)
train_speaq <- cached_eval(
  functionName = "method_pipeline",
  args = list(
    method = "speaq",
    train_paths = train_paths,
    install_deps = TRUE,
    verbose = FALSE
  ),
  expr = method_pipeline(train_paths, "speaq", install_deps = TRUE)
)
train_binning <- cached_eval(
  functionName = "method_pipeline",
  args = list(
    method = "binning",
    train_paths = train_paths,
    install_deps = FALSE,
    verbose = FALSE
  ),
  expr = method_pipeline(train_paths, "binning", install_deps = FALSE)
)
```

# 1. Deconvolution

In this chapter, we compare deconvolution behavior across metabodecon, speaq,
and binning. For each method we show five representative spectra selected from
the MSE distribution: worst, 0.25-quantile, median, 0.75-quantile, and best.

## Metabodecon

```{r decon-metabodecon-before, echo=TRUE}
idx <- pick_rep_indices(train_metabodecon$decon_mse)
plot_rep_grid(
  train_metabodecon$raw_specs,
  mse = train_metabodecon$decon_mse,
  idx = idx,
  title_prefix = "Metabodecon before deconvolution"
)
```

```{r decon-metabodecon-after, echo=TRUE}
idx <- pick_rep_indices(train_metabodecon$decon_mse)
plot_rep_grid(
  train_metabodecon$decon_specs,
  mse = train_metabodecon$decon_mse,
  idx = idx,
  title_prefix = "Metabodecon after deconvolution"
)
summary(train_metabodecon$decon_mse)
```

## Speaq

Speaq itself performs alignment and feature extraction; for a like-for-like
deconvolution view we use the same deconvolution step before speaq alignment.

```{r decon-speaq-before-after, echo=TRUE}
idx <- pick_rep_indices(train_speaq$decon_mse)
plot_rep_grid(train_speaq$raw_specs, train_speaq$decon_mse, idx, title_prefix = "Speaq pipeline before deconvolution")
plot_rep_grid(train_speaq$decon_specs, train_speaq$decon_mse, idx, title_prefix = "Speaq pipeline after deconvolution")
summary(train_speaq$decon_mse)
```

## Binning

Binning has no explicit Lorentzian deconvolution. We therefore interpret the
binned representation as a coarse approximation step and assess it with the same
MSE-based representative spectrum strategy.

```{r decon-binning-before-after, echo=TRUE}
idx <- pick_rep_indices(train_binning$align_mse)
plot_rep_grid(train_binning$raw_specs, train_binning$align_mse, idx,
  title_prefix = "Binning before approximation")
plot_rep_grid(train_binning$align_specs, train_binning$align_mse, idx,
  title_prefix = "Binning after approximation")
summary(train_binning$align_mse)
```

The MSE summaries quantify reconstruction quality for each strategy and provide
a direct way to identify outlier spectra requiring manual QC.

# 2. Alignment

This chapter evaluates alignment output for each method using:

1. 3x2 spectrum panels before alignment,
2. 3x2 spectrum panels after alignment,
3. a heatmap of aligned intensities.

## Metabodecon

```{r align-metabodecon, echo=TRUE}
idx <- pick_rep_indices(train_metabodecon$align_mse)
plot_rep_grid(train_metabodecon$decon_specs, train_metabodecon$align_mse, idx,
  title_prefix = "Metabodecon before alignment")
plot_rep_grid(train_metabodecon$align_specs, train_metabodecon$align_mse, idx,
  title_prefix = "Metabodecon after alignment")
plot_alignment_heatmap(train_metabodecon$align_mat, "Metabodecon aligned heatmap")
summary(train_metabodecon$align_mse)
```

## Speaq

```{r align-speaq, echo=TRUE}
idx <- pick_rep_indices(train_speaq$align_mse)
plot_rep_grid(train_speaq$decon_specs, train_speaq$align_mse, idx,
  title_prefix = "Speaq before alignment")
plot_rep_grid(train_speaq$align_specs, train_speaq$align_mse, idx,
  title_prefix = "Speaq after alignment")
plot_alignment_heatmap(train_speaq$align_mat, "Speaq aligned heatmap")
summary(train_speaq$align_mse)
```

## Binning

```{r align-binning, echo=TRUE}
idx <- pick_rep_indices(train_binning$align_mse)
plot_rep_grid(train_binning$raw_specs, train_binning$align_mse, idx,
  title_prefix = "Binning before alignment")
plot_rep_grid(train_binning$align_specs, train_binning$align_mse, idx,
  title_prefix = "Binning after alignment")
plot_alignment_heatmap(train_binning$align_mat, "Binning heatmap")
summary(train_binning$align_mse)
```

Heatmaps provide an overview of global peak consistency after alignment,
complementing the five-spectrum diagnostic panels.

# 3. Classification

We now perform AKI-vs-control classification for each method using elastic-net
logistic regression (`glmnet`) with a 2:1 train/test split.

```{r classification-run, echo=TRUE, message=FALSE}
cls <- cached_eval(
  functionName = "classify_aki_patients",
  args = list(
    data = paths,
    labels = labels,
    methods = c("metabodecon", "binning", "speaq"),
    split_ratio = c(2, 1),
    seed = 1,
    install_deps = TRUE,
    verbose = FALSE
  ),
  expr = metabodecon:::classify_aki_patients(
    data = paths,
    labels = labels,
    methods = c("metabodecon", "binning", "speaq"),
    split_ratio = c(2, 1),
    seed = 1,
    install_deps = TRUE,
    verbose = FALSE
  )
)

cls$metrics
```

```{r classification-visualize, echo=TRUE}
metrics <- cls$metrics
old <- par(no.readonly = TRUE)
on.exit(par(old), add = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

barplot(
  metrics$bal_acc,
  names.arg = metrics$method,
  ylim = c(0, 1),
  ylab = "Balanced accuracy",
  main = "Classification performance"
)

barplot(
  metrics$auc,
  names.arg = metrics$method,
  ylim = c(0, 1),
  ylab = "AUROC",
  main = "ROC area"
)
```

These plots summarize predictive performance per preprocessing method under the
same train/test split and model family.

# 4. Comparison

Finally, we compare the complete pipeline behavior across methods:

1. reconstruction/alignment error distributions,
2. classification metrics,
3. qualitative visual diagnostics from representative spectra and heatmaps.

```{r comparison-summary, echo=TRUE}
comparison <- data.frame(
  method = c("metabodecon", "speaq", "binning"),
  decon_mse_median = c(
    stats::median(train_metabodecon$decon_mse),
    stats::median(train_speaq$decon_mse),
    NA_real_
  ),
  align_mse_median = c(
    stats::median(train_metabodecon$align_mse),
    stats::median(train_speaq$align_mse),
    stats::median(train_binning$align_mse)
  )
)

if (exists("cls")) {
  comparison <- merge(comparison, cls$metrics,
    by.x = "method", by.y = "method", all.x = TRUE)
}

comparison[order(comparison$auc, decreasing = TRUE), ]
```

Overall, this chaptered workflow separates deconvolution quality,
alignment quality, and downstream classification performance while keeping all
three methods directly comparable.
