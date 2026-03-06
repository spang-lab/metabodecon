## Helper Functions

Based on these instructions we're gonna implement the following helper functions:

`bin_spectra(spectra, regions, binwidth)`\
For creating a matrix of binned features from a list of spectra. Binning will be
done by splitting the given `regions` into bins of width `binwidth` and summing
the intensities of all data points that fall into each bin. The resulting binned
features will be stored in a matrix, where rows correspond to spectra and
columns correspond to bins.

`get_quantile_reference(X)`\
For computing the reference distribution for quantile normalization from a
feature matrix `X`.

`quantile_normalize(X, ref)`\
For applying quantile normalization to a feature matrix `X` using a reference
distribution `ref`.

`train_model(spectra, cls, cost, gamma, nfeat, normalize)`\
For training a radial SVM classifier with margin regularization strength
`cost`, RBF kernel width `gamma` and `nfeat` number of features. Features are
obtained by binning the given input spectra using `bin_spectra()` and
selecting the top `nfeat` based on t-scores from a Welch t-test, comparing the
two classes in `cls`. If `normalize` is TRUE, quantile normalization is
applied to the binned features before feature selection and SVM training. The
feature-selection based on t-scores is described in 'Gronwald et al., 2011'.
The returned object will contain the trained SVM model, the
quantile-normalisation reference distribution and the indices of the selected
features.

`find_best_params(X, y, costs, gammas, nfeats, normalize, k)`\
For training multiple SVMs across a grid using k-fold cross-validation to
select the best parameters. In this vignette we use `costs = (1:10) / 10`,
`gammas = 2^(-12:-5)` and `nfeats = 2:3`. We limit the number of features to
such a narrow range, because the original paper states, that the average
optimum is 2.4 across validation folds with a standard deviation of 0.5. This
means, that only values of 2 and 3 can have been observed (e.g. 2,2,2,3,3) or
the SD would have been higher. Since the goal of this article is not to rerun
the whole analysis workflow, but to reproduce the reported results as closely
as possible, there is no need to search a wider range of feature counts.

## The Final Classification Pipeline

The final classification pipeline will then be implemented as follows:

```{r pipeline-cv-cleaned}
benchmark_aki_classification_cv_cleaned <- function() {
  # Read spectra and metadata
  path <- metabodecon::download_example_datasets()
  aki_path <- file.path(path, "bruker", "aki")
  spectra <- metabodecon::read_spectra(aki_path)
  meta <- read.csv(file.path(aki_path, "s_MTBLS24.txt"), sep = "\t")

  # Bin spectra and prepare y vector for classification
  X <- bin_spectra(spectra)
  typ <- meta$Factor.Value.Acute.Kidney.Injury. == "Acute Kidney Injury"
  y <- factor(typ, levels = c(0, 1), labels = c("Control", "AKI"))

  # Run classification across 5 folds and collect results
  test_ids <- get_test_ids(seq_len(nrow(meta)), nfolds = 5, y = y)
  res <- lapply(test_ids, benchmark_aki_classification_holdout, X, y)

  # Format results as 5 x 5 table
  results <- do.call(rbind, res)
  colnames(results) <- c("cost", "gamma", "nfeat", "auc", "acc")
  rownames(results) <- paste0("Fold", 1:5)

  # Add means and standard deviations across folds as additional rows
  means <- colMeans(results)
  sds <- apply(results, 2, sd)
  rbind(results, Mean = means, SD = sds)
  results
}
```

where `benchmark_aki_classification_holdout()` is defined as follows:

```{r pipeline-holdout}
benchmark_aki_classification_holdout <- function(te, X, y) {
  tr <- setdiff(seq_len(nrow(X)), te)
  costs <- (1:10) / 10
  gammas <- 2^(-12:-5)
  nfeats <- 2:3
  p <- find_best_params_cached(X[tr, ], y[tr], costs, gammas, nfeats, TRUE, 5)
  svm <- train_model(X[tr, ], y[tr], p$cost, p$gamma, p$nfeat)
  yhat <- predict_svm(svm, X[te, , drop = FALSE])
  cls <- factor(yhat$cls, levels = c(0, 1), labels = levels(y))
  auc <- auc(y[te], yhat$prob)
  acc <- mean(y[te] == cls)
  c(p$cost, p$gamma, p$nfeat, auc, acc)
}
```

### Affy quantile normalization

To mirror the affy-style quantile normalization workflow, we also create
`X_aqn` using quantile normalization across samples on the full matrix.

```{r aqn-heatmap-complex, out.width = "100%"}
#| fig.cap: |
#|   Affy-style quantile-normalized data (X aqn).
X_aqn <- t(preprocessCore::normalize.quantiles(t(as.matrix(X))))
rownames(X_aqn) <- rownames(X)
colnames(X_aqn) <- colnames(X)
attr(X_aqn, "bin_centers") <- attr(X, "bin_centers")
visualize_feature_matrix(X_aqn, y)
```

### PQN Normalization

Another common normalization method is Probabilistic Quotient Normalization
(PQN), which is based on the idea of normalizing each spectrum by a reference
spectrum that represents the typical sample in the dataset. The reference
spectrum is usually computed as the median spectrum across all samples. PQN is
particularly useful when there are large variations in the overall intensity of
the spectra, which is often the case in metabolomics data. At the end of this
section you can see a comparison of the results obtained with quantile
normalization and PQN.

```{r pqn-heatmap-complex, out.width = "100%"}
#| fig.cap: |
#|   PQN-normalized data: left sample summaries on a symlog2 x-scale,
#|   right heatmap over ppm bins.
ref_pqn <- get_pqn_reference(X)
X_pqn <- pqn_normalize(X, ref_pqn)
rownames(X_pqn) <- rownames(X)
visualize_feature_matrix(X_pqn, y)
```

### Creatinine normalization

Following the ADPKD supplementary methods, we normalize each spectrum by the
creatinine CH2 signal around 4.06 ppm.

```{r crn-heatmap-complex, out.width = "100%"}
#| fig.cap: |
#|   Creatinine-normalized data (X crn) using the bin nearest 4.06 ppm.
ppm <- as.numeric(colnames(X))
j_cr <- which.min(abs(ppm - 4.06))
X_crn <- sweep(X, 1, X[, j_cr], "/")
rownames(X_crn) <- rownames(X)
colnames(X_crn) <- colnames(X)
attr(X_crn, "bin_centers") <- attr(X, "bin_centers")
visualize_feature_matrix(X_crn, y)
```


### Hyperparameter Tuning Old

First, we train one model for fold 1.1 and inspect its predictions.

```{r tune-one}
v <- val1[[1]]
t <- setdiff(tr[[1]], v)
Xt <- X[t, , drop = FALSE]
yt <- y[t]
Xv <- X[v, , drop = FALSE]
yv <- y[v]
c0 <- 0.8
g0 <- 2^-9
p0 <- 3

m <- train_model(Xt, yt, cost = c0, gamma = g0, nfeat = p0, normalize = TRUE)
pr <- predict_svm(m, Xv); cl <- factor(pr$cls, levels = c(0, 1), labels = levels(y))
a <- auc(yv, pr$prob); ac <- mean(yv == cl)
print_svm_brief(m, X, cost = c0, gamma = g0, nfeat = p0, auc = a, acc = ac)
plot_prob_scatter_roc(pr$prob, as_binary01(yv), "Fold 1.1 predictions", "Fold 1.1 ROC")

```

Now we evaluate the full grid for fold 1.1 and then average AUCs across all
inner folds of outer fold 1.

```{r tune-grid}
#| fig.cap: |
#|   Mean inner-fold AUC surfaces for
#|   hyperparameter tuning at nfeat = 2 (left) and nfeat = 3
#|   (right), across cost and gamma.
costs <- (1:10) / 10
gammas <- 2^(-12:-5)
nfeats <- 2:3
grid_demo <- plot_inner_grid_surfaces(X, y, tr[[1]], val1, costs, gammas, nfeats)
grid_demo$best
```

In this vignette we follow the standard strategy: compute one AUC per inner
validation fold and average these AUCs across folds to select parameters.


### Performance Evaluation

Finally, we repeat the procedure for all outer folds, collect fold-wise metrics,
and pool all outer-fold predictions for one global AUC/accuracy estimate.

```{r outer-cv}
#| fig.cap: |
#|   Pooled outer-fold AKI probabilities
#|   from nested CV. Each point is one sample; color encodes class.
#|   The dashed line indicates the 0.5 classification threshold.
outer <- plot_outer_cv_predictions(
  X,
  y,
  te,
  costs,
  gammas,
  nfeats
)

print(round(outer$perf, 4))
c(global_auc = outer$auc_all, global_acc = outer$acc_all)
```

The resulting values can now be compared to the reported 83% AUC and 76%
accuracy from the original paper.

## Extending the original results

We now rerun the full nested CV benchmark for two model types (`svm`,
`lasso`) and six normalization choices (`quantile`, `pqn`, `creatinine`,
`sum`, `median`, `none`).

The benchmarking helper supports the normalization timepoint through
`norm_time`:

- `norm_time = "cv"`: normalize correctly inside each training split.
- `norm_time = "global"`: normalize once before CV, then use `none` inside CV.

The table reports mean outer-fold AUC ± SD for both settings.

```{r extending-results}
models <- c("svm", "lasso")
norms <- c("quantile", "pqn", "creatinine", "sum", "median", "none")
costs <- (1:10) / 10
gammas <- 2^(-12:-5)
nfeats <- 2:3

ext <- benchmark_extension_table(
  X = X,
  y = y,
  models = models,
  norms = norms,
  k = 5,
  costs = costs,
  gammas = gammas,
  nfeats = nfeats
)

tab <- data.frame(
  Method = paste(ext$model, ext$normalisation, sep = " + "),
  `AUC mean ± SD (CV-wise norm)` = ext$auc_cv_norm,
  `AUC mean ± SD (global pre-norm)` = ext$auc_global_norm,
  check.names = FALSE
)

print(tab, row.names = FALSE)
```
