# Metabodecon Models

**\[experimental\]**

**WARNING: These functions are experimental and must not be used in
production. Their API is very likely to change in
non-backwards-compatible ways over the next few weeks.**

Utilities for fitting, tuning and benchmarking 'metabodecon models'
(mdm).

A mdm is a essentially a glmnet lasso model, fitted on a feature matrix
X, obtained by deconvoluting and aligning spectra and snapping their
peaks to a shared reference spectrum. Lambda is chosen via
cross-validation on X. Deconvolution parameters (`npmax` or
`nfit`/`smit`/`smws`/`delta`), alignment parameter `maxShift`, and
peak-combining parameters (`maxCombine`, `combineMethod`) are tunable
hyperparameters.

When `npmax > 0`, deconvolution runs an internal grid search over
smoothing and fitting parameters and keeps at most `npmax` peaks. The
explicit `nfit`/`smit`/`smws`/`delta` values are ignored in that case.
When `npmax == 0`, the explicit parameters are used directly.

`fit_mdm()` deconvolutes spectra, aligns detected peaks, combines them
via
[`get_si_mat()`](https://spang-lab.github.io/metabodecon/reference/get_si_mat.md)
and fits one lasso model via
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html).
Lambda is selected internally by cross-validation but no performance
metrics are reported. Use `cv_mdm()` to tune preprocessing parameters
and `benchmark_mdm()` for unbiased performance estimates.

`cv_mdm()` evaluates a grid of preprocessing parameter combinations. For
each grid row it builds a feature matrix, runs
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
with a fixed fold assignment, and records the held-out accuracy and AUC
at `lambda.min`. Returns the model with the best AUC and the full
performance grid.

`benchmark_mdm()` wraps `cv_mdm()` in an outer cross-validation loop to
estimate end-to-end predictive performance. It returns the per-fold
models and held-out predictions.

## Usage

``` r
fit_mdm(
  spectra,
  y,
  sfr = NULL,
  use_rust = TRUE,
  nworkers = 1,
  verbosity = 2,
  seed = 1,
  cadir = decon_cachedir(),
  check = TRUE,
  npmax = 1000,
  nfit = 5,
  smit = 2,
  smws = 5,
  delta = 6.4,
  maxShift = 200,
  maxCombine = 50,
  combineMethod = 2,
  nfolds = 10
)

cv_mdm(
  spectra,
  y,
  sfr = NULL,
  use_rust = TRUE,
  nworkers = 1,
  verbosity = 2,
  seed = 1,
  cadir = decon_cachedir(),
  check = TRUE,
  combineMethod = 2,
  pgrid = get_pgrid("dynamic2"),
  ignore = NULL,
  warm_cache = TRUE,
  nfolds = 10
)

benchmark_mdm(
  spectra,
  y,
  sfr = NULL,
  use_rust = TRUE,
  verbosity = 2,
  seed = 1,
  cadir = decon_cachedir(),
  nwo = 1,
  nwi = half_cores(),
  combineMethod = 2,
  pgrid = get_pgrid(),
  nfo = 5,
  nfl = 10
)

get_pgrid(conf = "dynamic2")
```

## Arguments

- spectra:

  List-like spectra object with `cs` and `si` vectors.

- y:

  Factor vector with class labels for each spectrum.

- sfr:

  Signal free region. See
  [`deconvolute()`](https://spang-lab.github.io/metabodecon/reference/deconvolute.md)
  for details.

- use_rust:

  Logical. Whether to use the Rust backend.

- nworkers:

  Number of workers used by `fit_mdm()` and `cv_mdm()` for deconvolution
  and alignment.

- verbosity:

  Integer. Verbosity level; each nested call decrements by 1. Messages
  print when `verbosity >= 1`.

- seed:

  Random seed used for fold assignment in `benchmark_mdm()` and for
  inner cross-validation splits in `cv_mdm()` and `fit_mdm()`.

- cadir:

  Directory used to cache deconvolution results across grid points and
  processes. Defaults to
  [`decon_cachedir()`](https://spang-lab.github.io/metabodecon/reference/decon_cachedir.md).

- check:

  Logical. Whether to validate inputs at function entry.

- npmax:

  Maximum number of peaks to retain. When `npmax > 0`, deconvolution
  runs an internal grid search and the explicit
  `nfit`/`smit`/`smws`/`delta` are ignored. Set to 0 to use explicit
  deconvolution parameters instead.

- nfit:

  Number of Lorentz-curve fitting iterations (used when `npmax == 0`).

- smit:

  Number of smoothing iterations (used when `npmax == 0`).

- smws:

  Smoothing window size (used when `npmax == 0`).

- delta:

  Peak-filter threshold (used when `npmax == 0`).

- maxShift:

  Maximum alignment shift.

- maxCombine:

  Maximum peak-combining distance in datapoints, passed to
  [`get_si_mat()`](https://spang-lab.github.io/metabodecon/reference/get_si_mat.md).
  Interpretation depends on `combineMethod`.

- combineMethod:

  Peak-combining backend passed to
  [`get_si_mat()`](https://spang-lab.github.io/metabodecon/reference/get_si_mat.md).
  `1`: column merging via
  [`combine_peaks()`](https://spang-lab.github.io/metabodecon/reference/combine_peaks.md).
  `2` (default): nearest-neighbour snapping to reference peaks.

- nfolds:

  Number of folds for the `cv.glmnet()` call in `fit_mdm()` and
  `cv_mdm()`. Default 10.

- pgrid:

  Data frame of preprocessing parameter combinations as returned by
  `get_pgrid()`.

- ignore:

  Optional integer vector of sample indices to exclude.

- warm_cache:

  Logical. Whether to pre-populate the disk cache before starting the
  grid search.

- nwo:

  Number of workers for the outer cross-validation in `benchmark_mdm()`.
  Each outer worker holds a full copy of `spectra` and `y`.

- nwi:

  Number of workers used inside each outer fold for `cv_mdm()` and
  [`predict.mdm()`](https://spang-lab.github.io/metabodecon/reference/mdm_methods.md).

- nfo:

  Number of outer folds in `benchmark_mdm()`.

- nfl:

  Number of folds for the `cv.glmnet()` call inside `benchmark_mdm()`.
  Passed as `nfolds` to `fit_mdm()`.

- conf:

  Character string selecting a predefined parameter grid configuration.

## Value

`fit_mdm()` returns an object of class `mdm` with elements `model` (a
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
object), `ref` (the reference alignment spectrum) and `meta` (list of
preprocessing parameters).

`cv_mdm()` returns an `mdm` object with an additional element `pgrid`
containing the performance grid.

`benchmark_mdm()` returns a list with elements:

- `models`: List of `mdm` objects, one per outer fold.

- `predictions`: Data frame with columns `fold`, `true`, `link`, `prob`,
  `pred`.

## Details

### Disk caching across processes

`cv_mdm()` and `benchmark_mdm()` can share deconvolution results through
directory `cadir`. This on-disk cache lets parallel workers reuse
expensive preprocessing results while keeping RAM use manageable. Disk
caching can be disabled by setting `cadir = NULL`, but this will
increase runtime drastically (from a few minutes/hours to several days).

### RAM caching within processes

`fit_mdm()` caches the most recent deconvolution and alignment results
in RAM when the input `spectra` carries a `"hash"` attribute. This
avoids redundant preprocessing when `cv_mdm()` evaluates several
settings on the same spectra.

### Memory usage by parallel workers

Each outer worker in `benchmark_mdm()` keeps its own copy of `spectra`
and `y`, so memory use grows roughly with `nwo`. Large `nwo` values can
therefore require substantial RAM for large input datasets.

Inner workers (`nwi`) are used inside each outer worker for
deconvolution, alignment and prediction on held-out spectra. They
process only the spectra needed for the current task, so their memory
impact is much smaller. In RAM-constrained settings, use `nwo = 1` and
increase `nwi` instead.

The total amount of processes spawned is `nwo * nwi.`

## Examples

``` r
if (FALSE) { # \dontrun{



     # -~-~-~ Inputs -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

     aki <- read_aki_data()
     spectra <- aki$spectra
     attr(spectra, "hash") <- rlang::hash(spectra)
     y <- factor(aki$meta$type, levels = c("Control", "AKI"))
     names(y) <- rownames(aki$meta)
     sfr <- c(11, -2)
     cadir <- cachedir("deconvs", persistent = TRUE)



     # -~-~-~-~-~-~-~-~-~ Single -~-~-~-~-~-~-~-~-~

     # Best "simple" model from aki.R (0.797 [0.75 - 0.83])
     mdm <- fit_mdm(
         spectra, y, sfr=NULL,
         nfit=5, smit=3, smws=3, delta=3, npmax=0,
         maxShift=128, maxCombine=256, combineMethod=1,
         nworkers=half_cores(), cadir=cadir
     )

     mdm <- fit_mdm(
         spectra, y, sfr = NULL,
         nfit=5, smit=0, smws=0, delta=0, npmax=1000,
         maxShift=64, maxCombine=16, combineMethod=1,
         nworkers=half_cores(), cadir=cadir
     )



     # -~-~-~-~-~-~-~-~-~ Search -~-~-~-~-~-~-~-~-~

     mdm_grid_stat1 <- cv_mdm(
         spectra, y, pgrid=get_pgrid("static1"),
         nworkers=half_cores(), use_rust=TRUE, cadir=cadir
     )
     saveRDS(mdm_grid_stat1, "tmp/mdm_grid_stat1.rds")


     mdm_grid_stat2 <- cv_mdm(
         spectra, y, pgrid=get_pgrid("static2"),
         nworkers=half_cores(), use_rust=TRUE, cadir=cadir
     )
     saveRDS(mdm_grid_stat2, "tmp/mdm_grid_stat2.rds")
     #
     # Best static2: acc=0.82, auc=0.89
     # smit=2, smws=3, delta=6, nfit=7, npmax=0, maxShift=100, maxCombine=30


     mdm_grid_stat3 <- cv_mdm(
         spectra, y, pgrid=get_pgrid("static3"),
         nworkers=half_cores(), use_rust=TRUE, cadir=cadir
     )
     saveRDS(mdm_grid_stat3, "tmp/mdm_grid_stat3.rds")
     #
     # Best static3: acc=79.34% auc=0.8589
     # smit=2, smws=7, delta=8, nfit=10, npmax=0, maxShift=100, maxCombine=30


     mdm_grid_dyn2 <- cv_mdm(
         spectra, y, pgrid=get_pgrid("dynamic2"),
         nworkers=half_cores(), use_rust=TRUE, cadir=cadir
     )

     # -~-~-~ Full Benchmark -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
     bm <- benchmark_mdm(spectra, y, sfr, cadir=cadir)

     # -~-~-~ Interactive Development -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
     stub(fit_mdm, spectra=spectra, y=y, sfr=sfr)
     stub(cv_mdm, spectra=spectra, y=y, sfr=sfr)
     stub(benchmark_mdm, spectra=spectra, y=y, sfr=sfr)
} # }
```
