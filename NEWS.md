# metabodecon 2.0.0.18

* **New helper `harmonize_grid()`**: pre-aligns a corpus of spectra
  onto a single shared chemical-shift grid by integer-datapoint
  shifting. Each spectrum's `$si` is rolled by the integer offset
  from the target grid (default: median first-ppm anchor) and
  zero-padded on the vacated edge. Sub-datapoint residual is ~0.5 dp,
  two orders of magnitude smaller than typical Lorentzian widths, so
  the rounding error is invisible to downstream fits. On AKI (106
  spectra, 131k points, up to 100 dp / 0.016 ppm calibration drift),
  pre-harmonizing lifts the mean pairwise SI correlation from 0.42 to
  0.67 (all datapoints) / 0.53 to 0.65 (metabolite-only region).
  `read_aki_data()` / `cache_aki_data()` now call it automatically.

* **`cssh` is gone.** The alignment pipeline previously carried an
  auxiliary "shared grid" field `cssh` on each spectrum, populated by
  `clupa()` and read by `snap_to_ref()` / `si_mat()` / `peak_mat()` /
  `bin()`. With `harmonize_grid()` running upstream, every spectrum's
  own `$cs` IS the shared grid, so `cssh` was a redundant alias.
  Removed entirely: `ensure_cssh()`, `make_cssh()`, `pci_on_cssh()`,
  `bind_to_cssh()`. Replaced by `ensure_shared_cs()` (assertion only)
  and `pci_on_cs()`. `clupa()`, `snap_to_ref()`, `combine_peaks()`,
  `snap_nw()`, `build_consensus()`, `build_clupa_consensus()` now
  refuse to run when input spectra disagree on `$cs` — pre-call
  `harmonize_grid()` if your corpus isn't pre-aligned.

* **`supsh` is gone.** CluPA always uses the Lorentz superposition
  (`$sit$sup`) already attached at deconvolution time as the FFT
  input — this is the speaq-equivalent shape and matches v1.7.0's
  `get_sup_mat(decons2)` input to `dohCluster`. The experimental
  `supsh="triangle"` / `"rectangle"` / `"sparse"` / `"eiffel"` knob
  and the `shift_method` knob (slide vs. recompute) were removed.
  The bundled `hclust_align` remains byte-equivalent to
  `speaq::hClustAlign` (see `tests/testthat/test-speaq.R`).

# metabodecon 2.0.0.17

* `clupa()` now defaults to `supsh="lorentz"` (was `"triangle"`),
  restoring the historical speaq-CluPA input shape — the full Lorentz
  superposition evaluated on `cssh`. `align()` follows. The bundled
  `hclust_align` implementation remains byte-equivalent to
  `speaq::hClustAlign` (`use_speaq=FALSE` default; flip to `TRUE` for
  the actual speaq backend). `"triangle"` / `"rectangle"` / `"sparse"`
  / `"eiffel"` are still selectable as experimental shapes.
* The shared alignment grid `cssh` is now defined as `ref$cs` exactly,
  so the reference's own peak indices are its native integer positions
  (no convert_pos round-trip). `find_ref()` no longer constructs a
  shared `cssh` upfront — it compares candidate references in ppm
  space directly, since `find_ref_ind` only needs pairwise distances.
  `make_cssh()` has been removed; `ensure_cssh()` falls back to
  `x[[1]]$cs` for standalone callers that bypass `clupa()`.

# metabodecon 2.0.0.16

* **Breaking pick rule**: both `fit_mdm()` (grid winner) and
  `fit_lasso()` (`lambda*`) now select by **accuracy with AUC as
  tiebreaker**, replacing the AUC-only pick. This stabilises the
  reported metric pair under ranger's near-0.5 probability
  squeeze (small-sample OOB averaging shrinks predicted
  probabilities toward 0.5, which can leave a high-AUC cell at a
  low headline ACC). Existing callers see the same return shape;
  only the selected row / lambda may change.

# metabodecon 2.0.0.15

* `benchmark()` now accepts a `seed` vector → repeated k-fold CV.
  `get_test_ids()` returns a flat `length(seed) * nfolds` list with
  per-element `seed` / `fold` attributes; `benchmark()` iterates the
  flat list, stamps the originating seed onto `$performance` and
  `$predictions`, and exposes per-seed mean ± sd of acc/auc on
  `$overall$acc_seed_mean/sd` / `$overall$auc_seed_mean/sd`. Per-fold
  log lines become `[seed S, fold F/k] ...` when seeds were swept.
  Backwards compatible: scalar `seed` keeps the old shape.

# metabodecon 2.0.0.14

* `bin700()` is ~35x faster per spectrum. The per-range bin
  aggregator now uses a `cumsum` + run-boundary diff instead of
  `tapply`/`factor`; in profiling, `bin700` went from ~76% of a
  10-fold bin-baseline `benchmark()` to a negligible share.

# metabodecon 2.0.0.13

* `benchmark()` now emits a single line per outer-CV fold of the form
  `[fold i/k] acc=X.XXX auc=X.XXX | mean acc=X.XXX auc=X.XXX` (per-fold
  + running mean, 3 decimals). The inner `fit_mdm()` / `predict.mdm()`
  calls are silenced at the default `verbosity=2` — bump to
  `verbosity=3` to see the grid-search output again.

# metabodecon 2.0.0.12

* `bin700()` now picks the per-peak reconstruction position in
  `x0sn → x0al → x0` priority (falling back to the next column when
  `x0sn` is `NA`, e.g. for peaks left unmatched by `snap_to_ref` beyond
  `maxCombine`). This makes the bin baseline chainable across all
  preprocessing combinations — `decon`, `decon+align`,
  `decon+align+snap` — and lets a `snap_to_ref`-driven shift propagate
  to the bin sums.

# metabodecon 2.0.0.11

* `bin700()` is now a `feat_fun`, not a `snap_fun`. It returns a
  numeric matrix with one row per spectrum and 700 columns of bin
  sums (no peakPos filtering — the bin layout is hardcoded). Per-
  spectrum dispatch: raw spectra are binned from `$si`; deconvoluted
  / aligned spectra are reconstructed via `lorentz_sup(cs, lcpar)`
  (preferring `x0al`) and then binned. Lets snap and bin be chosen
  independently: e.g. `decon=deconvolute, align=clupa,
  snap=identity_snap, feat=bin700` runs the binning baseline on
  CluPA-aligned reconstructions; `snap=snap_to_ref, feat=peak_mat`
  remains the default mdm pipeline.

# metabodecon 2.0.0.10

* New `snap_fun` `bin700()` that discretizes each spectrum onto the
  Zacharias 2013 700-bin grid (300 bins 6.5-9.5 ppm + 400 bins
  0.5-4.5 ppm at 0.01 ppm width). Per-spectrum dispatch: raw spectra
  are binned from `$si`; deconvoluted / aligned spectra are first
  reconstructed via `lorentz_sup(cs, lcpar)` (using `x0al` when
  present) and then binned. Lets a binning baseline drop into
  `fit_mdm()` as `decon_fun=identity2, align_fun=identity_align,
  snap_fun=bin700` and lets the same `bin700` be reused after a
  CluPA-aligned deconvolution.

# metabodecon 2.0.0.9

* Tightened the `fit_fun` contract: `acc_se` / `auc_se` are now
  required in the returned list (use `NA_real_` when the backend
  produces a single point estimate, as `fit_ranger()` does).
  `fit_mdm()` no longer falls back to `NA_real_` when they are
  absent.

# metabodecon 2.0.0.8

* `fit_mdm()` / `benchmark()` replace the `mog` data-frame argument
  and the `get_mog()` helper with three scalar-or-vector pipeline
  parameters — `npmax`, `maxShift`, `maxCombine` (defaults `"auto"`,
  `"auto"`, `10`). The underlying `(nfit, smit, smws, delta)` tuple
  is always selected from each spectrum's `$deg` cache via `npmax`,
  so it never appears in the public API. When any of the three is a
  vector, `fit_mdm()` iterates over their cartesian product, reuses
  the deconvolution across rows that share `npmax` and the alignment
  across rows that share `(npmax, maxShift)`, and returns the `mdm`
  with the highest `auc`. The returned object carries the best row's
  scalar `acc` / `auc` / `acc_se` / `auc_se` directly and the full
  per-row table on `$mog`. `npmax="auto"` is resolved up front so the
  same integer is stored on the model for prediction-time replay.
* `read_aki_data()` and the new `cache_aki_data()` helper materialize
  an enriched AKI dataset whose spectra already carry their `$deg`
  grids, so `fit_mdm()` / `benchmark()` skip the slow
  `grid_deconvolute_spectra()` step. The cache filename encodes the
  backend (`R` vs `rust`) and an MD5 digest of `deg`, so caches built
  with different parameter grids or backends coexist.

# metabodecon 2.0.0.7

* `fit_mdm()` / `benchmark()` now accept `maxShift="auto"` with any
  CluPA-compatible `align_fun` (was: strict `identical(align_fun, clupa)`,
  which rejected thin wrappers like the paper's `clupa_speaq`). The auto
  sweep in `find_maxShift_dip()` uses the caller's `align_fun` so the
  picked maxShift comes from the same backend as the final alignment.
* Fixed `sprintf` crash in `fit_mdm()`'s grid-row and auto-pick log lines
  when `npmax` was a character (`"auto"` / `"intrinsic"`).

# metabodecon 2.0.0.6

* Renamed `fit_ranger500()` / `predict_ranger500()` to `fit_ranger()` /
  `predict_ranger()`; default `num.trees` bumped from 1000 to 5000 so
  OOB acc/AUC is well-converged out of the box for typical mdm sample
  sizes. Dropped the synthetic `ranger500` / `lasso` subclasses — the
  `coef()` / `plot()` dispatch in [predict.mdm] now reads the native
  `ranger` / `cv.glmnet` classes.
* `fit_lasso()` now averages per-lambda OOF performance across reps
  *before* picking the optimum (was: averaged each rep's own
  `lambda.min` performance), and exposes the chosen `lambda*` via
  `model$lambda.min` so [metabodecon::predict_lasso()] consumes it
  unchanged. All reps share a single lambda path discovered by rep 1.
  Stable lambda pick → faster convergence of the reported acc/AUC.
* Unified `fit_mdm()` around a `fit_fun(X, y) -> list(model, acc, auc)`
  contract: each learner now owns its own generalization-score estimate
  (OOB for `fit_ranger500()`, repeated `cv.glmnet` OOF for
  `fit_lasso()`). The flattened-CV scaffolding inside `fit_mdm` is gone;
  the pipeline collapses to `decon -> align -> snap -> feat -> fit`
  with all five stages pluggable. `snap_fun` is now an explicit stage
  (`snap_to_ref()` default, plus `combine_peaks()`, new
  `snap_nw_blind()` for Needleman-Wunsch with a label-blind consensus,
  and new `identity_snap()`). `fit_mdm2()`, `fit_mdm3()`, `fit_mdm4()`,
  `benchmark2()`, `bootstrap_mdm()` and `cv_mdm()` are removed; their
  use cases are covered by combinations of `fit_mdm()`'s pluggable
  stages and `benchmark()`.

# metabodecon 2.0.0.5

* `clupa()` now owns the shared `cssh` grid and per-spectrum `supsh`
  rather than `deconvolute()`. `decon2` objects no longer carry `cssh`
  or `sit$supsh`; both are attached on demand when alignment starts.
* New supsh shapes: `"triangle"` (default) and `"rectangle"` — narrow
  bounded-support shapes that are 10-80x faster than the Lorentzian
  on a 128k-point grid and that fix the FFT cross-correlator's
  tail-spillover sensitivity. `"lorentz"` and `"sparse"` remain
  available for backwards compatibility / coarse drift.
* `clupa()` rebuilds the supsh from the (shifted) peak list after
  each FFT shift instead of sliding the old vector and edge-padding
  the vacated columns. Controlled via the new `shift_method=` arg
  (default `"auto"` picks rebuild for sparse shapes, slide for
  lorentz). Speaq-backend byte-equivalence no longer holds for
  recompute mode (use `shift_method="slide"` to restore parity).
* `clupa()` / `align()` gain a `y=` argument for class-aware
  alignment: when supplied, the reference is a CluPA-aligned
  consensus built from one representative per class via
  `build_clupa_consensus()` (new private helper), so peaks present
  in only one class still have a matching reference column.

# metabodecon 2.0.0.4

* New `fit_mdm4()`: `deconvolute -> clupa -> snap_nw -> learner` pipeline,
  with joint inner-CV selection of `(maxShift, gap_tol, lambda)` and
  consensus reference built from CluPA-aligned positions. Replaces
  `fit_mdm2`'s `snap_to_ref` post-stage with the 1-to-1 NW snap. Both
  `glmnet` and `ranger` learners supported.
* `snap_nw()` and `build_consensus()` gain a `pos_field` argument so
  they can operate on `x0al` (post-CluPA aligned positions) instead of
  raw `x0`.
* `predict.mdm()` gains a `kind="clupa_nw"` branch that re-runs CluPA +
  consensus NW snap on new spectra at predict time.

# metabodecon 2.0.0.3

* New `snap_nw()` and `build_consensus()` alignment helpers.
  `snap_nw` is a Needleman-Wunsch drop-in for `snap_to_ref`; the
  pairwise DP runs in C (`src/align_dp.c`). `build_consensus` builds
  a class-aware consensus peak list from training spectra so future
  spectra can be aligned to a single fixed reference at predict time.
* New `fit_mdm3()`: `deconvolute -> consensus-NW -> glmnet` pipeline
  with internal repeated k-fold CV for joint (gap_tol, lambda)
  hyperparameter selection. Saves the entire CluPA stage relative to
  `fit_mdm2()`.

# metabodecon 2.0.0.2

* `sim2` regenerated with three discriminating peaks (factors
  `(1.25, 1.25, 1/1.25)`) instead of five, per-peak ppm jitter sd
  raised to 4 datapoints (~0.00060 ppm) and global ppm shift sd
  lowered to 8 datapoints (~0.00120 ppm). RNG seed changed from 42
  to 15 (chosen so that the supervised grid search on the
  training half places the three discriminating peaks in the top
  10 ranger features by permutation importance).
  `attr(sim2, "true_x0")` is now a 3-vector.

# metabodecon 2.0.0.1

* Default `grid_deconvolute_spectrum()` / `grid_deconvolute_spectra()` grid
  changed to `(nfit=10, smit=1:3, smws=c(3,5,7,9), delta=(1:5)*1.6)` (60 cells).

# metabodecon 1.7.0

* `deconvolute()` gained `npmax`, `igrs`, and `cachedir` parameters for limiting
  the number of peaks, ignoring ppm ranges, and caching results to disk.
* `deconvolute()` `use_rust` now accepts `0.5` to select an experimental new R
  backend that produces results identical to the Rust backend.
* `align()` gained a `method` parameter: `1` = speaq, `2` = built-in CluPA
  reimplementation (default), `3` = peak-based pairwise alignment directly on
  deconvoluted peak parameters (experimental).
* `align()` gained a `full` parameter. When `FALSE`, the aligned superposition
  is not reconstructed, saving time during grid searches.
* `get_si_mat()` gained `maxCombine`, `combineMethod`, `ref`, and `drop_zero`
  parameters for peak combining and reference-based snapping.
* `plot_spectra()` gained `foc_rgn`, `what`, `cols`, and `names` parameters for
  focus-region zooming, signal selection, custom colors, and legend labels.
* `draw_spectrum()` can now display true/false/missed peaks.
* `tree()` gained `show.counts`, `files.first`, and `max.entries` parameters.
  Added `tree_preview()` as a compact alias.
* Added `fit_mdm()` for fitting lasso models on deconvoluted NMR spectra,
  with built-in cross-validated grid search over a preprocessing grid (`mog`).
* Added `benchmark()` for nested-CV performance estimation over `fit_mdm()`.
* Added `identity2()` no-op decon function for skipping deconvolution in
  `fit_mdm()` (e.g. for binning baselines).
* Added `get_mog()` for predefined model fitting grids.
* Added S3 methods for `mdm` objects: `predict`, `print`, `coef`, `plot`,
  `summary`.
* Added `c()`, `format()`, and `summary()` methods for all public and private
  class families.
* Optimized `lorentz_sup()`: auto-selects a compiled C backend (~10x faster)
  when available, otherwise uses a vectorized R backend (~3x faster).
  Controllable via `options(metabodecon.lorentz_sup_version)`.
* Added AKI example dataset.

# metabodecon 1.6.3

* Improved exception for r-universe machine in `test-install_mdrb.R`.

# metabodecon 1.6.2

* Added an exception to `test-install_mdrb.R` so the test also works on
  r-universe machines.

# metabodecon 1.6.1

* Minor documentation improvements.

# metabodecon 1.6.0

* Updated Getting-Started vignette to use `deconvolute()` instead of
  `generate_lorentz_curves()`.
* Fixed a bug in `find_peaks()` that sometimes caused the borders of a peak to
  be chosen suboptimally. The new implementation is also about 100 times faster.

# metabodecon 1.5.2

* Removed all links to the old `TODOS.md` and `ARCHIVE.md` from the package
  documentation, as these files were removed with v1.5.1.

# metabodecon 1.5.1

* Enabled slow tests in automatic R-CMD-check workflow.
* Improved examples to also work in case the Bioconductor dependencies of
  'speaq' are not installed.
* Improved Github Actions to test and handle the case of missing Bioconductor
  dependencies.
* Disabled the automatic testing of the `check_mdrb_deps()` example, as it can
  take longer than 5 seconds on some systems.
* Removed `TODOS.md` and `ARCHIVE.md` from the package.

# metabodecon 1.5.0

* Soft-deprecated `generate_lorentz_curves()`
* Improved `install_mdrb()`. Installation is now done from R Universe instead of
  Github, allowing installation of pre-compiled binaries, which is way faster.
* Added argument `verbose` to `check_mdrb_deps()`.
* Added the possibility to provide plain `spectra` objects to `plot_spectra()`.
  Previously, only deconvoluted spectra were accepted, i.e. objects of class
  `decons0`, `decons1` or `decons2`.
* Improved `align()`. Integrals are now calculated as `A * pi`, representing the
  area under the Lorentzian curve as an improper integral, rather than being
  bounded by the signal range as a definite integral. That's faster and more
  accurate.

# metabodecon 1.4.3

* Improved questions asked by `deconvolute()` in interactive mode
* Added function `get_si_mat()` for extracting a matrix of aligned signal
  integrals from `aligns` objects
* Added author details to every function
* Added deprecation notes to `get_ppm_range()`, `gen_feat_mat()`,
  `speaq_align()`, `combine_peaks()`, `dohCluster()`,
  `calculate_lorentz_curves`, `generate_lorentz_curves()` and
  `generate_lorentz_curves_sim()`.

# metabodecon 1.4.2

- Fixed `aaa_Get_Started` entry in Manual
- Changed default value of argument `verbose` from FALSE to TRUE for function
  `align()` and `deconvolute()`.
- Added argument `install_deps` to `align()`. If non-CRAN dependencies required
  by `align()` are missing and `install_deps` is TRUE, these dependencies are
  now installed automatically. If `install_deps` is NULL (default), the user is
  asked interactively for confirmation before attempting the install.

# metabodecon 1.4.1

- Added `get_started()` and `metabodecon-package` to manual
- Improved `plot_spectrum()` default margins.
- Improved plots shown during `deconvolute()`: SFR and WSHW are now both shown
  as rectangles instead of lines.
- Improved `install_mdrb()` example
- Improved `Get_Started` article
- Included `Get_Started` article as vignette within the package

# metabodecon 1.4.0

- Improved Github Workflow (GWF) to test installation on a clean
  Windows/Linux/Mac OS with R pre-installed, but without R-tools and any
  packages.
- Added `use_rust` option to `deconvolute()`. If `use_rust` is TRUE, the
  deconvolution is done using the implementation from Rust package
  [metabodecon-rust](https://github.com/SombkeMaximilian/metabodecon-rust/tree/main).
  Using the Rust backend requires R package
  [mdrb](https://github.com/spang-lab/mdrb) (Metabodecon Rust Backend) to be
  installed first. For this purpose the following additional functions are
  provided:
    - `install_mdrb()`: Installs mdrb
    - `check_mdrb()`: Checks whether a suitable version of mdrb is already
      installed
    - `check_mdrb_deps()`: Checks whether all required system dependencies of
      mdrb are installed

# metabodecon 1.3.0

- Added Github Workflow (GWF) to test installation on a clean Windows/Linux/Mac
  OS with R pre-installed, but without R-tools and any packages. Closes Todo
  'Test Install on clean OS'.
- Fixed GWF for testing code coverage script
- Improved formatting for R-CMD-check-GWF and pkgdown-GWF
- Improved `align()`. The new implementation is faster and returns more
  information. In particular, the chemical shifts of the aligned peaks centers
  as well as the superposition of the aligned peaks are returned.
- Improved documentation and defaults for `speaq_align()` and `combine_peaks()`.
- Improved `draw_spectrum()`:
    - Added parameters `bt_text`, `lt_text`, `tp_text` and `rt_text` to
      `plot_spectrum()` to allow for full control over the text labels at the
      plot margins.
    - Added parameter `sf_vert` to `plot_spectrum()` to allow configuration of
      the height of the vertical lines drawn at the peak centers.
    - Added the option to fill the area under lorentzian curves with color.
    - Improved the legend of the plot.

# metabodecon 1.2.6

- Fixed a bug in `MetaboDecon1D()` that caused argument `file_path` to be
  interpreted as a relative path, even if it was an absolute path.

# metabodecon 1.2.5

- Fixed a bug in `read_spectrum()` that caused argument `raw` to not be passed
  on to `read_jcampdx()`.

# metabodecon 1.2.4

- Documentation updates

# metabodecon 1.2.3

- Documentation updates

# metabodecon 1.2.2

- Documentation updates

# metabodecon 1.2.1

- Documentation updates

# metabodecon 1.2.0

Finished the following tasks.

- CRAN-0: Omit "Functions for" in title
- CRAN-1: Omit "Functions for" in DESCRIPTION
- CRAN-2: Explain acronyms like NMR
- CRAN-3: Use correct reference format in DESCRIPTION
- CRAN-4: Explain return value in function docs
- CRAN-5: Remove examples from unexported functions
- CRAN-6: Fix vignettes
- CRAN-7: Check dontrun examples
- CRAN-8: Functions should not write to disk by default
- CRAN-9: Functions should not change working dir or global options
- FEATURE-01: Use temp dirs for example data
- FEATURE-02: Add minimal example dataset
- FEATURE-03: Batch Mode
- FEATURE-04: Parallelize
- FEATURE-05: Add test suite
- FEATURE-06: Return lambda in hertz
- FEATURE-07: Improve return value
- FEATURE-09 Implement `read_spectra()`
- FEATURE-11: Accept dataframes in GLC
- FEATURE-14: Provide simulated datasets
- FEATURE-15: Add lifecycle badges
- FEATURE-16: Improve multiprocessing
- FEATURE-17: Discard output
- FEATURE-18: Implement `plot_spectrum()`
- FEATURE-20: Implement `deconvolute_blood()`
- FIX-1: Prevent crashes for high smoothing
- REFACTOR-01: Combine load_spectrum functions
- REFACTOR-02: Improve Text Output (`-License`, `+Timestamps`)
- REFACTOR-04: Plotting speed
- REFACTOR-05: Speedup smoothing
- REFACTOR-06: Use a single unit as source of truth
- REFACTOR-07: Split monolithic functions into smaller parts
- REFACTOR-08: Improve docs for Metabodecon1D return value
- REFACTOR-09: Replace glc with `generate_lorentz_curves()`
- REFACTOR-10: Replace all md1d with `MetaboDecon1D()` calls
- REFACTOR-11: Implement `calc_prarp()`
- REFACTOR-12: Write compliance tests
- REFACTOR-13: Write PRARP tests

# metabodecon 1.1.1

API:

* Fixed a bug in `generate_lorentz_curves()` that caused the function to always
  use file format "bruker", even when file format "jcampdx" was specified.

Datasets:

* Fixed filenames of samples in blood dataset (renamed from `Bood_<nr>` to
  `blood_<nr>`).
* Renamed `example_datasets/jcampdx/urine/urine.dx` to
  `example_datasets/jcampdx/urine/urine_1.dx` and renamed
  `example_datasets/bruker/urine/urine/` to
  `example_datasets/bruker/urine/urine_1/`. This was done because `list.files`
  seems to return different orderings for `urine.dx` and `urine_2.dx` in
  different operating systems, whereas `urine_1.dx` and `urine_2.dx` are sorted
  the same way everywhere. This makes it easier to write clear and concise test
  cases, because we don't need to check for file ordering.

Documentation:

* Fixed broken image in
  [vignettes/FAQ.Rmd](https://github.com/spang-lab/metabodecon/blob/main/vignettes/FAQ.Rmd).

Testing:

* Added unit tests for `generate_lorentz_curves()`.
* Enabled parallel processing for unit tests.
* Created initial versions of
  `tests/testthat/test-generate_lorentz_curves-[1-4].R`.
* Added `generate_lorentz_curves_v2()` to
  `DESCRIPTION/Config/testthat/start-first`.
* Adjusted existing tests to use the updated version of `example_datasets`
  (sample `urine` was renamed to `urine_1`, as mentioned in above in section
  *Datasets*)

Internal:

* Added functions `%||%`, `msg()` and `msgf` to `R/util.R`.
* Added elements `range_water_signal_ppm` and `signal_free_region` to returned
  list of function `deconvolute_spectrum`.
* Function `with` now prints error messages to stderr even if the message stream
  is redirected.
* Copied function `deconvolution()` from `R/MetaboDecon1D.R` to `R/main_v2.R` as
  `.deconvolute_spectrum`.
* Fixed order of params in `deconvolution`.
* Fixed `download_example_datasets()`. Argument `overwrite` is passed correctly
  on to `cache_example_datasets()`.
* Changed URL of example datasets `xds$url` from
  `https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip`
  to
  `https://github.com/spang-lab/metabodecon/releases/download/v1.1.0/example_datasets.zip`.
* Improved `cache_example_datasets()`. Extraction now only is done if `extract
  == TRUE` AND the resulting folder does not yet exist (saves approx. 2-3s on
  each call). To overwrite a possible existing folder, argument `overwrite` can
  be set to TRUE.
* Fixed formatting of `test_helpers.R`
* Added linter config `.lintr`

# metabodecon 1.1.0

API:

* Improved function `download_example_datasets()` by adding caching and making
  it more stable
* Replaced function `get_data_dir()` with `datadir()` and its helper functions
  `datadir_persistent()`, `datadir_temp()` and `tempdir`
* Function `get_data_dir()` is now deprecated in favour of `datadir()`

Documentation:

* Added question about file structure to
  [vignettes/FAQ.Rmd](https://github.com/spang-lab/metabodecon/blob/main/vignettes/FAQ.Rmd)
* Created categories for function reference in
  [_pkgdown.yml](https://github.com/spang-lab/metabodecon/blob/main/_pkgdown.yml)

Datasets:

* Moved `misc/datasets` to `misc/example_datasets`
* Moved `misc/examples/usage_example.R` to `misc/code_examples/sage_example.R`

Internal:

* Added unit tests
* Removed script `check_package.R`
* Moved functions from `util.R` to `datadir.R`
* Added `grDevices`, `stats` and `utils` as internal imports
* Added lots of test helper functions in `R/test_helpers.R`
* Added function `generate_lorentz_curves_v2()` which will replace
  `generate_lorentz_curves()` as soon as we have new features AND 100% backwards
  compatibility
* Fixed bug in `with()` that caused `get_datadir_mock()` to be called after
  redirection took place causing unexpected message output
* Fixed bug in `datadir()` that caused the resulting path to end with a slash on
  Unix-like systems and without a slash on Windows, if `file` was not specified
* `RUN_SLOW_TESTS` is now set to TRUE for the CI pipeline

# metabodecon 1.0.3

API:

* Updated `get_data_dir()` to accept `"blood"` as new value for parameter
  `dataset_name`
* Updated `download_example_datasets()` to download the datasets from the github
  repo instead of the old spang-lab repo

Documentation:

* Removed table of contents from `README.md` as it's a bit overkill for approx.
  50 lines of text
* Improved documentation

Internal:

* Switched from MIT License to GPL-3 to match the license of the predecessor
  package `MetaboDecon1D`
* Added `docs` folder to `.gitignore`. Reason: we changed all vignettes to
  pkgdown articles which will be displayed only at our Github Pages website and
  can be regenerated from folder `vignettes` upon deployment.
* Created `TODOS.md` and added it to `.Rbuildignore` (Update 2025-09-14: TODOS
  are no longer tracked in `TODOS.md`, but outside of the repository. To
  retrieve the last actively maintained version of `TODOS.md`, checkout commit
  8b1f61b, i.e., v1.5.0.)
* Improved `.gitignore`

# metabodecon 1.0.2

* Minor URL and spelling adjustments to pass CRAN checks

# metabodecon 1.0.1

* Fixed some spelling errors.
* Removed unused `CONTRIBUTE.md` (instead a section within `README.md` is used)

# metabodecon 1.0.0

* Initial CRAN submission.
