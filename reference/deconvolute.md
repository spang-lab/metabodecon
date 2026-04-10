# Deconvolute one or more NMR spectra

Deconvolutes NMR spectra by modeling each detected signal within a
spectrum as Lorentz Curve.

## Usage

``` r
deconvolute(
  x,
  nfit = 3,
  smopts = c(2, 5),
  delta = 6.4,
  sfr = NULL,
  wshw = 0,
  ask = FALSE,
  force = FALSE,
  verbose = TRUE,
  nworkers = 1,
  use_rust = FALSE,
  npmax = 0,
  igrs = list(),
  cadir = decon_cachedir()
)
```

## Arguments

- x:

  A `spectrum` or `spectra` object as described in [Metabodecon
  Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

- nfit:

  Integer. Number of iterations for approximating the parameters for the
  Lorentz curves. See 'Details'.

- smopts:

  Numeric vector with two entries: the number of smoothing iterations
  and the number of data points to use for smoothing (must be odd). See
  'Details'.

- delta:

  Threshold for peak filtering. Higher values result in more peaks being
  filtered out. A peak is filtered if its score is below \\\mu + \sigma
  \cdot \delta\\, where \\\mu\\ is the average peak score in the
  signal-free region (SFR), and \\\sigma\\ is the standard deviation of
  peak scores in the SFR. See 'Details'.

- sfr:

  Numeric vector with two entries: the ppm positions for the left and
  right border of the signal-free region of the spectrum. See 'Details'.

- wshw:

  Half-width of the water artifact in ppm. See 'Details'.

- ask:

  Logical. Whether to ask for user input during the deconvolution
  process. If FALSE, the provided default values will be used.

- force:

  If FALSE, the function stops with an error message if no peaks are
  found in the signal free region (SFR), as these peaks are required as
  a reference for peak filtering. If TRUE, the function instead proceeds
  without peak filtering, potentially increasing runtime and memory
  usage significantly.

- verbose:

  Logical. Whether to print log messages during the deconvolution
  process.

- nworkers:

  Number of workers to use for parallel processing. If `"auto"`, the
  number of workers will be determined automatically. If a number
  greater than 1, it will be limited to the number of spectra.

- use_rust:

  Controls the deconvolution backend. Accepts `FALSE` / `0` (legacy R
  implementation, default), `TRUE` / `1` (Rust backend via
  [mdrb](https://github.com/spang-lab/mdrb)), `0.5` (experimental new R
  implementation that produces results identical to the Rust backend;
  planned to become the new default in a future release), or `NULL`
  (auto-detect: uses Rust if available, otherwise legacy R). When set to
  `TRUE` / `1` and mdrb is not installed, an error is thrown.

- npmax:

  Integer. Maximum number of peaks allowed in the result. If
  `npmax >= 1`, the `nfit`, `smopts` and `delta` arguments are ignored
  and a grid search over predefined parameter combinations is performed
  instead. The combination with the smallest residual area ratio and
  fewer than `npmax` peaks is selected. Grid search results are cached
  to disk by default; see `cadir`.

- igrs:

  Ignore regions. List of length-2 numeric vectors specifying the start
  and endpoints of the chemical shift regions to ignore during
  deconvolution. Peaks whose centers fall inside any ignore region are
  excluded from fitting.

- cadir:

  Directory for caching grid search results and deconvolution results
  for `npmax >= 1`. Defaults to
  [`decon_cachedir()`](https://spang-lab.github.io/metabodecon/reference/decon_cachedir.md).
  If `NULL`, caching is disabled. Pass a custom path to use a different
  cache directory.

## Value

A 'decon2' object as described in [Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

## Details

First, an automated curvature based signal selection is performed. Each
signal is represented by 3 data points to allow the determination of
initial Lorentz curves. These Lorentz curves are then iteratively
adjusted to optimally approximate the measured spectrum.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
## Deconvolute a single spectrum
spectrum <- sim[[1]]
decon <- deconvolute(spectrum)
#> 2026-04-10 07:15:46.97 Starting deconvolution of 1 spectrum using 1 worker
#> 2026-04-10 07:15:46.97 Starting deconvolution of sim_01 using R (legacy) backend
#> 2026-04-10 07:15:46.97 Removing water signal
#> 2026-04-10 07:15:46.97 Removing negative signals
#> 2026-04-10 07:15:46.97 Smoothing signals
#> 2026-04-10 07:15:46.97 Starting peak selection
#> 2026-04-10 07:15:46.97 Detected 314 peaks
#> 2026-04-10 07:15:46.97 Removing peaks with low scores
#> 2026-04-10 07:15:46.97 Removed 287 peaks
#> 2026-04-10 07:15:46.98 Initializing Lorentz curves
#> 2026-04-10 07:15:46.98 MSE at peak tiplet positions: 4.0838805770844048836921
#> 2026-04-10 07:15:46.98 Refining Lorentz Curves
#> 2026-04-10 07:15:46.98 MSE at peak tiplet positions: 0.1609359876216345797140
#> 2026-04-10 07:15:46.98 MSE at peak tiplet positions: 0.0228015051613790313556
#> 2026-04-10 07:15:46.98 MSE at peak tiplet positions: 0.0071638016610617799920
#> 2026-04-10 07:15:46.98 Formatting return object as decon2
#> 2026-04-10 07:15:46.98 Finished deconvolution of sim_01
#> 2026-04-10 07:15:46.98 Finished deconvolution of 1 spectrum in 0.017 secs

## Read multiple spectra from disk and deconvolute at once
spectra_dir <- metabodecon_file("sim_subset")
spectra <- read_spectra(spectra_dir)
decons <- deconvolute(spectra, sfr = c(3.55, 3.35))
#> 2026-04-10 07:15:46.99 Starting deconvolution of 2 spectra using 1 worker
#> 2026-04-10 07:15:46.99 Starting deconvolution of sim_01 using R (legacy) backend
#> 2026-04-10 07:15:46.99 Removing water signal
#> 2026-04-10 07:15:46.99 Removing negative signals
#> 2026-04-10 07:15:46.99 Smoothing signals
#> 2026-04-10 07:15:47.00 Starting peak selection
#> 2026-04-10 07:15:47.00 Detected 314 peaks
#> 2026-04-10 07:15:47.00 Removing peaks with low scores
#> 2026-04-10 07:15:47.00 Removed 287 peaks
#> 2026-04-10 07:15:47.00 Initializing Lorentz curves
#> 2026-04-10 07:15:47.00 MSE at peak tiplet positions: 4.0838805770844048836921
#> 2026-04-10 07:15:47.00 Refining Lorentz Curves
#> 2026-04-10 07:15:47.00 MSE at peak tiplet positions: 0.1609359876216345797140
#> 2026-04-10 07:15:47.01 MSE at peak tiplet positions: 0.0228015051613790313556
#> 2026-04-10 07:15:47.01 MSE at peak tiplet positions: 0.0071638016610617799920
#> 2026-04-10 07:15:47.01 Formatting return object as decon2
#> 2026-04-10 07:15:47.01 Finished deconvolution of sim_01
#> 2026-04-10 07:15:47.01 Starting deconvolution of sim_02 using R (legacy) backend
#> 2026-04-10 07:15:47.01 Removing water signal
#> 2026-04-10 07:15:47.01 Removing negative signals
#> 2026-04-10 07:15:47.01 Smoothing signals
#> 2026-04-10 07:15:47.02 Starting peak selection
#> 2026-04-10 07:15:47.02 Detected 316 peaks
#> 2026-04-10 07:15:47.02 Removing peaks with low scores
#> 2026-04-10 07:15:47.02 Removed 286 peaks
#> 2026-04-10 07:15:47.02 Initializing Lorentz curves
#> 2026-04-10 07:15:47.02 MSE at peak tiplet positions: 3.8338943428876719465848
#> 2026-04-10 07:15:47.02 Refining Lorentz Curves
#> 2026-04-10 07:15:47.02 MSE at peak tiplet positions: 0.1289481941626757499630
#> 2026-04-10 07:15:47.02 MSE at peak tiplet positions: 0.0135651899090413786964
#> 2026-04-10 07:15:47.03 MSE at peak tiplet positions: 0.0025556755331531087749
#> 2026-04-10 07:15:47.03 Formatting return object as decon2
#> 2026-04-10 07:15:47.03 Finished deconvolution of sim_02
#> 2026-04-10 07:15:47.03 Finished deconvolution of 2 spectra in 0.035 secs
```
