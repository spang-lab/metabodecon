# Deconvolute one or more NMR spectra

Deconvolutes NMR spectra by modeling each detected signal within a
spectrum as Lorentz Curve.

This function has been deprecated with metabodecon version v1.4.3 and
will be removed with version 2.0.0. Please use
[`deconvolute()`](https://spang-lab.github.io/metabodecon/reference/deconvolute.md)
instead.

**\[deprecated\]**

## Usage

``` r
generate_lorentz_curves(
  data_path,
  file_format = "bruker",
  make_rds = FALSE,
  expno = 10,
  procno = 10,
  raw = TRUE,
  nfit = 10,
  smopts = c(2, 5),
  delta = 6.4,
  sfr = c(11.44494, -1.8828),
  wshw = 0.1527692,
  ask = TRUE,
  force = FALSE,
  verbose = TRUE,
  nworkers = 1
)

generate_lorentz_curves_sim(
  data_path,
  file_format = "bruker",
  make_rds = FALSE,
  expno = 10,
  procno = 10,
  raw = TRUE,
  nfit = 10,
  smopts = c(2, 5),
  delta = 6.4,
  sfr = c(3.55, 3.35),
  wshw = 0,
  ask = FALSE,
  force = FALSE,
  verbose = FALSE,
  nworkers = 1
)
```

## Arguments

- data_path:

  The path of the file/folder containing the spectrum data. E.g.
  `"example_datasets/jcampdx/urine/urine_1.dx"` or
  `"example_datasets/bruker/urine/urine"`.

- file_format:

  The file_format of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.

- make_rds:

  Logical or character. If TRUE, stores results as an RDS file on disk.
  If a character string, saves the RDS file with the specified name.
  Should be set to TRUE if many spectra are evaluated to decrease
  computation time.

- expno, procno:

  The experiment/processing number for the file. E.g. `"10"`. Only
  relevant if `file_format` equals `"bruker"`. For details see section
  [File
  Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure)
  in the metabodecon FAQ.

- raw:

  If `FALSE`, scales the returned signal intensities based on
  information available in the spectrum metadata, in particular
  `NC_proc`. For details see `processing-reference.pdf`, available at
  <https://www.bruker.com/en.html> at section 'Services & Support \>
  Documentation & Manuals \> Magnetic Resonance \> Acquisition &
  Processing \> TopSpin Processing Commands and Parameters' (requires
  login).

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

  If `TRUE`, try to continue when encountering errors and print info
  messages instead. To hide these messages as well, set `silent = TRUE`.

- verbose:

  Logical. Whether to print log messages during the deconvolution
  process.

- nworkers:

  Number of workers to use for parallel processing. If `"auto"`, the
  number of workers will be determined automatically. If a number
  greater than 1, it will be limited to the number of spectra.

## Value

A 'decon1' object if a single spectrum was provided. A 'decons1' object
if multiple spectra were provided. See [Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html)
for details.

## Details

First, an automated curvature based signal selection is performed. Each
signal is represented by 3 data points to allow the determination of
initial Lorentz curves. These Lorentz curves are then iteratively
adjusted to optimally approximate the measured spectrum.

`generate_lorentz_curves_sim()` is identical to
`generate_lorentz_curves()` except for the defaults, which are optimized
for deconvoluting the 'sim' dataset, shipped with 'metabodecon'. The
'sim' dataset is a simulated dataset, which is much smaller than a real
NMR spectra and lacks a water signal. This makes it ideal for use in
examples. However, the default values for `sfr`, `wshw`, and `delta` in
the "normal" `generate_lorentz_curves()` function are not optimal for
this dataset. To avoid having to define the optimal parameters
repeatedly in examples, this function is provided to deconvolute the
"Sim" dataset with suitable parameters.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
## Define the paths to the example datasets we want to deconvolute:
## `sim_dir`: directory containing 16 simulated spectra
## `sim_01`: path to the first spectrum in the `sim` directory
## `sim_01_spec`: the first spectrum in the `sim` directory as a dataframe

sim_dir        <- metabodecon_file("sim_subset")
sim_1_dir      <- file.path(sim_dir, "sim_01")
sim_2_dir      <- file.path(sim_dir, "sim_02")
sim_1_spectrum <- read_spectrum(sim_1_dir)
sim_2_spectrum <- read_spectrum(sim_2_dir)
sim_spectra    <- read_spectra(sim_dir)


## Show that `generate_lorentz_curves()` and `generate_lorentz_curves_sim()`
## produce the same results:

sim_1_decon0 <- generate_lorentz_curves(
    data_path = sim_1_dir, # Path to directory containing spectra
    sfr = c(3.55, 3.35),   # Borders of signal free region (SFR) in ppm
    wshw = 0,              # Half width of water signal (WS) in ppm
    ask = FALSE,           # Don't ask for user input
    verbose = FALSE        # Suppress status messages
)
sim_1_decon1 <- generate_lorentz_curves_sim(sim_1_dir)
stopifnot(all.equal(sim_1_decon0, sim_1_decon1))


## Show that passing a spectrum produces the same results as passing the
## the corresponding directory:

decon_from_spectrum_dir <- generate_lorentz_curves_sim(sim_1_dir)
decon_from_spectrum_obj <- generate_lorentz_curves_sim(sim_1_spectrum)
decons_from_spectra_obj <- generate_lorentz_curves_sim(sim_spectra)
decons_from_spectra_dir <- generate_lorentz_curves_sim(sim_dir)

most.equal <- function(x1, x2) {
    ignore <- which(names(x1) %in% c("number_of_files", "filename"))
    equal <- all.equal(x1[-ignore], x2[-ignore])
    invisible(stopifnot(isTRUE(equal)))
}

all.equal(  decon_from_spectrum_dir, decon_from_spectrum_obj     )
#> [1] TRUE
all.equal(  decons_from_spectra_dir, decons_from_spectra_obj     )
#> [1] TRUE
most.equal( decon_from_spectrum_dir, decons_from_spectra_obj[[1]])
most.equal( decon_from_spectrum_dir, decons_from_spectra_dir[[1]])
```
