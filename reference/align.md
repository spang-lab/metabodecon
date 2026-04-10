# Align Spectra

Align signals across a list of deconvoluted spectra using the 'CluPA'
algorithm from the 'speaq' package, described in Beirnaert et al. (2018)
<doi:10.1371/journal.pcbi.1006018> and Vu et al. (2011)
<doi:10.1186/1471-2105-12-405> plus the additional peak combination
described in
[`combine_peaks()`](https://spang-lab.github.io/metabodecon/reference/combine_peaks.md).

## Usage

``` r
align(
  x,
  maxShift = 50,
  maxCombine = 0,
  verbose = TRUE,
  method = 2,
  install_deps = NULL,
  nworkers = 1,
  ref = NULL,
  full = TRUE
)
```

## Arguments

- x:

  An object of type `decons1` or `decons2` as described in [Metabodecon
  Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).
  To align `decons0` objects (as returned by the now deprecated
  [MetaboDecon1D](https://spang-lab.github.io/metabodecon/reference/MetaboDecon1D.md)),
  you can use
  [`as_decons2()`](https://spang-lab.github.io/metabodecon/reference/as_metabodecon_class.md)
  to convert it to a `decons2` object first.

- maxShift:

  Maximum number of points along the "ppm-axis" a value can be moved by
  the 'speaq' package. 50 is a suitable starting value for plasma
  spectra with a digital resolution of 128K. Note that this parameter
  has to be individually optimized depending on the type of analyzed
  spectra and the digital resolution. For urine which is more prone to
  chemical shift variations this value most probably has to be
  increased. Passed as argument `maxShift` to
  [`speaq_align()`](https://spang-lab.github.io/metabodecon/reference/speaq_align.md).

- maxCombine:

  Amount of adjacent columns which may be combined for improving the
  alignment during the CluPA step. We recommend setting this to 0 and
  instead relying on the peak snapping implemented in
  [`get_si_mat()`](https://spang-lab.github.io/metabodecon/reference/get_si_mat.md).
  Since version 1.7.0 the default is 0. Non-zero values are deprecated
  and support will be removed in a future version.

- verbose:

  Whether to print additional information during the alignment process.

- method:

  Alignment backend. `1`: use the original implementation from the
  'speaq' package. `2` (default): use metabodecon's built-in
  reimplementation of the CluPA algorithm. `3`: use metabodecon's
  peak-based pairwise alignment, which works directly on the
  deconvoluted peak parameters (`x0`, `lambda`, `A`). **Method 3 is
  experimental and must not be used in production. It is very likely to
  change in non-backwards-compatible ways over the next few weeks.**

- install_deps:

  Only used when `method = 1`. 'speaq' relies on the 'MassSpecWavelet'
  and 'impute' packages. Both, 'MassSpecWavelet' and 'impute' are not
  available on CRAN, but can be installed from
  [Bioconductor](https://www.bioconductor.org/) or
  [R-Universe](https://r-universe.dev/). If `install_deps=TRUE`, these
  packages will be automatically installed from R-Universe without
  asking for confirmation. If `install_deps=NULL` (default), the user
  will be asked for confirmation before installing missing dependencies.
  If asking for confirmation is not possible or `install_deps=FALSE`,
  the function will raise an error if the packages are not installed.

- nworkers:

  Number of parallel workers for the alignment. Default is 1 (no
  parallelism).

- ref:

  Optional reference spectrum of type `align` or `decon2`. When
  supplied, all spectra in `x` are aligned towards this reference. The
  reference is prepended to `x` internally and removed from the result.
  If `NULL` (default), the reference is chosen automatically via
  [`speaq::findRef()`](https://rdrr.io/pkg/speaq/man/findRef.html).

- full:

  If `TRUE` (default), store the full aligned Lorentz-curve
  superposition in each returned spectrum. If `FALSE`, skip that
  reconstruction step to save time and memory. Skipping reconstruction
  of the superposition means the plotting routines will not work as
  expected, so only use `full=FALSE` if you are sure you don't need the
  superposition.

## Value

An object of type `align` as described in [Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
if (interactive()) {
    # Example requires an interactive R session, because in case of missing
    # dependencies the user will be asked for confirmation to install them.
    decons <- deconvolute(sim[1:2], sfr = c(3.55, 3.35))
    aligned <- align(decons)
}
```
