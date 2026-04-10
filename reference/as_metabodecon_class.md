# Convert to a Metabodecon Object

Convert a object to a Metabodecon object.

## Usage

``` r
as_spectrum(x, sf = c(1000, 1e+06))

as_decon0(x, sf = NULL, spectrum = NULL, optional = TRUE)

as_decon1(
  x,
  sf = c(1000, 1e+06),
  spectrum = NULL,
  sfr = NULL,
  wshw = NULL,
  bwc = 2
)

as_decon2(
  x,
  sf = c(1000, 1e+06),
  spectrum = NULL,
  sfr = NULL,
  wshw = NULL,
  bwc = 2
)

as_spectra(
  x,
  file_format = "bruker",
  expno = 10,
  procno = 10,
  raw = FALSE,
  silent = TRUE,
  force = FALSE
)

as_decons0(x, sfs = list(c(1000, 1e+06)), spectra = list(NULL), nworkers = 1)

as_decons1(
  x,
  sfs = list(c(1000, 1e+06)),
  spectra = list(NULL),
  sfrs = list(NULL),
  wshws = list(NULL),
  bwc = 2,
  nworkers = 1
)

as_decons2(
  x,
  sfs = list(c(1000, 1e+06)),
  spectra = list(NULL),
  sfrs = list(NULL),
  wshws = list(NULL),
  bwc = 2,
  nworkers = 1
)
```

## Arguments

- x:

  The object to convert.

- sf:

  Scale factor used during Only required if `x` is a decon0 object.

- spectrum, spectra:

  The `spectrum`/`spectra` object corresponding to `x` as returned by
  [`read_spectrum()`](https://spang-lab.github.io/metabodecon/reference/read_spectrum.md)
  /
  [read_spectra](https://spang-lab.github.io/metabodecon/reference/read_spectra.md).
  Only required if `x` is a decon0 object.

- optional:

  Logical. If `TRUE`, the two optional elements `signal_free_region` and
  `range_water_signal_ppm` are included in the returned `decon0` object.

- sfr, sfrs:

  `sfr` should be a vector specifying the borders of the signal free
  region. `sfrs` should be a list of such vectors. Only required if `x`
  is a `decon0` object where element `signal_free_region` is missing (or
  a `decons0` objected containing such `decon0` objects).

- wshw, wshws:

  `wshw` should specify the half width of the water signal region.
  `wshws` should be a list of such values. Only required if `x` is a
  `decon0` object where element `range_water_signal_ppm` is missing (or
  a `decons0` objected containing such `decon0` objects).

- bwc:

  Level of backwards compatibility. If `bwc == 0`, bug fixes introduced
  after version 0.2.2 of Metabodecon are not used. If `bwc == 1`, new
  features introduced after version 0.2.2 of Metabodecon (e.g. faster
  algorithms) are not used. If `bwc == 2`, all bug fixes and features
  introduced after version 0.2.2 are used. Support for `bwc == 0` will
  be removed in 'metabodecon v2.0'.

- file_format:

  The file_format of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.

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

- silent:

  If `TRUE`, no output will be printed to the console.

- force:

  If `TRUE`, try to continue when encountering errors and print info
  messages instead. To hide these messages as well, set `silent = TRUE`.

- sfs:

  List of scale factors. Only required if `x` is a list of decon0
  objects.

- nworkers:

  Number of workers for parallel processing.

## Value

An object of the specified class.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
dirpath <- metabodecon_file("sim_subset")
spectra <- read_spectra(dirpath)
spectrum <- spectra[[1]]
decons1 <- generate_lorentz_curves_sim(spectra)
decon1 <- generate_lorentz_curves_sim(spectrum)
decon2 <- as_decon2(decon1)
```
