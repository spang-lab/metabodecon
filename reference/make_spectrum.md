# Create a Spectrum Object

Creates a spectrum object from the provided signal intensities,
frequencies and chemical shifts.

## Usage

``` r
make_spectrum(
  si,
  cs_max,
  cs_width,
  fq_ref,
  fq_width = NULL,
  force = FALSE,
  silent = FALSE,
  name = NULL,
  path = NULL,
  type = NULL,
  simpar = NULL,
  mfs = NULL
)
```

## Arguments

- si:

  Numeric vector of signal intensities, ordered from highest to lowest
  corresponding chemical shift.

- cs_max:

  The highest chemical shift value in ppm, usually shown as left end of
  the spectrum.

- cs_width:

  The width of the spectrum in ppm.

- fq_ref:

  The reference frequency in Hz.

- fq_width:

  The width of the spectrum in Hz. Only used to check whether the values
  calculated from `cs_max`, `cs_width` and `fq_ref` match the provided
  value. If `NULL`, this check will be skipped.

- force:

  If `TRUE`, the function will not raise an error in case of
  discrepancies between the calculated and the provided spectrum width
  in Hz, but will print a info message instead. To hide this message as
  well, set `silent = TRUE`.

- silent:

  If `TRUE`, no output will be printed to the console.

- name:

  The name of the spectrum, e.g. "Blood 1" or "Urine Mouse X23D".

- path:

  The path to the spectrum file, e.g.
  "/example_datasets/bruker/urine/urine_1".

- type:

  The type of experiment, e.g. "H1 CPMG" or "H1 NOESY".

- simpar:

  The simulation parameters used to generate the spectrum.

- mfs:

  The magnetic field strength in Tesla.

## Value

A `spectrum` object as described in [Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
si <- c(1, 1, 3, 7, 8, 3, 8, 5, 2, 1)
cs_max <- 14.8
cs_width <- 20.0
fq_ref <- 600.25 * 1e6
fq_width <- 12005
spectrum <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width)
spectrum2 <- make_spectrum(si, cs_max, cs_width, fq_ref, fq_width = 12010, force = FALSE)
#> 2026-04-10 05:50:38.79 Calculated spectrum width in Hz (12005) does not match the provided value (12010). Continuing anyways, because `force` equals `TRUE`. Please note that all downstream calculations using frequencies might be wrong, so be sure to double check the results.
```
