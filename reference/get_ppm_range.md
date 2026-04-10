# Get PPM Range covered by Spectra

Helper function of
[`align()`](https://spang-lab.github.io/metabodecon/reference/align.md).
Should not be called directly by the user.

Returns the ppm range across all peaks of the provided deconvoluted
spectra.

Direct usage of this function has been deprecated with metabodecon
version 1.4.3 and will be removed with metabodecon version 2.0.0.

**\[deprecated\]**

## Usage

``` r
get_ppm_range(spectrum_data, full_range = FALSE)
```

## Arguments

- spectrum_data:

  A list of deconvoluted spectra as returned by
  [`generate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/generate_lorentz_curves.md).

- full_range:

  If TRUE, the full range of the spectra is returned. If FALSE, only the
  range from the lowest to the highest peak center is returned.

## Value

A vector containing the lowest and highest ppm value over all peaks of
the provided deconvoluted spectra.

## Author

2021-2024 Wolfram Gronwald: initial version.  
2024-2025 Tobias Schmidt: refactored initial version.

## Examples

``` r
spectrum_data <- generate_lorentz_curves(
    data_path = sim[1:2],
    nfit = 3,
    sfr = c(3.55, 3.35),
    wshw = 0,
    ask = FALSE,
    verbose = FALSE
)
ppm_rng <- get_ppm_range(spectrum_data)
print(ppm_rng)
#> [1] 3.36275 3.51725
```
