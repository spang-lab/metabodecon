# Combine Peaks

Helper function of
[`align()`](https://spang-lab.github.io/metabodecon/reference/align.md).
Should not be called directly by the user.

Even after calling
[`speaq_align()`](https://spang-lab.github.io/metabodecon/reference/speaq_align.md),
the alignment of individual signals is not always perfect, as 'speaq'
performs a segment-wise alignment i.e. groups of signals are aligned.
For further improvements, partly filled neighboring columns are merged.
See 'Details' for an illustrative example.

Direct usage of this function has been deprecated with metabodecon
version 1.4.3 and will be removed with metabodecon version 2.0.0.

**\[deprecated\]**

## Usage

``` r
combine_peaks(
  shifted_mat,
  range = 5,
  lower_bound = 1,
  spectrum_data = NULL,
  data_path = NULL
)
```

## Arguments

- shifted_mat:

  The matrix returned by
  [`speaq_align()`](https://spang-lab.github.io/metabodecon/reference/speaq_align.md).

- range:

  Amount of adjacent columns which are permitted to be used for
  improving the alignment.

- lower_bound:

  Minimum amount of non-zero elements per column to trigger the
  alignment improvement.

- spectrum_data:

  The list of deconvoluted spectra as returned by
  [`generate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/generate_lorentz_curves.md)
  that was used to generate `shifted_mat`. No longer required since
  version 1.2 of Metabodecon.

- data_path:

  If not NULL, the returned dataframes `long` and `short` are written to
  `data_path` as "aligned_res_long.csv" and "aligned_res_short.csv".

## Value

A list containing two data frames `long` and `short`. The first data
frame contains one column for each data point in the original spectrum.
The second data frame contains only columns where at least one entry is
non-zero.

## Details

Example of what the function does:

    |            | 1    | 2    | 3    | 4    | 5    |
    |----------- |------|------|------|------|------|
    | Spectrum 1 | 0.13 | 0    | 0    | 0.11 | 0    |
    | Spectrum 2 | 0    | 0.88 | 0    | 0.12 | 0    |
    | Spectrum 3 | 0.07 | 0.56 | 0.30 | 0    | 0    |
    | Spectrum 4 | 0.08 | 0    | 0.07 | 0    | 0.07 |
    | Spectrum 5 | 0.04 | 0    | 0    | 0.04 | 0    |

becomes

    |            | 1    | 2    | 3    | 4    | 5    |
    |----------- |------|------|------|------|------|
    | Spectrum 1 | 0.13 | 0    | 0    | 0.11 | 0    |
    | Spectrum 2 | 0    | 0.88 | 0    | 0.12 | 0    |
    | Spectrum 3 | 0.07 | 0.56 | 0    | 0.30 | 0    |
    | Spectrum 4 | 0.08 | 0    | 0    | 0.07 | 0.07 |
    | Spectrum 5 | 0.04 | 0    | 0    | 0.04 | 0    |

I.e.

1.  Column 1 and 2 get NOT merged, because they have a common non-zero
    entry.

2.  Column 3 and 4 get merged, because they are in `range` of each other
    and have no common non-zero entries.

3.  Column 4 and 5 get NOT merged, because it is more beneficial to
    merge column 3 and 4, as they have more mergeable entries and after
    merging column 3 and 4, column 4 and 5 have a common non-zero entry.

## Author

2021-2024 Wolfram Gronwald: initial version.  
2024-2025 Tobias Schmidt: refactored initial version.

## Examples

``` r
deps <- c("MassSpecWavelet", "impute")
deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
if (all(deps_installed)) {
    # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
    sim_subset <- metabodecon_file("bruker/sim_subset")
    spectrum_data <- generate_lorentz_curves_sim(sim_subset)
    shifted_mat <- speaq_align(spectrum_data = spectrum_data, verbose = FALSE)
    range <- 5
    lower_bound <- 1
    obj <- combine_peaks(shifted_mat, range, lower_bound)
    str(obj)
}
#> List of 2
#>  $ short: num [1:2, 1:35] 0.000497 0.011831 0.190205 0.186338 0.003531 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:35] "3.5081" "3.50165" "3.497" "3.49475" ...
#>  $ long : num [1:2, 1:2048] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2048] "3.59" "3.58985" "3.5897" "3.58955" ...
```
