# Generate Feature Matrix.

Helper function of
[`align()`](https://spang-lab.github.io/metabodecon/reference/align.md).
Should not be called directly by the user.

Generates a list of elements required by
[`speaq_align()`](https://spang-lab.github.io/metabodecon/reference/speaq_align.md).
See 'Value' for a detailed description of the list elements.

Direct usage of this function has been deprecated with metabodecon
version 1.4.3 and will be removed with metabodecon version 2.0.0.

**\[deprecated\]**

## Usage

``` r
gen_feat_mat(
  data_path,
  ppm_range = get_ppm_range(data_path),
  si_size_real_spectrum = length(data_path$y_values),
  scale_factor_x = 1000,
  warn = TRUE
)
```

## Arguments

- data_path:

  A list of deconvoluted spectra as returned by
  [`generate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/generate_lorentz_curves.md).
  In older versions, this could also be the path passed to
  [`generate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/generate_lorentz_curves.md),
  but this is deprecated and will trigger a warning. See 'Details' for
  more information.

- ppm_range:

  The ppm range over which your signals are distributed.

- si_size_real_spectrum:

  Number of data points in your spectra.

- scale_factor_x:

  The x scale factor used during the deconvolution.

- warn:

  Whether to print a warning in case a file path is passed to
  `data_path` instead of a list of deconvoluted spectra.

## Value

A list with the following elements:

`data_matrix`: A data.frame where each row corresponds to one spectrum
and each column to one data point, i.e. for 10 input spectra with 131072
data points each `data_matrix` would have dimensions 10 x 131072.

`peakList`: A list of vectors, where each vector contains the indices of
the peaks in the corresponding spectrum. The indices increase from left
to right, i.e. the smallest index corresponds to the highest ppm value,
as the ppm values decrease from left to right.

`w`: A list of vectors where each vector contains the "position
parameter" of the peaks in the corresponding spectrum.

`A`: A list of vectors where each vector contains the "area parameter"
of the peaks in the corresponding spectrum.

`lambda`: A list of vectors where each vector contains the "width
parameter" of the peaks in the corresponding spectrum.

## Details

Before version 1.2 of metabodecon, the deconvolution functions
`generate_lorentz_curves` and `MetaboDecon1D` wrote their output
partially as txt files to their input folder. Back then,
`gen_feat_mat()` used those txt files as input to generate the feature
matrix. Since version 1.2 these txt files are no longer created by
default, to prevent accidental modifications of the input folders.
Therefore, the recommended way to pass the required information to
`gen_feat_mat()` is to directly pass the output of
[`generate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/generate_lorentz_curves.md)
to `gen_feat_mat()`. However, to stay backwards compatible, the name of
parameter `data_path` was not changed and passing an actual path to
`data_path` is still possible, but will result in a warning (unless
`warn` is set to `FALSE`).

## Author

2021-2024 Wolfram Gronwald: initial version.  
2024-2025 Tobias Schmidt: refactored initial version.

## Examples

``` r
sim_subset <- metabodecon_file("sim_subset")
decons <- generate_lorentz_curves_sim(sim_subset)
obj <- gen_feat_mat(decons)
str(obj, 2, give.attr = FALSE)
#> List of 5
#>  $ data_matrix: num [1:2, 1:2048] 0.0122 0.0105 0.0122 0.0106 0.0123 ...
#>  $ peakList   :List of 2
#>   ..$ sim_01: num [1:27] 486 529 560 630 702 ...
#>   ..$ sim_02: num [1:30] 552 590 636 691 766 ...
#>  $ w          :List of 2
#>   ..$ sim_01: num [1:27] 1.56 1.52 1.49 1.42 1.35 ...
#>   ..$ sim_02: num [1:30] 1.5 1.46 1.41 1.36 1.28 ...
#>  $ A          :List of 2
#>   ..$ sim_01: num [1:27] -0.000158 -0.060544 -0.001124 -0.155302 -0.019516 ...
#>   ..$ sim_02: num [1:30] -0.00377 -0.05931 -0.00255 -0.14508 -0.01957 ...
#>  $ lambda     :List of 2
#>   ..$ sim_01: num [1:27] -0.00713 -0.00763 -0.008 -0.00838 -0.00732 ...
#>   ..$ sim_02: num [1:30] -0.01002 -0.00794 -0.01318 -0.00856 -0.00757 ...
```
