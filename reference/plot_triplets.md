# Plot peak triplets for variable range

Plots the peak triplets for each peak detected by
[`MetaboDecon1D()`](https://spang-lab.github.io/metabodecon/reference/MetaboDecon1D.md)
and stores the plots as png at `outdir`.

Superseded by
[`plot_spectrum()`](https://spang-lab.github.io/metabodecon/reference/plot_spectrum.md)
since metabodecon v1.2.0. Will be replaced with v2.

**\[deprecated\]**

## Usage

``` r
plot_triplets(
  deconv_result,
  x_range = c(),
  y_range = c(),
  out_dir = ".",
  ask = TRUE
)
```

## Arguments

- deconv_result:

  Saved result of the MetaboDecon1D() function

- x_range:

  Row vector with two entries consisting of the ppm start and the ppm
  end value to scale the range of the x-axis (optional)

- y_range:

  Row vector with two entries consisting of the ppm start and the ppm
  end value to scale the range of the y-axis (optional)

- out_dir:

  Directory to save the png files (optional)

- ask:

  Logical value to ask the user if the png files should be saved in the
  specified directory (optional)

## Value

No return value, called for side effect of plotting.

## See also

[`MetaboDecon1D()`](https://spang-lab.github.io/metabodecon/reference/MetaboDecon1D.md),
[`calculate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/calculate_lorentz_curves.md),
[`plot_lorentz_curves_save_as_png()`](https://spang-lab.github.io/metabodecon/reference/plot_lorentz_curves_save_as_png.md),
[`plot_spectrum_superposition_save_as_png()`](https://spang-lab.github.io/metabodecon/reference/plot_spectrum_superposition_save_as_png.md)

## Author

2020-2021 Martina Haeckl: initial version.  
2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks.

## Examples

``` r
sim <- metabodecon_file("bruker/sim_subset")
sim_decon <- generate_lorentz_curves_sim(sim)
png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
plot_triplets(sim_decon, out_dir = png_dir, ask = FALSE)
#> Plot triplets of sim_01
#> Plot triplets of sim_02
dir(png_dir, full.names = TRUE)
#> [1] "/tmp/RtmpoVTN85/metabodecon/sim_decon/pngs/sim_01_lorentz_curves.png"    
#> [2] "/tmp/RtmpoVTN85/metabodecon/sim_decon/pngs/sim_01_peak_triplets.png"     
#> [3] "/tmp/RtmpoVTN85/metabodecon/sim_decon/pngs/sim_01_sum_lorentz_curves.png"
#> [4] "/tmp/RtmpoVTN85/metabodecon/sim_decon/pngs/sim_02_lorentz_curves.png"    
#> [5] "/tmp/RtmpoVTN85/metabodecon/sim_decon/pngs/sim_02_peak_triplets.png"     
#> [6] "/tmp/RtmpoVTN85/metabodecon/sim_decon/pngs/sim_02_sum_lorentz_curves.png"
```
