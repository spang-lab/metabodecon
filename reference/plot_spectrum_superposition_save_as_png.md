# Plot spectrum approx for variable range

Plots the original spectrum and the superposition of all generated
Lorentz curves and saves the result as png under the specified filepath.

Superseded by
[`plot_spectrum()`](https://spang-lab.github.io/metabodecon/reference/plot_spectrum.md)
since metabodecon v1.2.0. Will be replaced with v2.

**\[deprecated\]**

## Usage

``` r
plot_spectrum_superposition_save_as_png(
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

  Path to the directory where the png files should be saved. Default is
  the current working directory.

- ask:

  Logical value. Whether to ask for confirmation from the user before
  writing files to disk.

## Value

NULL, called for side effects.

## See also

[`MetaboDecon1D()`](https://spang-lab.github.io/metabodecon/reference/MetaboDecon1D.md),
[`calculate_lorentz_curves()`](https://spang-lab.github.io/metabodecon/reference/calculate_lorentz_curves.md),
[`plot_triplets()`](https://spang-lab.github.io/metabodecon/reference/plot_triplets.md),
[`plot_lorentz_curves_save_as_png()`](https://spang-lab.github.io/metabodecon/reference/plot_lorentz_curves_save_as_png.md)

## Author

2020-2021 Martina Haeckl: initial version.

## Examples

``` r
sim <- metabodecon_file("bruker/sim_subset")
sim_decon <- generate_lorentz_curves_sim(sim)
png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
plot_spectrum_superposition_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
#> Plot superposition of sim_01
#> Plot superposition of sim_02
dir(png_dir, full.names = TRUE)
#> [1] "/tmp/RtmpTeNkY0/metabodecon/sim_decon/pngs/sim_01_lorentz_curves.png"    
#> [2] "/tmp/RtmpTeNkY0/metabodecon/sim_decon/pngs/sim_01_sum_lorentz_curves.png"
#> [3] "/tmp/RtmpTeNkY0/metabodecon/sim_decon/pngs/sim_02_lorentz_curves.png"    
#> [4] "/tmp/RtmpTeNkY0/metabodecon/sim_decon/pngs/sim_02_sum_lorentz_curves.png"
```
