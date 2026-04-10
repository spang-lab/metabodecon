# Plot lorentz curves for variable range

Plots the original spectrum and all generated Lorentz curves and save
the result as png under the filepath.

Superseded by
[`plot_spectrum()`](https://spang-lab.github.io/metabodecon/reference/plot_spectrum.md)
since metabodecon v1.2.0. Will be replaced with metabodecon v2.

**\[deprecated\]**

## Usage

``` r
plot_lorentz_curves_save_as_png(
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
  writing files to disk. Default is TRUE.

## Value

NULL, called for side effects.

## See also

[`MetaboDecon1D()`](https://spang-lab.github.io/metabodecon/reference/MetaboDecon1D.md),
[`plot_triplets()`](https://spang-lab.github.io/metabodecon/reference/plot_triplets.md),
[`plot_spectrum_superposition_save_as_png()`](https://spang-lab.github.io/metabodecon/reference/plot_spectrum_superposition_save_as_png.md)

## Author

2020-2021 Martina Haeckl: initial version.  
2024-2025 Tobias Schmidt: Minor updates to pass CRAN checks

## Examples

``` r
sim <- metabodecon_file("bruker/sim_subset")
sim_decon <- generate_lorentz_curves_sim(sim)
png_dir <- tmpdir("sim_decon/pngs", create = TRUE)
plot_lorentz_curves_save_as_png(sim_decon, out_dir = png_dir, ask = FALSE)
#> Plot Lorentz curves of sim_01
#> Plot Lorentz curves of sim_02
dir(png_dir, full.names = TRUE)
#> [1] "/tmp/RtmpZMz2EC/metabodecon/sim_decon/pngs/sim_01_lorentz_curves.png"
#> [2] "/tmp/RtmpZMz2EC/metabodecon/sim_decon/pngs/sim_02_lorentz_curves.png"
```
