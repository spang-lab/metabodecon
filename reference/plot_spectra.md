# Plot Spectra

Plot a set of deconvoluted spectra.

## Usage

``` r
plot_spectra(
  obj,
  ...,
  foc_rgn = NULL,
  what = c("si", "sup", "supal"),
  sfy = 1e+06,
  cols = NULL,
  names = NULL,
  xlab = "Chemical Shift [ppm]",
  ylab = paste("Signal Intensity [au] /", sfy),
  mar = c(4.1, 4.1, 1.1, 0.1),
  lgd = list()
)
```

## Arguments

- obj:

  An object of type `decons0`, `decons1` or `decons2`. For details see
  [Metabodecon
  Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

- ...:

  Additional arguments passed to the conversion function.

- foc_rgn:

  Numeric vector of length 2 specifying the focus region in ppm (e.g.
  `c(3.55, 3.35)`). If NULL (default), the full spectrum is shown.

- what:

  Which signal to plot: `"supal"` (aligned superposition, default with
  fallback to `"sup"` then `"si"`), `"sup"` (superposition) or `"si"`
  (raw).

- sfy:

  Scaling factor for the y-axis.

- cols:

  Character vector of colors, one per spectrum. Defaults to
  `rainbow(n)`.

- names:

  Character vector of legend labels. Defaults to spectrum names.

- xlab:

  Label for the x-axis.

- ylab:

  Label for the y-axis.

- mar:

  A numeric vector of length 4, which specifies the margins of the plot.

- lgd:

  Logical or list. If TRUE, a legend is drawn at "topright" with
  `cex = 0.8`. If a list, its elements are passed to
  [`legend()`](https://rdrr.io/r/graphics/legend.html) to override
  position, size, etc. Set `show = FALSE` inside the list (or pass
  `lgd = FALSE`) to hide.

## Value

A plot of the deconvoluted spectra.

## See also

[`plot_spectrum()`](https://spang-lab.github.io/metabodecon/reference/plot_spectrum.md)
for a much more sophisticated plotting routine suitable for plotting a
single spectrum.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
obj <- deconvolute(sim[1:4], sfr = c(3.55, 3.35))
#> 2026-07-13 18:21:31.59 Starting deconvolution of 4 spectra using 1 worker
#> 2026-07-13 18:21:31.59 Starting deconvolution of sim_01 using R (legacy) backend
#> 2026-07-13 18:21:31.59 Removing water signal
#> 2026-07-13 18:21:31.59 Removing negative signals
#> 2026-07-13 18:21:31.59 Smoothing signals
#> 2026-07-13 18:21:31.59 Starting peak selection
#> 2026-07-13 18:21:31.60 Detected 314 peaks
#> 2026-07-13 18:21:31.60 Removing peaks with low scores
#> 2026-07-13 18:21:31.60 Removed 287 peaks
#> 2026-07-13 18:21:31.60 Initializing Lorentz curves
#> 2026-07-13 18:21:31.60 MSE at peak tiplet positions: 4.0838805770844048836921
#> 2026-07-13 18:21:31.60 Refining Lorentz Curves
#> 2026-07-13 18:21:31.60 MSE at peak tiplet positions: 0.1609359876216345797140
#> 2026-07-13 18:21:31.60 MSE at peak tiplet positions: 0.0228015051613790278862
#> 2026-07-13 18:21:31.60 MSE at peak tiplet positions: 0.0071638016610617982066
#> 2026-07-13 18:21:31.60 Formatting return object as decon2
#> 2026-07-13 18:21:31.61 Finished deconvolution of sim_01
#> 2026-07-13 18:21:31.61 Starting deconvolution of sim_02 using R (legacy) backend
#> 2026-07-13 18:21:31.61 Removing water signal
#> 2026-07-13 18:21:31.61 Removing negative signals
#> 2026-07-13 18:21:31.61 Smoothing signals
#> 2026-07-13 18:21:31.61 Starting peak selection
#> 2026-07-13 18:21:31.61 Detected 316 peaks
#> 2026-07-13 18:21:31.61 Removing peaks with low scores
#> 2026-07-13 18:21:31.61 Removed 286 peaks
#> 2026-07-13 18:21:31.61 Initializing Lorentz curves
#> 2026-07-13 18:21:31.61 MSE at peak tiplet positions: 3.8338943428876719465848
#> 2026-07-13 18:21:31.62 Refining Lorentz Curves
#> 2026-07-13 18:21:31.62 MSE at peak tiplet positions: 0.1289481941626757499630
#> 2026-07-13 18:21:31.62 MSE at peak tiplet positions: 0.0135651899090413925741
#> 2026-07-13 18:21:31.62 MSE at peak tiplet positions: 0.0025556755331531126781
#> 2026-07-13 18:21:31.62 Formatting return object as decon2
#> 2026-07-13 18:21:31.62 Finished deconvolution of sim_02
#> 2026-07-13 18:21:31.62 Starting deconvolution of sim_03 using R (legacy) backend
#> 2026-07-13 18:21:31.82 Removing water signal
#> 2026-07-13 18:21:31.82 Removing negative signals
#> 2026-07-13 18:21:31.82 Smoothing signals
#> 2026-07-13 18:21:31.83 Starting peak selection
#> 2026-07-13 18:21:31.83 Detected 333 peaks
#> 2026-07-13 18:21:31.83 Removing peaks with low scores
#> 2026-07-13 18:21:31.83 Removed 308 peaks
#> 2026-07-13 18:21:31.83 Initializing Lorentz curves
#> 2026-07-13 18:21:31.83 MSE at peak tiplet positions: 1.4917065120183621296235
#> 2026-07-13 18:21:31.83 Refining Lorentz Curves
#> 2026-07-13 18:21:31.83 MSE at peak tiplet positions: 0.0569971157280155932279
#> 2026-07-13 18:21:31.83 MSE at peak tiplet positions: 0.0065629979536274835050
#> 2026-07-13 18:21:31.83 MSE at peak tiplet positions: 0.0013913916281140697225
#> 2026-07-13 18:21:31.84 Formatting return object as decon2
#> 2026-07-13 18:21:31.84 Finished deconvolution of sim_03
#> 2026-07-13 18:21:31.84 Starting deconvolution of sim_04 using R (legacy) backend
#> 2026-07-13 18:21:31.84 Removing water signal
#> 2026-07-13 18:21:31.84 Removing negative signals
#> 2026-07-13 18:21:31.84 Smoothing signals
#> 2026-07-13 18:21:31.84 Starting peak selection
#> 2026-07-13 18:21:31.84 Detected 325 peaks
#> 2026-07-13 18:21:31.84 Removing peaks with low scores
#> 2026-07-13 18:21:31.84 Removed 299 peaks
#> 2026-07-13 18:21:31.84 Initializing Lorentz curves
#> 2026-07-13 18:21:31.84 MSE at peak tiplet positions: 2.2382155282476525748336
#> 2026-07-13 18:21:31.85 Refining Lorentz Curves
#> 2026-07-13 18:21:31.85 MSE at peak tiplet positions: 0.0843491698981613108321
#> 2026-07-13 18:21:31.85 MSE at peak tiplet positions: 0.0101688144550079323514
#> 2026-07-13 18:21:31.85 MSE at peak tiplet positions: 0.0031861616084395993388
#> 2026-07-13 18:21:31.85 Formatting return object as decon2
#> 2026-07-13 18:21:31.85 Finished deconvolution of sim_04
#> 2026-07-13 18:21:31.85 Finished deconvolution of 4 spectra in 0.266 secs
plot_spectra(obj)
```
