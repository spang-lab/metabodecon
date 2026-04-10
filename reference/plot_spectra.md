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
#> 2026-04-10 05:59:10.80 Starting deconvolution of 4 spectra using 1 worker
#> 2026-04-10 05:59:10.80 Starting deconvolution of sim_01 using R (legacy) backend
#> 2026-04-10 05:59:10.80 Removing water signal
#> 2026-04-10 05:59:10.80 Removing negative signals
#> 2026-04-10 05:59:10.80 Smoothing signals
#> 2026-04-10 05:59:10.81 Starting peak selection
#> 2026-04-10 05:59:10.81 Detected 314 peaks
#> 2026-04-10 05:59:10.81 Removing peaks with low scores
#> 2026-04-10 05:59:10.81 Removed 287 peaks
#> 2026-04-10 05:59:10.81 Initializing Lorentz curves
#> 2026-04-10 05:59:10.81 MSE at peak tiplet positions: 4.0838805770844048836921
#> 2026-04-10 05:59:10.81 Refining Lorentz Curves
#> 2026-04-10 05:59:10.81 MSE at peak tiplet positions: 0.1609359876216345797140
#> 2026-04-10 05:59:10.81 MSE at peak tiplet positions: 0.0228015051613790313556
#> 2026-04-10 05:59:10.81 MSE at peak tiplet positions: 0.0071638016610617799920
#> 2026-04-10 05:59:10.82 Formatting return object as decon2
#> 2026-04-10 05:59:10.82 Finished deconvolution of sim_01
#> 2026-04-10 05:59:10.82 Starting deconvolution of sim_02 using R (legacy) backend
#> 2026-04-10 05:59:10.82 Removing water signal
#> 2026-04-10 05:59:10.82 Removing negative signals
#> 2026-04-10 05:59:10.82 Smoothing signals
#> 2026-04-10 05:59:10.82 Starting peak selection
#> 2026-04-10 05:59:10.82 Detected 316 peaks
#> 2026-04-10 05:59:10.82 Removing peaks with low scores
#> 2026-04-10 05:59:10.82 Removed 286 peaks
#> 2026-04-10 05:59:10.83 Initializing Lorentz curves
#> 2026-04-10 05:59:10.83 MSE at peak tiplet positions: 3.8338943428876719465848
#> 2026-04-10 05:59:10.83 Refining Lorentz Curves
#> 2026-04-10 05:59:10.83 MSE at peak tiplet positions: 0.1289481941626757499630
#> 2026-04-10 05:59:10.83 MSE at peak tiplet positions: 0.0135651899090413786964
#> 2026-04-10 05:59:10.85 MSE at peak tiplet positions: 0.0025556755331531087749
#> 2026-04-10 05:59:10.86 Formatting return object as decon2
#> 2026-04-10 05:59:10.86 Finished deconvolution of sim_02
#> 2026-04-10 05:59:10.86 Starting deconvolution of sim_03 using R (legacy) backend
#> 2026-04-10 05:59:10.86 Removing water signal
#> 2026-04-10 05:59:10.86 Removing negative signals
#> 2026-04-10 05:59:10.86 Smoothing signals
#> 2026-04-10 05:59:10.86 Starting peak selection
#> 2026-04-10 05:59:10.86 Detected 333 peaks
#> 2026-04-10 05:59:10.86 Removing peaks with low scores
#> 2026-04-10 05:59:10.86 Removed 308 peaks
#> 2026-04-10 05:59:10.86 Initializing Lorentz curves
#> 2026-04-10 05:59:10.87 MSE at peak tiplet positions: 1.4917065120183623516681
#> 2026-04-10 05:59:10.87 Refining Lorentz Curves
#> 2026-04-10 05:59:10.87 MSE at peak tiplet positions: 0.0569971157280155932279
#> 2026-04-10 05:59:10.87 MSE at peak tiplet positions: 0.0065629979536274852397
#> 2026-04-10 05:59:10.87 MSE at peak tiplet positions: 0.0013913916281140725414
#> 2026-04-10 05:59:10.87 Formatting return object as decon2
#> 2026-04-10 05:59:10.87 Finished deconvolution of sim_03
#> 2026-04-10 05:59:10.87 Starting deconvolution of sim_04 using R (legacy) backend
#> 2026-04-10 05:59:10.87 Removing water signal
#> 2026-04-10 05:59:10.87 Removing negative signals
#> 2026-04-10 05:59:10.87 Smoothing signals
#> 2026-04-10 05:59:10.88 Starting peak selection
#> 2026-04-10 05:59:10.88 Detected 325 peaks
#> 2026-04-10 05:59:10.88 Removing peaks with low scores
#> 2026-04-10 05:59:10.88 Removed 299 peaks
#> 2026-04-10 05:59:10.88 Initializing Lorentz curves
#> 2026-04-10 05:59:10.88 MSE at peak tiplet positions: 2.2382155282476525748336
#> 2026-04-10 05:59:10.88 Refining Lorentz Curves
#> 2026-04-10 05:59:10.88 MSE at peak tiplet positions: 0.0843491698981613247099
#> 2026-04-10 05:59:10.89 MSE at peak tiplet positions: 0.0101688144550079444944
#> 2026-04-10 05:59:10.89 MSE at peak tiplet positions: 0.0031861616084395915326
#> 2026-04-10 05:59:10.89 Formatting return object as decon2
#> 2026-04-10 05:59:10.89 Finished deconvolution of sim_04
#> 2026-04-10 05:59:10.89 Finished deconvolution of 4 spectra in 0.093 secs
plot_spectra(obj)
```
