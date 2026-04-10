# S3 Methods for Printing Metabodecon Objects

S3 Methods for printing metabodecon objects as described in the
[Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/).

## Usage

``` r
# S3 method for class 'spectrum'
print(x, name = FALSE, ...)

# S3 method for class 'decon1'
print(x, name = FALSE, ...)

# S3 method for class 'decon2'
print(x, name = FALSE, ...)

# S3 method for class 'align'
print(x, name = FALSE, ...)

# S3 method for class 'spectra'
print(x, ...)

# S3 method for class 'decons1'
print(x, ...)

# S3 method for class 'decons2'
print(x, ...)

# S3 method for class 'aligns'
print(x, ...)
```

## Arguments

- x:

  The object to print.

- name:

  Logical. If TRUE, the name of the object is printed before the object.

- ...:

  Not used. Only accepted to comply with generic
  [`base::print()`](https://rdrr.io/r/base/print.html).

## Value

NULL, called for side effect of printing to the standard output device.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
print(sim[[1]])
#> spectrum object (2048 dp, 3.6 to 3.3 ppm)
print(sim[[1]], name = TRUE)
#> sim_01: spectrum object (2048 dp, 3.6 to 3.3 ppm)
print(sim)
#> spectra object consisting of 16 spectrum objects:
#> sim_01 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_02 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_03 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_04 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_05 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_06 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_07 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_08 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_09 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_10 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_11 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_12 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_13 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_14 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_15 (2048 datapoints from 3.28 - 3.59 ppm)
#> sim_16 (2048 datapoints from 3.28 - 3.59 ppm)
decon <- deconvolute(sim[[1]], sfr = c(3.55, 3.35))
#> 2026-04-10 07:15:56.84 Starting deconvolution of 1 spectrum using 1 worker
#> 2026-04-10 07:15:56.84 Starting deconvolution of sim_01 using R (legacy) backend
#> 2026-04-10 07:15:56.84 Removing water signal
#> 2026-04-10 07:15:56.84 Removing negative signals
#> 2026-04-10 07:15:56.84 Smoothing signals
#> 2026-04-10 07:15:56.84 Starting peak selection
#> 2026-04-10 07:15:56.84 Detected 314 peaks
#> 2026-04-10 07:15:56.84 Removing peaks with low scores
#> 2026-04-10 07:15:56.84 Removed 287 peaks
#> 2026-04-10 07:15:56.84 Initializing Lorentz curves
#> 2026-04-10 07:15:56.84 MSE at peak tiplet positions: 4.0838805770844048836921
#> 2026-04-10 07:15:56.85 Refining Lorentz Curves
#> 2026-04-10 07:15:56.85 MSE at peak tiplet positions: 0.1609359876216345797140
#> 2026-04-10 07:15:56.85 MSE at peak tiplet positions: 0.0228015051613790313556
#> 2026-04-10 07:15:56.85 MSE at peak tiplet positions: 0.0071638016610617799920
#> 2026-04-10 07:15:56.85 Formatting return object as decon2
#> 2026-04-10 07:15:56.85 Finished deconvolution of sim_01
#> 2026-04-10 07:15:56.85 Finished deconvolution of 1 spectrum in 0.017 secs
print(decon)
#> decon2 object (2048 dp, 3.6 to 3.3 ppm, 27 peaks)
```
