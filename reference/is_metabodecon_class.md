# Is an Object from a Metabodecon Class?

Check if an object is an instance of a specific 'Metabodecon Class'. See
[Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html)
for a list of classes.

## Usage

``` r
is_spectrum(x, check_class = TRUE, check_contents = FALSE)

is_decon0(x)

is_decon1(x)

is_decon2(x)

is_align(x)

is_spectra(
  x,
  check_class = TRUE,
  check_contents = FALSE,
  check_child_classes = FALSE
)

is_decons0(x)

is_decons1(x)

is_decons2(x)

is_aligns(x)
```

## Arguments

- x:

  The object to check.

- check_class:

  Logical indicating whether to check the class of the object.

- check_contents:

  Logical indicating whether to check the contents of the object.

- check_child_classes:

  Logical indicating whether to check the class of each element of the
  object.

## Value

TRUE if the object is an instance of the specified class, otherwise
FALSE.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
ss <- sim[1:2]
s1 <- sim[[1]]
is_spectra(ss) # TRUE
#> [1] TRUE
is_spectrum(s1) # TRUE
#> [1] TRUE
is_spectrum(s1, check_contents = TRUE) # TRUE
#> [1] TRUE

dd <- deconvolute(ss, sfr = c(3.55, 3.35))
#> 2026-07-22 11:53:16.75 Starting deconvolution of 2 spectra using 1 worker
#> 2026-07-22 11:53:16.75 Starting deconvolution of sim_01 using R (legacy) backend
#> 2026-07-22 11:53:16.75 Removing water signal
#> 2026-07-22 11:53:16.75 Removing negative signals
#> 2026-07-22 11:53:16.75 Smoothing signals
#> 2026-07-22 11:53:16.76 Starting peak selection
#> 2026-07-22 11:53:16.76 Detected 314 peaks
#> 2026-07-22 11:53:16.76 Removing peaks with low scores
#> 2026-07-22 11:53:16.76 Removed 287 peaks
#> 2026-07-22 11:53:16.76 Initializing Lorentz curves
#> 2026-07-22 11:53:16.76 MSE at peak tiplet positions: 4.0838805770844048836921
#> 2026-07-22 11:53:16.76 Refining Lorentz Curves
#> 2026-07-22 11:53:16.76 MSE at peak tiplet positions: 0.1609359876216345797140
#> 2026-07-22 11:53:16.76 MSE at peak tiplet positions: 0.0228015051613790278862
#> 2026-07-22 11:53:16.76 MSE at peak tiplet positions: 0.0071638016610617982066
#> 2026-07-22 11:53:16.76 Formatting return object as decon2
#> 2026-07-22 11:53:16.77 Finished deconvolution of sim_01
#> 2026-07-22 11:53:16.77 Starting deconvolution of sim_02 using R (legacy) backend
#> 2026-07-22 11:53:16.77 Removing water signal
#> 2026-07-22 11:53:16.77 Removing negative signals
#> 2026-07-22 11:53:16.77 Smoothing signals
#> 2026-07-22 11:53:16.77 Starting peak selection
#> 2026-07-22 11:53:16.77 Detected 316 peaks
#> 2026-07-22 11:53:16.77 Removing peaks with low scores
#> 2026-07-22 11:53:16.77 Removed 286 peaks
#> 2026-07-22 11:53:16.77 Initializing Lorentz curves
#> 2026-07-22 11:53:16.77 MSE at peak tiplet positions: 3.8338943428876719465848
#> 2026-07-22 11:53:16.78 Refining Lorentz Curves
#> 2026-07-22 11:53:16.78 MSE at peak tiplet positions: 0.1289481941626757499630
#> 2026-07-22 11:53:16.78 MSE at peak tiplet positions: 0.0135651899090413925741
#> 2026-07-22 11:53:16.78 MSE at peak tiplet positions: 0.0025556755331531126781
#> 2026-07-22 11:53:16.78 Formatting return object as decon2
#> 2026-07-22 11:53:16.78 Finished deconvolution of sim_02
#> 2026-07-22 11:53:16.78 Finished deconvolution of 2 spectra in 0.036 secs
d1 <- dd[[1]]
is_decons0(dd) # FALSE
#> [1] FALSE
is_decons1(dd) # FALSE
#> [1] FALSE
is_decons2(dd) # TRUE
#> [1] TRUE
is_decon0(d1) # FALSE
#> [1] FALSE
is_decon1(d1) # FALSE
#> [1] FALSE
is_decon2(d1) # TRUE
#> [1] TRUE

if (interactive()) {
    # Example requires an interactive R session, because in case of missing
    # dependencies the user will be asked for confirmation to install them.
    aa <- align(dd)
    a1 <- aa[[1]]
    is_align(a1) # TRUE
    is_aligns(aa) # TRUE
}
```
