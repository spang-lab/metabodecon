# Get Default Cache Directory for Deconvolution

Returns the temporary session-scoped cache directory used by
deconvolution functions such as
[`deconvolute()`](https://spang-lab.github.io/metabodecon/reference/deconvolute.md)
and the internal `grid_deconvolute_spectrum()` helper.

## Usage

``` r
decon_cachedir()
```

## Value

Path to the shared temporary cache directory for deconvolution.

## Author

2026 Tobias Schmidt: initial version.

## Examples

``` r
decon_cachedir()
#> [1] "/tmp/RtmpTeNkY0/metabodecon/cache/deconvs"
```
