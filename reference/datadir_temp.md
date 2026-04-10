# Temporary Data Directory

Returns the path to the temporary data directory where metabodecon's
data sets are stored. This directory equals subdirectory 'data' of
metabodecons temporary session directory
[`tmpdir()`](https://spang-lab.github.io/metabodecon/reference/tmpdir.md)
plus additional path normalization.

## Usage

``` r
datadir_temp()
```

## Value

Returns the path to the temporary data directory.

## See also

[`tmpdir()`](https://spang-lab.github.io/metabodecon/reference/tmpdir.md),
[`datadir()`](https://spang-lab.github.io/metabodecon/reference/datadir.md),
[`datadir_persistent()`](https://spang-lab.github.io/metabodecon/reference/datadir_persistent.md)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
datadir_temp()
#> [1] "/tmp/RtmpTeNkY0/metabodecon/data"
```
