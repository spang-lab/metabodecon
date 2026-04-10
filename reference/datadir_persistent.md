# Persistent Data Directory

Returns the path to the persistent data directory where metabodecon's
data sets are stored. This directory equals the data directory returned
by [`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) plus
additional path normalization.

## Usage

``` r
datadir_persistent()
```

## Value

Path to the persistent data directory.

## See also

[`datadir()`](https://spang-lab.github.io/metabodecon/reference/datadir.md),
[`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
datadir_persistent()
#> [1] "/home/runner/.local/share/R/metabodecon"
```
