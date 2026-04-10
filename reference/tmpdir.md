# Temporary Session Directory

Returns the path to metabodecon's temporary session directory. This
directory equals subdirectory 'metabodecon' of R's temporary session
directory [`base::tempdir()`](https://rdrr.io/r/base/tempfile.html) plus
additional path normalization.

## Usage

``` r
tmpdir(subdir = NULL, create = FALSE)
```

## Arguments

- subdir:

  Optional subdirectory within the temporary session directory.

- create:

  Whether to create the directory if it does not yet exist.

## Value

Returns the path to the temporary session directory.

## See also

[`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md)
[`datadir_persistent()`](https://spang-lab.github.io/metabodecon/reference/datadir_persistent.md)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
tmpdir()
#> [1] "/tmp/RtmpTeNkY0/metabodecon"
tmpdir("simulate_spectra")
#> [1] "/tmp/RtmpTeNkY0/metabodecon/simulate_spectra"
```
