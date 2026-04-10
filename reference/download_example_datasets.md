# Download metabodecon Example Datasets

Downloads example datasets that can be used to test the functionality of
the metabodecon package. These datasets are not included in the package
by default due to size constraints. The datasets are downloaded as zip
file and extracted automatically, unless extraction is disabled by the
user.

## Usage

``` r
download_example_datasets(
  dst_dir = NULL,
  extract = TRUE,
  persistent = NULL,
  overwrite = FALSE,
  silent = FALSE
)
```

## Arguments

- dst_dir:

  The destination directory where the downloaded datasets will be
  stored. If NULL, the function will return the path to the cached zip
  file.

- extract:

  Logical. If TRUE, the downloaded zip file will be extracted.

- persistent:

  Logical. If TRUE, the downloaded datasets will be cached at
  [`datadir_persistent()`](https://spang-lab.github.io/metabodecon/reference/datadir_persistent.md)
  to speed up future calls to `download_example_datasets()`. If FALSE,
  the datasets will be cached at
  [`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md).
  If NULL, the function will check both paths for the cached datasets
  but will return
  [`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md)
  if the cached file does not yet exist.

- overwrite:

  Logical. If TRUE, existing files with the same name in the destination
  directory will be overwritten.

- silent:

  Logical. If TRUE, no output will be printed to the console.

## Value

The path to the downloaded (and possibly extracted) datasets.

## See also

[`datadir()`](https://spang-lab.github.io/metabodecon/reference/datadir.md)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
if (interactive()) {
     zip <- download_example_datasets(extract = FALSE, persistent = FALSE)
     dir <- download_example_datasets(extract = TRUE)
}
```
