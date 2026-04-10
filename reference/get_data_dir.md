# Retrieve directory path of an example dataset

Returns the path to the directory storing the example files shipped with
metabodecon.

Deprecated since metabodecon v1.2.0. Please use
[`datadir()`](https://spang-lab.github.io/metabodecon/reference/datadir.md)
instead. See examples below for usage.

**\[deprecated\]**

## Usage

``` r
get_data_dir(
  dataset_name = c("", "blood", "test", "urine", "aki"),
  warn = TRUE
)
```

## Arguments

- dataset_name:

  Either `""`, `"test"`, `"blood"`, `"urine"` or `"aki"`.

- warn:

  Whether to print a warning message when the example folders do not yet
  exist, i.e.
  [`download_example_datasets()`](https://spang-lab.github.io/metabodecon/reference/download_example_datasets.md)
  has not been called yet.

## Value

Path to the directory storing the example files.

## See also

[`download_example_datasets()`](https://spang-lab.github.io/metabodecon/reference/download_example_datasets.md)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
x <- get_data_dir("urine")                     # Deprecated
#> Warning: /tmp/RtmpoVTN85/metabodecon/data does not exist. Please call `download_example_datasets()` first.
#> Warning: /tmp/RtmpoVTN85/metabodecon/data/example_datasets/bruker/urine does not exist. Please call `download_example_datasets(extract = TRUE)` first.
y <- datadir("example_datasets/bruker/urine")  # Preferred
#> Warning: /tmp/RtmpoVTN85/metabodecon/data/example_datasets/bruker/urine does not exist. Please call `download_example_datasets()` first.
cat(x, y, sep = "\n")
#> /tmp/RtmpoVTN85/metabodecon/data/example_datasets/bruker/urine
#> /tmp/RtmpoVTN85/metabodecon/data/example_datasets/bruker/urine
```
