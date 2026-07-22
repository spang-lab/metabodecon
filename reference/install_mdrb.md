# Install Rust Backend

Installs metabodecon's Rust backend
[mdrb](https://github.com/spang-lab/mdrb) from
[R-Universe](https://spang-lab.r-universe.dev/mdrb).

lifecycle::badge("experimental")

## Usage

``` r
install_mdrb(ask = TRUE, ...)
```

## Arguments

- ask:

  Whether to ask for confirmation before attempting installation.
  Default is TRUE.

- ...:

  Additional arguments to pass to
  [`utils::install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  when attempting installation of mdrb.

## Value

NULL. Called for side effect of installing the Rust backend.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
if (interactive()) try(install_mdrb())
```
