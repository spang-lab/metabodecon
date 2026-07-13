# Install Rust Backend

Installs metabodecon's Rust backend
[mdrb](https://github.com/spang-lab/mdrb) from
[R-Universe](https://spang-lab.r-universe.dev/mdrb).

**\[experimental\]**

The Rust backend is entirely optional; metabodecon's pure-R backend is
the default and always available. R-Universe provides pre-built `mdrb`
binaries only for the **two most recent R releases**. On those,
installation is a plain binary download and requires no toolchain. On
older R versions (still `>= 4.2`) or platforms without a pre-built
binary, `mdrb` must be built from source, which requires a Rust
toolchain (`cargo` and `rustc` `>= 1.80`; check with
[`check_mdrb_deps()`](https://spang-lab.github.io/metabodecon/reference/check_mdrb.md)).
If automatic installation fails, `install_mdrb()` does not error: it
prints guidance and returns `FALSE`, pointing to
<https://github.com/spang-lab/mdrb> for manual installation.

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

Invisibly returns `TRUE` if `mdrb` is available after the call, else
`FALSE`. Called mainly for the side effect of installing the Rust
backend.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
if (interactive()) try(install_mdrb())
```
