# Working with the Rust bindings

1. Install Rust and Cargo via [rustup](https://www.rust-lang.org/tools/install)
2. Install [rextendr](https://github.com/extendr/rextendr)
3. Run `rextendr::document()` and `devtools::load_all()`

## Important

When editing the Rust code, `rextendr::document()` can be used to generate the wrapper functions.
I currently don't have the time to look into how to disable automatic exports of Rust struct
bindings, so NAMESPACE must be edited by hand to remove the export statements.

## Rust dependency

The Metabodecon Rust crate is currently a git dependency, meaning that there can sometimes be issues
with having an outdated local copy of the crate. Running `rextendr::clean()` or manually deleting
the build artifacts should usually make it download the latest version. If it doesn't, try deleting
the `Cargo.lock` file in `src/rust`. Normally, this could be fixed by running `cargo update`,
however, this will default to the virtual manifest in the root of the project, while rextendr uses a
different manifest in the `src/rust` directory.

# Examples

## Using the wrapper functions

```R
spectra <- read_spectra("misc/example_datasets/bruker/blood", "bruker", 10, 10)
rust_deconvolution <- deconvolute_rust(spectra[[1]], sfr = c(-2.2, 11.8), nfit = 10, smopts = c(2, 5), delta = 6.4, ignore_regions = c(4.7, 4.9), parallel = TRUE, optimize_settings = FALSE)
rust_deconvolutions <- multi_deconvolute_rust(spectra, sfr = c(-2.2, 11.8), nfit = 10, smopts = c(2, 5), delta = 6.4, ignore_regions = c(4.7, 4.9), parallel = TRUE, optimize_settings = FALSE)
```

## Using the structs directly

```R
# Read the spectra
spectra <- read_spectra("misc/example_datasets/bruker/blood", "bruker", 10, 10)
spectra <- lapply(spectra, function(s) Spectrum$new(s$cs, s$si, c(-2.2, 11.8)))

# Configure the Deconvoluter
deconvoluter <- Deconvoluter$new()
deconvoluter$set_moving_average_smoother(4, 3)
deconvoluter$add_ignore_region(4.7, 4.9)

# Deconvolute the spectra
deconvolutions <- deconvoluter$deconvolute_spectra(spectra)
deconvolutions <- deconvoluter$par_deconvolute_spectra(spectra)

# Getting the Lorentzian parameters
lorentzians <- lapply(deconvolutions, function(d) d$lorentzians())

# Compute the superposition of the Lorentzians for the first spectrum
superposition_internal <- deconvolutions[[1]]$par_superposition_vec(spectra[[1]]$chemical_shifts())

# Alternative method
A <- lorentzians[[1]]$A
lambda <- lorentzians[[1]]$lambda
x0 <- lorentzians[[1]]$x0
superposition_parameters <- Lorentzian$par_superposition_vec(spectra[[1]]$chemical_shifts(), A, lambda, x0)
```
