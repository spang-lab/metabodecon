# Compile bindings locally

1. Install Rust and Cargo via [rustup](https://www.rust-lang.org/tools/install)
2. Install [rextendr](https://github.com/extendr/rextendr)
3. Run `rextendr::document()` (this will generate and modify some files I didn't commit yet)
4. Run `devtools::load_all()`
5. The file `R/extendr-wrappers.R` will contain the wrapper functions.

# Example

## Using the bindings

```R
rextendr::document()
devtools::load_all()

# Read the spectra
spectra <- Spectrum$from_bruker_set("misc/example_datasets/bruker/blood", 10, 10, c(-2.2, 11.8))

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

## Using the wrapper functions

```R
rextendr::document()
devtools::load_all()

spectra <- metabodecon::read_spectra("misc/example_datasets/bruker/blood", 10, 10, c(-2.2, 11.8))
rust_deconvolution <- deconvolute_rust(spectra[[1]], sfr = c(-2.2, 11.8), nfit = 10, smopts = c(2, 5), delta = 6.4, ignore_regions = c(4.7, 4.9), parallel = TRUE, optimize_settings = FALSE)
rust_deconvolutions <- multi_deconvolute_rust(spectra, sfr = c(-2.2, 11.8), nfit = 10, smopts = c(2, 5), delta = 6.4, ignore_regions = c(4.7, 4.9), parallel = TRUE, optimize_settings = FALSE)
```
