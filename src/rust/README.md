# Compile bindings locally

1. Install Rust and Cargo via [rustup](https://www.rust-lang.org/tools/install)
2. Install [rextendr](https://github.com/extendr/rextendr)
3. Run `rextendr::document()` (this will generate and modify some files I didn't commit yet)
4. The file `R/extendr-wrappers.R` will contain the wrapper functions.

# Examples

## Private functions

```R
spec <- metabodecon::read_spectrum("misc/example_datasets/bruker/blood/blood_01")
spec$sfr <- c(-2.2, 11.8)
lorentzians <- backend_deconvolute_spectrum(spec, c(4.7, 4.9))
lorentzians <- backend_par_deconvolute_spectrum(spec, c(4.7, 4.9))
```

## Struct bindings

```R
spectra <- Spectrum$from_bruker_set("misc/example_datasets/bruker/blood", 10, 10, c(-2.2, 11.8))
deconvoluter <- Deconvoluter$new()
deconvoluter$add_ignore_region(4.7, 4.9)
deconvolutions <- deconvoluter$deconvolute_spectra(spectra)
deconvolutions <- deconvoluter$par_deconvolute_spectra(spectra)
lorentzians <- lapply(deconvolutions, function(d) d$lorentzians())
```
