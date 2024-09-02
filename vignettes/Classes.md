---
title: "Classes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Classes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Metabodecon introduces a set of different S3 classes with corresponding methods to handle the data. The main classes are:

- `rawSpectrum`
- `rawSpectra`
- `deconvolutedSpectrum`
- `deconvolutedSpectra`
- `alignedSpectrum`
- `alignedSpectra`

The data flow is as follows:

1. To get an `alignedSpectra` object, you need to pass a `deconvolutedSpectra` object to `align_spectra()`
2. To get a `deconvolutedSpectra` object, you need a `glcSpectra` object

3. Create one or more spectrum objects using `read_spectrum()`, `simulate_spectrum()`, or `make_spectrum()`.
4. Create a spectra object using `make_spectra()` or `read_spectra()`.