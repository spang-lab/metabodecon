---
title: "Purpose"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Purpose}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r knitr-setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)
```

Most statistical methods, such as linear models, statistical tests, clustering methods, or machine learning algorithms, expect a two-dimensional data structure where each row represents a sample and each column represents a feature.

In the case of metabolomic fingerprinting, the naive input would be a matrix with one row per tissue sample and one column per metabolite. However, depending on the measurement technology used, the desired data may not be directly measurable and must first be estimated from raw data using various methods.

The required methods depend on the technology used to measure the metabolites and the associated physical phenomena. The physical phenomena that influence the transformation from metabolite concentration to measured NMR signal are described in Listing 1 and graphically represented in Figure 1.

1. In NMR, metabolite concentrations are not measured directly; instead, concentrations of H1 nuclei in specific chemical environments are measured. A molecule can contain multiple H1 nuclei in different chemical environments, resulting in zero, one, or multiple signals per molecule. The chemical environment is determined by the ratio of magnetic field strength to the absorption frequency of the H1 nucleus. Instead of a vector (Mol1: 10, Mol2: 20, Mol3: 30), you might get a vector (CS1: 5, CS2: 10, CS3: 5, CS4: 20, CS5: 10, CS6: 10). This splitting of molecule concentrations into signal intensities can be understood as a mapping from R^m to R^n, where m is the number of molecules in the sample and n is the number of measurable chemical environments in the sample. We refer to this as molecule-signal mapping. The reverse transformation from signal intensities to molecule concentrations is called signal-molecule mapping.
2. Due to differences in temperature, magnetic field strength, etc., slightly different frequencies or magnetic field strengths are required to excite H1 nuclei in a specific chemical environment in different NMR measurements. This variation in energy between measurements is called signal shift. The reverse transformation, i.e., correcting the signal shift, is called signal alignment.
3. NMR spectra do not measure a sharp signal per chemical environment but rather a predefined number of signals within a specific energy range (essentially a "signal vector"). The problem, that all the signal vectors from the different chemical environments cover the same energy range is called signal overlap. The reverse transformation, i.e., separating the overlapped signals for different chemical environments, is called signal deconvolution.

Metabodecon provides a set of functions for reversing the phenomena described in points 2 and 3, in reverse order, i.e.:

1. Deconvoluting the overlapped signals
2. Aligning the deconvoluted signals

This results in a matrix of features where each row corresponds to a sample and each column to a concentration of H1 nuclei with a specific absorption energy.

Tools for solving the mapping of signal intensities to molecule concentrations include:

- TODO
- TODO
- TODO
