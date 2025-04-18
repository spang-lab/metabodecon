<!-- badges: start -->
[![R-CMD-check](https://github.com/spang-lab/metabodecon/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/spang-lab/metabodecon/actions/workflows/R-CMD-check.yaml)
[![Install-Check](https://github.com/spang-lab/metabodecon/actions/workflows/test-install.yaml/badge.svg)](https://github.com/spang-lab/metabodecon/actions/workflows/test-install.yaml)
[![Codecov test coverage](https://codecov.io/gh/spang-lab/metabodecon/branch/main/graph/badge.svg)](https://app.codecov.io/gh/spang-lab/metabodecon?branch=main)
[![GitHub version](https://img.shields.io/github/v/release/spang-lab/metabodecon?label=GitHub&color=blue)](https://github.com/spang-lab/metabodecon/releases)
[![CRAN version](https://img.shields.io/cran/v/metabodecon?label=CRAN&color=blue)](https://cran.r-project.org/package=metabodecon)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/metabodecon)](https://cranlogs.r-pkg.org/badges/grand-total/metabodecon)
<!-- badges: end -->

# metabodecon <img src="man/figures/logo.svg" alt="man/figures/logo.svg" align="right" height="138" />

A framework for deconvolution, alignment and postprocessing of 1D NMR spectra, resulting in a data matrix of aligned signal integrals. The deconvolution part uses the algorithm described in [Koh et al. (2009)](https://doi.org/10.1016/j.jmr.2009.09.003). The alignment part is based on functions from the 'speaq' package, described in [Beirnaert et al. (2018)](https://doi.org/10.1371/journal.pcbi.1006018) and [Vu et al. (2011)](https://doi.org/doi:10.1186/1471-2105-12-405). A detailed description and evaluation of an early version of the package, 'MetaboDecon1D v0.2.2', can be found in [Haeckl et al. (2021)](https://doi.org/doi:10.3390/metabo11070452).

## Installation

To install the **stable version** from [CRAN](https://cran.r-project.org/), including all [Bioconductor](https://www.bioconductor.org/) dependencies, paste the following commands in a running R session (e.g. in RStudio):

```R
install.packages("pak")
pak::pkg_install("metabodecon")
```

Alternatively, if you prefer installing via the traditional `install.packages()` function, you can do so by running the following commands:

```R
# Install Bioconductor dependencies
install.packages("BiocManager")
BiocManager::install(c("MassSpecWavelet", "impute"))

# Install metabodecon
install.packages("metabodecon")
```

To install the **development version** from [GitHub](https://github.com/spang-lab/metabodecon/) use:

```R
install.packages("pak")
pak::pkg_install("spang-lab/metabodecon")
```

## Usage

At [Get Started](https://spang-lab.github.io/metabodecon/articles/metabodecon.html) you can see an example how metabodecon can be used to deconvolute an existing data set, followed by alignment of the data and some additional postprocessing steps, resulting in a data matrix of aligned signal integrals.

At [Function Reference](https://spang-lab.github.io/metabodecon/reference/index.html) you get an overview of all functions provided by metabodecon.

## Documentation

metabodecon's documentation is available at [spang-lab.github.io/metabodecon](https://spang-lab.github.io/metabodecon/). It includes pages about

- [Getting Started](https://spang-lab.github.io/metabodecon/articles/metabodecon.html)
- [Contribution Guidelines](https://spang-lab.github.io/metabodecon/articles/Contributing.html)
- [Function Reference](https://spang-lab.github.io/metabodecon/reference/index.html)
