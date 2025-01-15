<!-- badges: start -->
[![R CMD check](https://github.com/spang-lab/metabodecon/workflows/R-CMD-check/badge.svg)](https://github.com/spang-lab/metabodecon/actions)
<!-- badges: end -->

# metabodecon <img src="man/figures/logo.svg" alt="man/figures/logo.svg" align="right" height="138" />

A framework for deconvolution, alignment and postprocessing of 1D NMR spectra, resulting in a data matrix of aligned signal integrals. The deconvolution part uses the algorithm described in [Koh et al. (2009)](https://doi.org/10.1016/j.jmr.2009.09.003). The alignment part is based on functions from the 'speaq' package, described in [Beirnaert et al. (2018)](https://doi.org/10.1371/journal.pcbi.1006018) and [Vu et al. (2011)](https://doi.org/doi:10.1186/1471-2105-12-405). A detailed description and evaluation of an early version of the package, 'MetaboDecon1D v0.2.2', can be found in [Haeckl et al. (2021)](https://doi.org/doi:10.3390/metabo11070452).

## Installation

Copy paste the following command in a running R session (e.g. in RStudio):

```R
if (!"devtools" %in% installed.packages()[, "Package"]) {
    install.packages("devtools", repos = c(CRAN = "https://cloud.r-project.org"))
}
devtools::install_github("spang-lab/metabodecon", build_manual = TRUE, build_vignettes = TRUE)
```

## Usage

At [Get Started](https://spang-lab.github.io/metabodecon/articles/metabodecon.html) you can see an example how metabodecon can be used to deconvolute an existing data set, followed by alignment of the data and some additional postprocessing steps, resulting in a data matrix of aligned signal integrals.

At [Function Reference](https://spang-lab.github.io/metabodecon/reference/index.html) you get an overview of all functions provided by metabodecon.

## Documentation

metabodecon's documentation is available at [spang-lab.github.io/metabodecon](https://spang-lab.github.io/metabodecon/). It includes pages about

- [Getting Started](https://spang-lab.github.io/metabodecon/articles/metabodecon.html)
- [Contribution Guidelines](https://spang-lab.github.io/metabodecon/articles/Contributing.html)
- [Function Reference](https://spang-lab.github.io/metabodecon/reference/index.html)
