<!-- badges: start -->
[![R CMD check](https://github.com/spang-lab/metabodecon/workflows/R-CMD-check/badge.svg)](https://github.com/spang-lab/metabodecon/actions)
<!-- badges: end -->

# metabodecon <img src="man/figures/logo.svg" alt="man/figures/logo.svg" align="right" height="138" />

A package for deconvolution of 1D NMR spectra using `MetaboDecon1D` followed by alignment of the data using the `speaq` package and some additional post-processing, resulting in a data matrix of aligned signal integrals.

## Installation

Copy paste the following command in a running R session (e.g. in RStudio):

```R
if (!"devtools" %in% installed.packages()[, "Package"]) {
    install.packages("devtools", repos = c(CRAN = "https://cloud.r-project.org"))
}
devtools::install_github("spang-lab/metabodecon", build_manual = TRUE, build_vignettes = TRUE)
```

## Usage

At [Get Started](https://spang-lab.github.io/metabodecon/articles/metabodecon.html) you can see an example how `metabodecon` can be used to deconvolute an existing data set, followed by alignment of the data and some additional postprocessing steps, resulting in a data matrix of aligned signal integrals.

At [Function Reference](https://spang-lab.github.io/metabodecon/reference/index.html) you get an overview of all functions provided by metabodecon.

## Documentation

MetaboDecon1D's documentation is available at [spang-lab.github.io/MetaboDecon1D](https://spang-lab.github.io/MetaboDecon1D/). It includes pages about

- [Getting Started](https://spang-lab.github.io/metabodecon/articles/metabodecon.html)
- [Contribution Guidelines](https://spang-lab.github.io/MetaboDecon1D/articles/Contributing.html)
- [Function Reference](https://spang-lab.github.io/MetaboDecon1D/reference/index.html)
