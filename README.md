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

## Contribute

Things you can update, are:

1. Function code in folder [R](R)
2. Function documentation in folder [R](R)
3. Package documentation in folder `vignettes`
4. Test cases in folder [tests](tests)
5. Dependencies in file [DESCRIPTION](DESCRIPTION)
6. Authors in file [DESCRIPTION](DESCRIPTION)

Whenever you update any of those things, you should run the below commands to check that everything is still working as expected

```R
devtools::test() # Execute tests from tests folder
devtools::document() # Build files in man folder
devtools::check() # Check package formalities
devtools::install() # Install as required by next command
pkgdown::build_site() # Build website in docs folder
```

After doing these steps, you can push your changes to Github and then use the following commands to release the package to CRAN:

```R
devtools::check(remote = TRUE, manual = TRUE)# Slower, but more realistic test than plain devtools::check()
devtools::spell_check() # Check spelling. Add false positives to inst/WORDLIST
revdepcheck::revdep_check(num_workers = 8) # Reverse dependency check
# See https://r-pkgs.org/release.html#sec-release-revdep-checks for details
devtools::submit_cran() # Submits the package to CRAN
```

Above steps are based on: <https://r-pkgs.org/release.html>
