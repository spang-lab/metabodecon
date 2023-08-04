# metabodecon <img src="man/figures/logo.svg" alt="man/figures/logo.svg" align="right" height="138" />

A package for deconvolution of 1D NMR spectra using `MetaboDecon1D` followed by alignment of the data by the `speaq` package and some additional postprocessing resulting in an aligned data matrix of aligned signal integrals.

- [Installation](#installation)
- [Usage](#usage)
- [Contribute](#contribute)
- [Open Issues](#open-issues)

## Installation

Copy paste the following command in a running R session (e.g. in RStudio):

```R
if (!"devtools" %in% installed.packages()[, "Package"]) {
    install.packages(
        pkgs = "devtools",
        repos = c(CRAN = "https://cloud.r-project.org")
    )
}
devtools::install_gitlab(
    repo = "grw28475/metabodecon",
    host = "https://gitlab.spang-lab.de/",
    auth_token = "glpat-ndxyfy5Ty7yksgy9MAFs",
    build_manual = TRUE,
    build_vignettes = TRUE
)
```

## Usage

The file [vignettes/metabodecon.Rmd](vignettes/metabodecon.Rmd) shows an example how `metabodecon` can be used to deconvolute an existing data set, followed by alignment of the data and some additional postprocessing steps, resulting in a data matrix of aligned signal integrals. To get a nice rendering of the document, it is highly recommended to download the repository and then open `docs/index.html` in a browser of your choice. Then you will get a view like this:

![man/figures/pkgdown_preview.png](man/figures/pkgdown_preview.png)

## Contribute

See file [CONTRIBUTE.md](CONTRIBUTE.md) for instructions on how to contribute to this project.

## Open Issues

1. In function `generate_lorentz_curves` it should be possible to specify the parameters via function arguments (this makes batch execution and/or repeated execution a lot easier).
2. Input prompting in function `generate_lorentz_curves` is broken on my Windows 11 system with R4.2.3. This should be checked.
3. Timestamps in `generate_lorentz_curves` would be nice, so you can see how long it will take for the function to finish.
