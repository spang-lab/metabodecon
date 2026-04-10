# Get URL of Metabodecon "Get Started" Page

`get_started` and `aaa_Get_Started` both return (and optionally open)
the URL of the "Get Started" page of the metabodecon documentation. The
`aaa_Get_Started` version exists, because functions are listed
alphabetically in the reference manual and we want `get_started` to be
shown at the top of the list (i.e., it needs to start with an 'a').

## Usage

``` r
aaa_Get_Started(open_browser = interactive())

get_started(open_browser = interactive())
```

## Arguments

- open_browser:

  If TRUE, the "Get Stated" page is opened in the default browser.

## Value

A character string containing the URL of the "Get Started" page.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
get_started(open_browser = FALSE)
#> [1] "https://spang-lab.github.io/metabodecon/articles/Get_Started.html"
get_started()
#> [1] "https://spang-lab.github.io/metabodecon/articles/Get_Started.html"
```
