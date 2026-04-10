# Read one or more spectra from Disk

`read_spectrum()` reads a single spectrum from disk and returns it as
`spectrum` object.
[`read_spectra()`](https://spang-lab.github.io/metabodecon/reference/read_spectra.md)
can be used to read multiple spectra at once and returns a `spectra`
object.

## Usage

``` r
read_spectrum(
  data_path = metabodecon_file("bruker/sim/sim_01"),
  file_format = "bruker",
  expno = 10,
  procno = 10,
  raw = FALSE,
  silent = TRUE,
  force = FALSE
)
```

## Arguments

- data_path:

  The path of the file/folder containing the spectrum data. E.g.
  `"example_datasets/jcampdx/urine/urine_1.dx"` or
  `"example_datasets/bruker/urine/urine"`.

- file_format:

  The file_format of the spectrum file. E.g. `"bruker"` or `"jcampdx"`.

- expno, procno:

  The experiment/processing number for the file. E.g. `"10"`. Only
  relevant if `file_format` equals `"bruker"`. For details see section
  [File
  Structure](https://spang-lab.github.io/metabodecon/articles/FAQ.html#file-structure)
  in the metabodecon FAQ.

- raw:

  If `FALSE`, scales the returned signal intensities based on
  information available in the spectrum metadata, in particular
  `NC_proc`. For details see `processing-reference.pdf`, available at
  <https://www.bruker.com/en.html> at section 'Services & Support \>
  Documentation & Manuals \> Magnetic Resonance \> Acquisition &
  Processing \> TopSpin Processing Commands and Parameters' (requires
  login).

- silent:

  If `TRUE`, no output will be printed to the console.

- force:

  If `TRUE`, try to continue when encountering errors and print info
  messages instead. To hide these messages as well, set `silent = TRUE`.

## Value

A `spectrum` object as described in [Metabodecon
Classes](https://spang-lab.github.io/metabodecon/articles/Classes.html).

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
relpath <- "example_datasets/bruker/urine"
urine <- system.file(relpath, package = "metabodecon")
urine_1 <- file.path(urine, "urine_1")
urine_2 <- file.path(urine, "urine_2")
x1 <- read_spectrum(urine_1)
x2 <- read_spectrum(urine_2)
xx <- read_spectra(urine)
str(xx)
#> List of 2
#>  $ urine_1:List of 3
#>   ..$ si  : num [1:131072] 316.2 250.8 26.2 -234.2 -265.5 ...
#>   ..$ cs  : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...
#>   ..$ meta:List of 6
#>   .. ..$ fq    : num [1:131072] 6e+08 6e+08 6e+08 6e+08 6e+08 ...
#>   .. ..$ name  : chr "urine_1"
#>   .. ..$ path  : chr "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/urine/urine_1"
#>   .. ..$ type  : NULL
#>   .. ..$ simpar: NULL
#>   .. ..$ mfs   : NULL
#>   ..- attr(*, "class")= chr "spectrum"
#>  $ urine_2:List of 3
#>   ..$ si  : num [1:131072] -1544 -1464 -1416 -1436 -1398 ...
#>   ..$ cs  : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...
#>   ..$ meta:List of 6
#>   .. ..$ fq    : num [1:131072] 6e+08 6e+08 6e+08 6e+08 6e+08 ...
#>   .. ..$ name  : chr "urine_2"
#>   .. ..$ path  : chr "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/urine/urine_2"
#>   .. ..$ type  : NULL
#>   .. ..$ simpar: NULL
#>   .. ..$ mfs   : NULL
#>   ..- attr(*, "class")= chr "spectrum"
#>  - attr(*, "class")= chr "spectra"
str(x1)
#> List of 3
#>  $ si  : num [1:131072] 316.2 250.8 26.2 -234.2 -265.5 ...
#>  $ cs  : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ...
#>  $ meta:List of 6
#>   ..$ fq    : num [1:131072] 6e+08 6e+08 6e+08 6e+08 6e+08 ...
#>   ..$ name  : chr "urine_1"
#>   ..$ path  : chr "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/urine/urine_1"
#>   ..$ type  : NULL
#>   ..$ simpar: NULL
#>   ..$ mfs   : NULL
#>  - attr(*, "class")= chr "spectrum"
stopifnot(all.equal(x1, xx$urine_1))
```
