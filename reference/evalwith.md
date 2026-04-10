# Evaluate an expression with predefined global state

Evaluates an expression with a predefined global state, including the:

- working directory (set via
  [`setwd()`](https://rdrr.io/r/base/getwd.html))

- global options (set via
  [`options()`](https://rdrr.io/r/base/options.html))

- graphical parameters (set via
  [`par()`](https://rdrr.io/r/graphics/par.html))

In addition to that, `evalwith` allows to:

- Redirect or capture the output and/or message stream via
  [`sink()`](https://rdrr.io/r/base/sink.html)

- Measure the runtime of the evaluated expression via
  [`system.time()`](https://rdrr.io/r/base/system.time.html)

- Creating a temporary test directory (inside
  [`tmpdir()`](https://spang-lab.github.io/metabodecon/reference/tmpdir.md))
  and populating it with input files according to `inputs`

- Predefine answers for calls to
  [`readline()`](https://rdrr.io/r/base/readline.html) happening during
  evaluation of `expr`

- Caching the result of the expression

All changes to the global state are reverted after the expression has
been evaluated.

## Usage

``` r
evalwith(
  expr,
  testdir = NULL,
  answers = NULL,
  output = NULL,
  message = NULL,
  plot = NULL,
  datadir_temp = c("default", "missing", "empty", "filled")[1],
  datadir_persistent = c("default", "missing", "empty", "filled")[1],
  inputs = character(),
  opts = NULL,
  pars = NULL,
  cache = FALSE,
  overwrite = FALSE
)
```

## Arguments

- expr:

  Expression to be evaluated.

- testdir:

  ID of the test directory. E.g. `"xyz/2"`. Will be created and
  populated with `inputs`. To clear, use `clear(testdir("xyz/2"))`.

- answers:

  Answers to be returned by readline().

- output:

  Path to the file where output stream should be redirected to. Use
  `"captured"` to capture the output.

- message:

  Path to the file where message stream be redirected to. Use
  `"captured"` to capture the messages.

- plot:

  An expression opening a device, the string "captured" or a path ending
  in ".pdf", ".svg", or ".png". Examples: `svg("tmp.svg")`,
  `quote(pdf("tmp.pdf"))`, `"captured"`, `"tmp.png"`. Passing
  `"captured"` is equivalent to passing `tempfile(fileext = ".png")`.

- datadir_temp:

  State of the mocked temporary data directory. See details section.

- datadir_persistent:

  State of the mocked persistent data directory. See details section.

- inputs:

  Paths to be copied to the test directory before evaluating `expr`.

- opts:

  Named list of options to be set. See
  [`options()`](https://rdrr.io/r/base/options.html).

- pars:

  Named list of parameters to be set. See
  [`par()`](https://rdrr.io/r/graphics/par.html).

- cache:

  Logical indicating whether to cache the result of the expression.

- overwrite:

  Logical indicating whether to overwrite the cache file if it already
  exists.

## Value

A list containing with following elements:

- `rv`: The return value of the expression.

- `runtime`: The "elapsed" runtime of the expression in seconds.
  Measured with
  [`system.time()`](https://rdrr.io/r/base/system.time.html).

- `output`: The captured output.

- `message`: The captured messages.

- `plot`: The path to the saved plot.

- `testdir`: The path to the test directory.

- `inputs`: The paths to the copied input files.

## Details

The `datadir_temp` and `datadir_persistent` arguments accept values
"missing", "filled" and "empty". Setting a value unequal NULL causes the
functions
[`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md)
and/or
[`datadir_persistent()`](https://spang-lab.github.io/metabodecon/reference/datadir_persistent.md)
to be replaced with mock functions pointing to fake directories.
Functions depending on these functions will then use the fake
directories instead of the real ones. When set to "missing" the returned
mock directory does not exist. When set to "empty" it exists and is
guaranteed to be empty. When set to "filled", it is populated with
example datasets.

Attention: the mocked functions, i.e.
[`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md)
and
[`datadir_persistent()`](https://spang-lab.github.io/metabodecon/reference/datadir_persistent.md)
cannot be used directly inside `expr` when called via
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html).
I'm not sure why, but it seems as if devtools and/or testthat have their
own copies of the functions which are used when the expression is
evaluated.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
x1 <- evalwith(output = "captured", cat("Helloworld\n"))
str(x1)
#> List of 7
#>  $ rv     : NULL
#>  $ runtime: num 0
#>  $ output : chr "Helloworld"
#>  $ message: chr(0) 
#>  $ plot   : NULL
#>  $ testdir: NULL
#>  $ inputs : chr(0) 

x2 <- evalwith(datadir_persistent = "missing", message = "captured", datadir())
#> Warning: /tmp/RtmpoVTN85/metabodecon/data does not exist. Please call `download_example_datasets()` first.
str(x2)
#> List of 7
#>  $ rv     : chr "/tmp/RtmpoVTN85/metabodecon/data"
#>  $ runtime: num 0.001
#>  $ output : chr(0) 
#>  $ message: chr(0) 
#>  $ plot   : NULL
#>  $ testdir: NULL
#>  $ inputs : chr(0) 

x3 <- evalwith(testdir = "dummy", inputs = "bruker/urine/urine_1", dir())
str(x3)
#> List of 7
#>  $ rv     : chr "urine_1"
#>  $ runtime: num 0
#>  $ output : chr(0) 
#>  $ message: chr(0) 
#>  $ plot   : NULL
#>  $ testdir: chr "dummy"
#>  $ inputs : chr "bruker/urine/urine_1"

x4 <- evalwith(Sys.sleep(0.02))
str(x4)
#> List of 7
#>  $ rv     : NULL
#>  $ runtime: num 0.02
#>  $ output : chr(0) 
#>  $ message: chr(0) 
#>  $ plot   : NULL
#>  $ testdir: NULL
#>  $ inputs : chr(0) 
```
