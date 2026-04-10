# S3 methods for mdm objects

**\[experimental\]**

**WARNING: These methods are experimental and must not be used in
production. Their API is very likely to change in
non-backwards-compatible ways over the next few weeks.**

S3 methods for objects of class `mdm` and `summary.mdm`.

`predict.mdm()` predicts probabilities, classes, link scores, or all
three from an `mdm` object. When `newdata` is a spectra object, the
spectra are deconvoluted, aligned and snapped to the reference stored in
the model before prediction. When `newdata` is a numeric matrix, it is
used directly as the feature matrix.

`print.mdm()` prints a compact model summary.

`coef.mdm()` returns lasso coefficients.

`plot.mdm()` plots the lasso path.

`summary.mdm()` builds a compact summary list.

`print.summary.mdm()` prints formatted output for `summary.mdm` objects.

## Usage

``` r
# S3 method for class 'mdm'
predict(
  object,
  newdata,
  type = c("all", "prob", "class", "link"),
  s = "lambda.min",
  nworkers = 1,
  verbosity = 1,
  ...
)

# S3 method for class 'mdm'
print(x, ...)

# S3 method for class 'mdm'
coef(object, ...)

# S3 method for class 'mdm'
plot(x, ...)

# S3 method for class 'mdm'
summary(object, ...)

# S3 method for class 'summary.mdm'
print(x, ...)
```

## Arguments

- object, x:

  A fitted `mdm` object (for `predict`, `coef`, `summary`, `print` and
  `plot`) or a `summary.mdm` object (for `print.summary.mdm`).

- newdata:

  Spectra object or numeric feature matrix.

- type:

  Prediction type, one of `"all"`, `"prob"`, `"class"`, `"link"`.

- s:

  Regularization value for lasso predictions.

- nworkers:

  Number of workers to deconvolute and align `newdata`.

- verbosity:

  Integer verbosity level.

- ...:

  Passed to underlying methods where applicable.

## Value

- `predict`: numeric vector of probabilities, classes, and/or link
  scores.

- `print`: invisibly returns `x`.

- `coef`: coefficient object from `glmnet`.

- `plot`: invisibly returns `NULL`.

- `summary`: object of class `summary.mdm`.

- `print.summary.mdm`: invisibly returns `x`.

## Examples

``` r
m <- structure(
  list(
    model = NULL,
    ref = NULL,
    meta = list(npmax = 1000, nfit = 3, smit = 2, smws = 5,
                delta = 6.4, maxShift = 100, maxCombine = 50,
                combineMethod = 2)
  ),
  class = "mdm"
)
print(m)
#> metabodecon model (mdm)
#>   npmax:         1000
#>   nfit:          3
#>   smit:          2
#>   smws:          5
#>   delta:         6.4
#>   maxShift:      100
#>   maxCombine:    50
#>   combineMethod: 2
summary(m)
#> Summary of mdm
#>   npmax:         1000
#>   nfit:          3
#>   smit:          2
#>   smws:          5
#>   delta:         6.4
#>   maxShift:      100
#>   maxCombine:    50
#>   combineMethod: 2
#>   n_peaks:       0
#>   grid_rows:     0

if (FALSE) { # \dontrun{
  m <- cv_mdm(spectra, y, sfr = c(11, -2))
  predict(m, test_spectra, type = "prob")
  coef(m)
  plot(m)
} # }
```
