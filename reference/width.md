# Calculate the Width of a Numeric Vector

Calculates the width of a numeric vector by computing the difference
between the maximum and minimum values in the vector.

## Usage

``` r
width(x)
```

## Arguments

- x:

  A numeric vector.

## Value

The width of the vector, calculated as the difference between its
maximum and minimum values.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
vec <- c(1, 3, 5, 7, 9)
width(vec)
#> [1] 8
```
