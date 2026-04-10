# Show head and tail rows of a matrix-like object

Returns the first and last `n` rows of a matrix or data frame. If the
input has fewer than `2*n` rows, overlapping rows are returned only
once.

## Usage

``` r
headtail(x, n = 6)
```

## Arguments

- x:

  A matrix or data frame.

- n:

  Number of rows to take from the top and bottom.

## Value

A subset of `x` containing head and tail rows.

## Author

2024-2026 Tobias Schmidt: initial version.

## Examples

``` r
x <- matrix(seq_len(30), nrow = 10)
headtail(x, n = 2)
#>      [,1] [,2] [,3]
#> [1,]    1   11   21
#> [2,]    2   12   22
#> [3,]    9   19   29
#> [4,]   10   20   30
```
