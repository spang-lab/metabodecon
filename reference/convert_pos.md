# Convert from unit A to unit B

Converts positions/widths from unit A to unit B. If the direction of
units A and B is reversed, the width's sign will be reversed as well. To
keep widths strictly positive, wrap the result with
[`abs()`](https://rdrr.io/r/base/MathFun.html).

## Usage

``` r
convert_pos(xa, ya, yb)

convert_width(xa, ya, yb)
```

## Arguments

- xa:

  A numeric vector specifying widths/positions in unit A.

- ya, yb:

  A numeric vector specifying the positions of at least two points in
  unit A / unit B.

## Value

A numeric vector of values converted from unit A to unit B.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
ya <- c(244, 246, 248, 250, 252)
yb <- c(15, 10, 5, 0, -5)
convert_width(c(2, 4, 8), ya, yb)
#> [1]  -5 -10 -20
convert_pos(c(247, 249), ya, yb)
#> [1] 7.5 2.5
```
