# Make transparent

Make a color transparent by adding an alpha channel.

## Usage

``` r
transp(col = "violet", alpha = 0.08)
```

## Arguments

- col:

  Character string specifying the color to make transparent.

- alpha:

  Numeric value between 0 and 1 specifying the transparency level.

## Value

A character string representing the color with an alpha channel.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
transp("violet", 0.08)
#> [1] "#EE82EE14"
transp("black", 0.5)
#> [1] "#00000080"
```
