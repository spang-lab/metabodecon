# Returns the names of a metabodecon collection object.

Returns the names of a metabodecon collection object.

## Usage

``` r
get_names(x, default = "spectrum_%d")
```

## Arguments

- x:

  A metabodecon collection object.

- default:

  Default names if no names are found. Passed on to `get_default_names`.

## Value

A character vector of names.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
s1 <- list()
s2 <- list(name = "foo")
s3 <- list(name = "foo", meta = list(name = "bar"))

get_names(list(s1, s1)) # c("spectrum_1", "spectrum_2")
#> [1] "spectrum_1" "spectrum_2"
get_names(list(s1, myspec = s1)) # c("spectrum_1", "myspec")
#> [1] "spectrum_1" "myspec"    
get_names(list(s1, myspec = s2)) # c("spectrum_1", "foo")
#> [1] "spectrum_1" "foo"       
get_names(list(s1, myspec = s3)) # c("spectrum_1", "bar")
#> [1] "spectrum_1" "bar"       
```
