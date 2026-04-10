# Return Path to File or Directory in metabodecon Package

Recursively searches for files or directories within the 'metabodecon'
package that match the given name.

## Usage

``` r
metabodecon_file(name = "sim_01")
```

## Arguments

- name:

  The name to search for.

## Value

The file or directory path.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
# Unambiguous paths
metabodecon_file("urine_1")
#> [1] "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/urine/urine_1"
metabodecon_file("urine_1.dx")
#> [1] "/home/runner/work/_temp/Library/metabodecon/example_datasets/jcampdx/urine/urine_1.dx"
metabodecon_file("sim/sim_01")
#> [1] "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/sim/sim_01"

# Ambiguous paths (i.e. multiple matches)
metabodecon_file("sim")
#> [1] "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/sim"
metabodecon_file("urine")
#> [1] "/home/runner/work/_temp/Library/metabodecon/example_datasets/bruker/urine" 
#> [2] "/home/runner/work/_temp/Library/metabodecon/example_datasets/jcampdx/urine"

# Non-existing paths (i.e. a character vector of length zero gets returned)
metabodecon_file("asdfasdf")
#> character(0)
```
