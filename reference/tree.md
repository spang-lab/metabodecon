# Print the Structure of a Directory Tree

Prints the structure of a directory tree up to a specified maximum level
of depth. It lists all files and directories under the specified path,
displaying them in a tree-like structure.

## Usage

``` r
tree(
  path,
  max.level = 2,
  max.entries = Inf,
  show.counts = FALSE,
  files.first = FALSE,
  level = 0,
  prefix = ""
)

tree_preview(
  path,
  max.level = 1,
  max.entries = 9,
  show.counts = TRUE,
  files.first = TRUE
)
```

## Arguments

- path:

  The root path from which to start listing the directory structure.

- max.level:

  The maximum depth of directories to list.

- max.entries:

  Maximum number of children to print per directory. If a directory has
  more entries than this limit, the first `ceiling(max.entries / 2)` and
  last `floor(max.entries / 2)` children are shown, with `... [N]` in
  between, where `N` is the number of omitted children.

- show.counts:

  Logical. If `TRUE`, prints the number of direct children (files +
  subdirectories) in brackets after each directory name. Disabled by
  default.

- files.first:

  Logical. If `TRUE`, files are listed before subdirectories. If `FALSE`
  (default), subdirectories are listed first.

- level:

  Internal parameter used for recursion, indicating the current level of
  depth.

- prefix:

  Internal parameter used for formatting the printed tree structure.

## Value

NULL, called for its side effect of printing the directory structure.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
metabodecon_dir <- system.file(package = "metabodecon")
tree(metabodecon_dir, max.level = 1)
#> /home/runner/work/_temp/Library/metabodecon
#> ├── Meta/
#> ├── R/
#> ├── data/
#> ├── example_datasets/
#> ├── help/
#> ├── html/
#> ├── DESCRIPTION
#> ├── INDEX
#> ├── NAMESPACE
#> ├── NEWS.md
#> └── WORDLIST
```
