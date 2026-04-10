# Return path to metabodecon's data directory

Returns the path to the directory where
[`download_example_datasets()`](https://spang-lab.github.io/metabodecon/reference/download_example_datasets.md)
stores metabodecon's example data sets or any file within that
directory. By default this directory is a subdirectory of R's temporary
session directory. If `persistent` is set to `TRUE`, the directory
equals the data directory returned by
[`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) instead.

## Usage

``` r
datadir(file = NULL, warn = TRUE, persistent = NULL)
```

## Arguments

- file:

  Relative path to a file within the data directory.

- warn:

  Print a warning message when the requested path does not yet exist?

- persistent:

  Return the path to the persistent data directory instead of the
  temporary one?

## Value

Path to the data directory or a file within it.

## Details

The decision to use a temporary data dir as default and a persistent one
only optionally was made to conform to CRAN package policies, which
state that:

*Packages should not write in the user's home filespace (including*
*clipboards), nor anywhere else on the file system apart from the R*
*session's temporary directory \[...\] Limited exceptions may be
allowed* *in interactive sessions if the package obtains confirmation
from the* *user. For R version 4.0 or later \[...\] packages may store
user-specific* *data, configuration and cache files in their respective
user directories* *obtained from
[`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) \[...\].*

Source:
[cran.r-project.org/web/packages/policies](https://cran.r-project.org/web/packages/policies.html).

## See also

[`download_example_datasets()`](https://spang-lab.github.io/metabodecon/reference/download_example_datasets.md),
[`datadir_persistent()`](https://spang-lab.github.io/metabodecon/reference/datadir_persistent.md),
[`datadir_temp()`](https://spang-lab.github.io/metabodecon/reference/datadir_temp.md)

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
# Get temporary datadir and persistent datadir
datadir(persistent = FALSE, warn = FALSE)
#> [1] "/tmp/RtmpoVTN85/metabodecon/data"
datadir(persistent = TRUE,  warn = FALSE)
#> [1] "/home/runner/.local/share/R/metabodecon"

# Get persistent datadir if existing else temp datadir. Set `warn = TRUE`
# to raise a warning if none of the directories exist yet.
datadir(warn = FALSE)
#> [1] "/tmp/RtmpoVTN85/metabodecon/data"
if (interactive()) datadir()

# Get PERSISTENT_DATADIR/bruker if existing else TEMP_DATADIR/bruker
datadir(file = "bruker/urine", warn = FALSE)
#> [1] "/tmp/RtmpoVTN85/metabodecon/data/bruker/urine"
```
