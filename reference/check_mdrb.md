# Check Rust Backend Requirements

`check_mdrb()` returns a boolean indicating whether a suitable version
of the metabodecon Rust backend
[mdrb](https://github.com/spang-lab/mdrb) is currently installed.

`check_mdrb_deps()` returns a list with information about the
installation status of mdrb system dependencies.

## Usage

``` r
check_mdrb(stop_on_fail = FALSE)

check_mdrb_deps(verbose = FALSE)
```

## Arguments

- stop_on_fail:

  If TRUE, an error is thrown if the check fails, providing instructions
  on how to install or upgrade mdrb.

- verbose:

  If TRUE, additional information is printed during the check process.

## Value

`check_mdrb()` returns TRUE if a suitable version of mdrb is installed,
else FALSE.

`check_mdrb_deps()` returns a data.frame as follows:

                check           passed  comment
        r       R >= 4.2        TRUE    Current: R 4.4.2
        rtools  Rtools exist    TRUE    Tested with: pkgbuild::has_build_tools()
        cargo   cargo >= 1.80   TRUE    Current: cargo 1.84.1 (66221abde 2024-11-19)
        rustc   rustc >= 1.80   TRUE    Current: rustc 1.84.1 (e71f9a9a9 2025-01-27)

Column `check` is a string describing the performed check.  
Column `passed` is a boolean indicating whether the check passed.  
Column `comment` is a string string describing the check result.

The rownames of the dataframe one-word descriptions of the performed
checks.

## Author

2024-2025 Tobias Schmidt: initial version.

## Examples

``` r
check_mdrb()
#> [1] TRUE

# \donttest{
# Checking dependencies might take more than 5 seconds, as it
# requires the compilation of a small test program as well as
# running `cargo --version` and `rustc --version`, which,
# depending on your system, might involve updating or installing
# Rust toolchain components.
check_mdrb_deps(verbose = TRUE)
#> 2026-04-10 01:12:09.12 Checking R version...
#> 2026-04-10 01:12:09.12 Checking if buildtools exist...
#> Trying to compile a simple C file
#> Running /opt/R/4.5.3/lib/R/bin/R CMD SHLIB foo.c
#> using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0’
#> gcc -std=gnu2x -I"/opt/R/4.5.3/lib/R/include" -DNDEBUG -DR_NO_REMAP  -I/usr/local/include    -fpic  -g -O2  -c foo.c -o foo.o
#> gcc -std=gnu2x -shared -L/opt/R/4.5.3/lib/R/lib -L/usr/local/lib -o foo.so foo.o -L/opt/R/4.5.3/lib/R/lib -lR
#> 
#> 2026-04-10 01:12:09.30 Checking cargo version...
#> 2026-04-10 01:12:09.31 Checking rustc version...
#> 2026-04-10 01:12:09.33 Done
#>                check passed                                      comment
#> r           R >= 4.2   TRUE                             Current: R 4.5.3
#> rtools  Rtools exist   TRUE        Testcall: pkgbuild::has_build_tools()
#> cargo  cargo >= 1.80   TRUE Current: cargo 1.94.1 (29ea6fb6a 2026-03-24)
#> rustc  rustc >= 1.80   TRUE Current: rustc 1.94.1 (e408947bf 2026-03-25)
# }
```
