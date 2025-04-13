# metabodecon 1.4.2

- Fixed `aaa_Get_Started` entry in Manual
- Changed default value of argument `verbose` from FALSE to TRUE for function `align()` and `deconvolute()`.
- Added argument `install_deps` to `align()`. If non-CRAN dependencies required by `align()` are missing and `install_deps` is TRUE, these dependencies are now installed automatically. If `install_deps` is NULL (default), the user is asked interactively for confirmation before attempting the install.

# metabodecon 1.4.1

- Added `get_started()` and `metabodecon-package` to manual
- Improved `plot_spectrum()` default margins.
- Improved plots shown during `deconvolute()`: SFW and WSHW are now both shown as rectangles instead of lines.
- Improved `install_mdrb()` example
- Improved `Get_Started` article
- Included `Get_Started` article as vignette within the package

# metabodecon 1.4.0

- Improved Github Workflow (GWF) to test installation on a clean Windows/Linux/Mac OS with R pre-installed, but without R-tools and any packages.
- Added `use_rust` option to `deconvolute()`. If `use_rust` is TRUE, the deconvolution is done using the implementation from Rust package [metabodecon-rust](https://github.com/SombkeMaximilian/metabodecon-rust/tree/main). Using the Rust backend requires R package [mdrb](https://github.com/spang-lab/mdrb) (Metabodecon Rust Backend) to be installed first. For this purpose the following additional functions are provided:
    - `install_mdrb()`: Installs mdrb
    - `check_mdrb()`: Checks whether a suitable version of mdrb is already installed
    - `check_mdrb_deps()`: Checks whether all required system dependencies of mdrb are installed

# metabodecon 1.3.0

- Added Github Workflow (GWF) to test installation on a clean Windows/Linux/Mac OS with R pre-installed, but without R-tools and any packages. Closes Todo [Test Install on clean OS].
- Fixed GWF for testing code coverage script
- Improved formatting for R-CMD-check-GWF and pkgdown-GWF
- Improved `align()`. The new implementation is faster and returns more information. In particular, the chemical shifts of the aligned peaks centers as well as the superposition of the aligned peaks are returned.
- Improved documentation and defaults for `speaq_align()` and `combine_peaks()`.
- Improved `draw_spectrum()`:
    - Added parameters `bt_text`, `lt_text`, `tp_text` and `rt_text` to `plot_spectrum()` to allow for full control over the text labels at the plot margins.
    - Added parameter `sf_vert` to `plot_spectrum()` to allow configuration of the height of the vertical lines drawn at the peak centers.
    - Added the option to fill the area under lorentzian curves with color.
    - Improved the legend of the plot.

[Test Install on clean OS]: TODOS.md#test-install-on-clean-os

# metabodecon 1.2.6

- Fixed a bug in `MetaboDecon1D()` that caused argument `file_path` to be interpreted as a relative path, even if it was an absolute path.

# metabodecon 1.2.5

- Fixed a bug in `read_spectrum()` that caused argument `raw` to not be passed on to `read_jcampdx()`.

# metabodecon 1.2.4

- Documentation updates

# metabodecon 1.2.3

- Documentation updates

# metabodecon 1.2.2

- Documentation updates

# metabodecon 1.2.1

- Documentation updates

# metabodecon 1.2.0

Finished the following tasks. For details about each task, see
[TODOS.md](https://github.com/spang-lab/metabodecon/blob/main/TODOS.md).

- CRAN-0: Omit "Functions for" in title
- CRAN-1: Omit "Functions for" in DESCRIPTION
- CRAN-2: Explain acronyms like NMR
- CRAN-3: Use correct reference format in DESCRIPTION
- CRAN-4: Explain return value in function docs
- CRAN-5: Remove examples from unexported functions
- CRAN-6: Fix vignettes
- CRAN-7: Check dontrun examples
- CRAN-8: Functions should not write to disk by default
- CRAN-9: Functions should not change working dir or global options
- FEATURE-01: Use temp dirs for example data
- FEATURE-02: Add minimal example dataset
- FEATURE-03: Batch Mode
- FEATURE-04: Parallelize
- FEATURE-05: Add test suite
- FEATURE-06: Return lambda in hertz
- FEATURE-07: Improve return value
- FEATURE-09 Implement `read_spectra()`
- FEATURE-11: Accept dataframes in GLC
- FEATURE-14: Provide simulated datasets
- FEATURE-15: Add lifecycle badges
- FEATURE-16: Improve multiprocessing
- FEATURE-17: Discard output
- FEATURE-18: Implement `plot_spectrum()`
- FEATURE-20: Implement `deconvolute_blood()`
- FIX-1: Prevent crashes for high smoothing
- REFACTOR-01: Combine load_spectrum functions
- REFACTOR-02: Improve Text Output (`-License`, `+Timestamps`)
- REFACTOR-04: Plotting speed
- REFACTOR-05: Speedup smoothing
- REFACTOR-06: Use a single unit as source of truth
- REFACTOR-07: Split monolithic functions into smaller parts
- REFACTOR-08: Improve docs for Metabodecon1D return value
- REFACTOR-09: Replace glc with `generate_lorentz_curves()`
- REFACTOR-10: Replace all md1d with `MetaboDecon1D()` calls
- REFACTOR-11: Implement `calc_prarp()`
- REFACTOR-12: Write compliance tests
- REFACTOR-13: Write PRARP tests

# metabodecon 1.1.1

API:

* Fixed a bug in `generate_lorentz_curves()` that caused the function to always use file format "bruker", even when file format "jcampdx" was specified.

Datasets:

* Fixed filenames of samples in blood dataset (renamed from `Bood_<nr>` to `blood_<nr>`).
* Renamed `example_datasets/jcampdx/urine/urine.dx` to `example_datasets/jcampdx/urine/urine_1.dx` and renamed `example_datasets/bruker/urine/urine/` to `example_datasets/bruker/urine/urine_1/`. This was done because `list.files` seems to return different orderings for `urine.dx` and `urine_2.dx` in different operating systems, whereas `urine_1.dx` and `urine_2.dx` are sorted the same way everywhere. This makes it easier to write clear and concise test cases, because we don't need to check for file ordering.

Documentation:

* Fixed broken image in [vignettes/FAQ.Rmd](https://github.com/spang-lab/metabodecon/blob/main/vignettes/FAQ.Rmd).

Testing:

* Added unit tests for `generate_lorentz_curves()`.
* Enabled parallel processing for unit tests.
* Created initial versions of `tests/testthat/test-generate_lorentz_curves-[1-4].R`.
* Added `generate_lorentz_curves_v2()` to `DESCRIPTION/Config/testthat/start-first`.
* Adjusted existing tests to use the updated version of `example_datasets` (sample `urine` was renamed to `urine_1`, as mentioned in above in section *Datasets*)

Internal:

* Added functions `%||%`, `msg()` and `msgf` to `R/util.R`.
* Added elements `range_water_signal_ppm` and `signal_free_region` to returned list of function `deconvolute_spectrum`.
* Function `with` now prints error messages to stderr even if the message stream is redirected.
* Copied function `deconvolution()` from `R/MetaboDecon1D.R` to `R/main_v2.R` as `.deconvolute_spectrum`.
* Fixed order of params in `deconvolution`.
* Fixed `download_example_datasets()`. Argument `overwrite` is passed correctly on to `cache_example_datasets()`.
* Changed URL of example datasets `xds$url` from `https://github.com/spang-lab/metabodecon/releases/download/v1.0.2/example_datasets.zip` to `https://github.com/spang-lab/metabodecon/releases/download/v1.1.0/example_datasets.zip`.
* Improved `cache_example_datasets()`. Extraction now only is done if `extract == TRUE` AND the resulting folder does not yet exist (saves approx. 2-3s on each call). To overwrite a possible existing folder, argument `overwrite` can be set to TRUE.
* Fixed formatting of `test_helpers.R`
* Added linter config `.lintr`

# metabodecon 1.1.0

API:

* Improved function `download_example_datasets()` by adding caching and making it more stable
* Replaced function `get_data_dir()` with `datadir()` and its helper functions `datadir_persistent()`, `datadir_temp()` and `tempdir`
* Function `get_data_dir()` is now deprecated in favour of `datadir()`

Documentation:

* Added question about file structure to [vignettes/FAQ.Rmd](https://github.com/spang-lab/metabodecon/blob/main/vignettes/FAQ.Rmd)
* Created categories for function reference in [_pkgdown.yml](https://github.com/spang-lab/metabodecon/blob/main/_pkgdown.yml)

Datasets:

* Moved `misc/datasets` to `misc/example_datasets`
* Moved `misc/examples/usage_example.R` to `misc/code_examples/sage_example.R`

Internal:

* Added unit tests
* Removed script `check_package.R`
* Moved functions from `util.R` to `datadir.R`
* Added `grDevices`, `stats` and `utils` as internal imports
* Added lots of test helper functions in `R/test_helpers.R`
* Added function `generate_lorentz_curves_v2()` which will replace `generate_lorentz_curves()` as soon as we have new features AND 100% backwards compatibility
* Fixed bug in `with()` that caused `get_datadir_mock()` to be called after redirection took place causing unexpected message output
* Fixed bug in `datadir()` that caused the resulting path to end with a slash on Unix-like systems and without a slash on Windows, if `file` was not specified
* `RUN_SLOW_TESTS` is now set to TRUE for the CI pipeline

# metabodecon 1.0.3

API:

* Updated `get_data_dir()` to accept `"blood"` as new value for parameter `dataset_name`
* Updated `download_example_datasets()` to download the datasets from the github repo instead of the old spang-lab repo

Documentation:

* Removed table of contents from `README.md` as it's a bit overkill for approx. 50 lines of text
* Improved documentation

Internal:

* Switched from MIT License to GPL-3 to match the license of the predecessor package `MetaboDecon1D`
* Added `docs` folder to `.gitignore`. Reason: we changed all vignettes to pkgdown articles which will be displayed only at our Github Pages website and can be regenerated from folder `vignettes` upon deployment.
* Created `TODOS.md` and added it to `.Rbuildignore`
* Improved `.gitignore`

# metabodecon 1.0.2

* Minor URL and spelling adjustments to pass CRAN checks

# metabodecon 1.0.1

* Fixed some spelling errors.
* Removed unused `CONTRIBUTE.md` (instead a section within `README.md` is used)

# metabodecon 1.0.0

* Initial CRAN submission.
