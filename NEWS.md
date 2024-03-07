# metabodecon 1.2.0

Development Branch:

* `batch-mode`

API:

* Changed `download_example_datasets`. If `extract=TRUE`, the path to the extracted folder is now returned, instead of the path to the zip file.

Documentation:

* Added `acqus` and `1r` file to section *File Structure* of [FAQ](vignettes/FAQ.Rmd)
* Added info about bruker namings to [FAQ](vignettes/FAQ.Rmd)

Internal:

* Added helper functions `load_jcampdx_spectrum()` and `load_bruker_spectrum()` and corresponding test cases. These functions are more of less directly extract from `deconvolute_spectrum()`.
* Added improved versions of `load_jcampdx_spectrum()` and `load_bruker_spectrum()` named `load_jcampdx_spectrum_v2()` and `load_bruker_spectrum_v2()` incl. test cases
* Added util function `normPath()`.
* Added test cases for `deconvolute_spectrum()`
* Moved helper function `.deconvolute_spectrum()` into seperate file `deconvolute_spectrum_v2` for development and renamed it to `deconvolute_spectrum_v2()`.
* Added helper function `get_signal_free_region_in_su()`
* Disabled `object_length_linter`
* Moved `deconvolution` function `MetaboDecon1D.R` into seperate file `MetaboDecon1D_deconvolution.R`.
* Moved `str2`, `get_num_input`, `get_str_input`, `get_yn_input`, `collapse` from `main_v2.R` into `util.R`
* Added functions `ppm_to_dp()`, `ppm_to_sdp()`, `dp_to_ppm()` to `util.R`
* Added `GLOSSARY.md`
* Extracted code parts from `deconvolute_spectrum_v2()` and moved them into functions `determine_water_signal`, `determine_signal_free_region`, `plot_sfr` and `plot_ws`.
* Added test cases for `determine_water_signal` and `determine_signal_free_region`
* Moved function `deconvolution()` from `MetaboDecon1D.R` into `MetaboDecon1D_deconvolution.R`

# metabodecon 1.1.1

API:

* Fixed a bug in `generate_lorentz_curves()` that caused the function to always use file format "bruker", even when file format "jcampdx" was specified.

Datasets:

* Fixed filenames of samples in blood dataset (renamed from `Bood_<nr>` to `blood_<nr>`).
* Renamed `example_datasets/jcampdx/urine/urine.dx` to `example_datasets/jcampdx/urine/urine_1.dx` and renamed `example_datasets/bruker/urine/urine/` to `example_datasets/bruker/urine/urine_1/`. This was done because `list.files` seems to return different orderings for `urine.dx` and `urine_2.dx` in different operating systems, whereas `urine_1.dx` and `urine_2.dx` are sorted the same way everywhere. This makes it easier to write clear and concise test cases, because we don't need to check for file ordering.

Documentation:

* Fixed broken image in [vignettes/FAQ.Rmd](vignettes/FAQ.Rmd).

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

* Added question about file structure to [vignettes/FAQ.Rmd](vignettes/FAQ.Rmd)
* Created categories for function reference in [_pkgdown.yml](_pkgdown.yml)

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
