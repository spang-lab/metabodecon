# metabodecon 1.1.0

Currently under development in branch `improve_tests`.

* Improved function `download_example_datasets()`
* Replaced function `get_data_dir()` with `datadir()` and its helper functions `datadir_persistent()`, `datadir_temp()` and `tempdir`
* Added question about file structure to [vignettes/FAQ.Rmd](vignettes/FAQ.Rmd)
* Created categories for function reference in [_pkgdown.yml](_pkgdown.yml)
* Moved `misc/datasets` to `misc/example_datasets`
* Moved `misc/examples/usage_example.R` to `misc/code_examples/sage_example.R`
* Function `get_data_dir()` is now deprecated in favour of `datadir()`

Only relevant for developers:

* Added unit tests
* Removed script `check_package.R`
* Moved functions from `util.R` to `datadir.R`
* Added `grDevices`, `stats` and `utils` as internal imports
* Added lots of test helper functions in `R/test_helpers.R`
* Added issue `6. Fix: generate_lorentz_curves should not write to input folders by default` to `TODOS.md`
* Added function `generate_lorentz_curves_v2()` which will replace `generate_lorentz_curves()` as soon as we have new features AND 100% backwards compatibility
* Bumped version to `1.1.0` in `DESCRIPTION`
* Updated `NEWS.md`
* Fixed bug in `with()` that caused `get_datadir_mock()` to be called after redirection took place causing unexpected message output
* Fixed bug in `datadir()` that caused the resulting path to end with a slash on Unix-like systems and without a slash on Windows, if `file` was not specified
* `RUN_SLOW_TESTS` is now set to TRUE for the CI pipeline

# metabodecon 1.0.3

* Improved documentation
* Switched from MIT License to GPL-3 to match the license of the predecessor package `MetaboDecon1D`
* Updated `get_data_dir()` to accept `"blood"` as new value for parameter `dataset_name`
* Updated `download_example_datasets()` to download the datasets from the github repo instead of the old spang-lab repo

Only relevant for developers:

* Removed table of contents from `README.md` as it's a bit overkill for approx. 50 lines of text
* Added `docs` folder to `.gitignore`. Reason: we changed all vignettes to pkgdown articles which will be displayed only at our Github Pages website and can be regenerated from folder `vignettes` upon deployment.
* Added lots of todos to `TODOS.md`
* Added `TODOS.md` to `.Rbuildignore`
* Improved `.gitignore`

# metabodecon 1.0.2

* Minor URL and spelling adjustments to pass CRAN checks

# metabodecon 1.0.1

* Fixed some spelling errors.
* Removed unused `CONTRIBUTE.md` (instead a section within `README.md` is used)

# metabodecon 1.0.0

* Initial CRAN submission.
