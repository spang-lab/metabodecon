# metabodecon 1.0.4

Currently under development in branch `batch`.

* Added `urine.dx` test file to `inst` folder to be able to write easy examples without requiring the user to download the data first. Solves todo `Feature: add minimal example dataset`.
* Added parameter `ask` to function `generate_lorentz_curves()` and parameters for all values currently obtained interactively from the user, i.e. `file_name`, `number_iterations`, `range_water_signal_ppm`, `signal_free_region`, `smoothing_param`, `delta` and, `scale_factor`

Only relevant for developers:

* Bumped version to 1.0.4
* Added the following todo to `TODOS.md`: *generate_lorentz_curves should not write to input folders by default*

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
