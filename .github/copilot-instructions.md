# Copilot Instructions for MetaboDecon

## Coding Guidelines

- All roxygen comments should start with a tag, in particular title and
  description should be formatted as `#' @title ...` and `#' @description ...`.
- Character limit is 80. Prefer short variable names like `x` and `y` to achieve
  that. If it's not clear from the function docs what a variable means, use a
  comment to describe it upon first use.
- Prefer the use of helper variables instead of function nesting to reduce line
  length and improve readability. E.g. `x <- f(a); y <- g(x)` instead of `y <-
  g(f(a))`. Function nesting is ok if everything still fits in 80 chars and the
  function names are short and readable.
- Do not use pipe operators `%>%` or `|>`.
- Always use fully qualified names for functions from other packages, e.g.
  `ggplot2::ggplot()`. Exceptions are functions from R's standard library like
  `sum()`, `mean()`, etc.

## Project Structure

- `_pkgdown.yml`: Configuration file for the pkgdown website.
- `ARCHIVE.md`: Archive of old project notes or documentation.
- `cran-comments.md`: Comments for CRAN submission.
- `CRAN-SUBMISSION`: Details about the CRAN submission process.
- `DESCRIPTION`: Metadata about the R package.
- `Dockerfile`: Instructions to build a Docker image for the project.
- `LICENSE.md`: Licensing information for the project.
- `NAMESPACE`: Defines the exported and imported functions for the package.
- `NEWS.md`: Changelog for the project.
- `README.md`: Overview and instructions for the project.
- `TODOS.md`: List of tasks and future improvements.
- `data/`: Contains example datasets (e.g., `sap.rda`, `sim.rda`).
- `docs/`: Generated documentation for the package.
- `inst/`: Additional files to be included in the package.
- `man/`: Documentation files for R functions.
- `misc/`: Miscellaneous files, including code examples and sketches.
- `pkgdown/`: Assets for the pkgdown website.
- `R/`: Contains the R scripts for the package.
- `tests/`: Unit tests for the package.
- `vignettes/`: Long-form documentation and tutorials.

## Vignettes

- `Classes.Rmd`: Documentation for the `Classes` module.
- `Compatibility.Rmd`: Details on compatibility with other tools.
- `Contributing.Rmd`: Guidelines for contributing to the project.
- `Datasets.Rmd`: Information about the datasets included in the package.
- `FAQ.Rmd`: Frequently asked questions.
- `Get_Started.Rmd`: A guide to getting started with the package.
- `MetaboDecon1D.Rmd`: Documentation for the `MetaboDecon1D` module.

## Modules

### align.R

Functions for aligning deconvoluted spectra.

- (exported) `align`: Aligns spectra data using the 'CluPA' algorithm from the 'speaq' package.
- (internal) `get_ppm_range`: Returns the ppm range covered by spectra.
- (internal) `gen_feat_mat`: Generates a feature matrix.
- (internal) `speaq_align`: Aligns signals using the 'speaq' package.
- (internal) `combine_peaks`: Combines peaks across spectra.
- (internal) `dohCluster`: Performs cluster-based peak alignment.
- (private) `get_decon_params`: Extracts deconvolution parameters.
- (private) `read_decon_params`: Reads deconvolution parameters from files.
- (private) `check_decon_params`: Checks the validity of deconvolution parameters.
- (private) `rm_zero_width_peaks`: Removes peaks with zero width.
- (private) `is_decon_obj`: Checks if an object is a deconvolution object.
- (private) `is_decon_list`: Checks if a list contains deconvolution objects.
- (private) `get_peak_indices`: Retrieves peak indices from a deconvolution object.
- (private) `get_sup_mat`: Constructs a signal intensity matrix.

Internal means, the function is currently exported, but should no longer be used directly and will be removed or made private in the future.

### class.R

Class definitions and methods for the classes used by the package.

- (exported) `print.${PUBLIC_CLASS}`: Prints a spectrum object.
- (exported) `is_${PUBLIC_CLASS}`: Checks if an object is a spectrum.
- (exported) `as_${PUBLIC_CLASS}`: Converts an object to a spectrum.
- (private) `is_{PRIVATE_CLASS}`: Checks if an object is an ispec.
- (private) `is_spectrum_or_spectra`: Checks if an object is a spectrum or spectra.
- (private) `get_name`: Retrieves the name of a MetaboDecon object.
- (private) `get_names`: Retrieves the names of a collection of MetaboDecon objects.
- (private) `set_names`: Sets the names of a collection of MetaboDecon objects.

With:

- PUBLIC_CLASS in: spectrum, decon1, decon2, align, spectra, decons1, decons2, aligns
- PRIVATE_CLASS in: ispec, idecon, rdecon

### data.R

Function for creating, updating and downloading example datasets.

- (exported) `download_example_datasets`: Downloads example datasets for testing.
- (exported) `metabodecon_file`: Returns the path to a file or directory in the package.
- (exported) `datadir`: Returns the path to the data directory.
- (exported) `datadir_persistent`: Returns the path to the persistent data directory.
- (exported) `datadir_temp`: Returns the path to the temporary data directory.
- (exported) `tmpdir`: Returns the path to the temporary session directory.
- (exported) `get_data_dir`: Deprecated function to retrieve the directory path of an example dataset.
- (prviate)`cache_example_datasets`: Caches example datasets.
- (prviate)`extract_example_datasets`: Extracts example datasets from a zip file.
- (prviate)`download_example_datasets_zip`: Downloads the example datasets zip file.
- (prviate)`zip_temp`: Returns the path to the temporary zip file.
- (prviate)`zip_persistent`: Returns the path to the persistent zip file.
- (prviate)`tmpfile`: Creates a temporary file.
- (prviate)`testdir`: Returns the path to a test directory.
- (prviate)`mockdir`: Returns the path to a mock directory.
- (prviate)`cachedir`: Creates and returns a cache directory.
- (prviate)`make_sap`: Creates the SAP dataset.
- (prviate)`update_sap`: Updates the SAP dataset.
- (prviate)`make_sim`: Creates the Sim dataset.
- (prviate)`update_sim`: Updates the Sim dataset.
- (prviate)`deconvolute_blood`: Deconvolutes the Blood dataset.
- (prviate)`get_sim_params`: Retrieves simulation parameters from a deconvolution object.

### decon.R

Functions for deconvoluting NMR spectra.

- (exported) `deconvolute`: Deconvolutes NMR spectra by modeling signals as Lorentz curves.
- (exported) `generate_lorentz_curves`: Generates Lorentz curves for spectra.
- (exported) `generate_lorentz_curves_sim`: Optimized for the "Sim" dataset.
- (private) `deconvolute_spectra`: Internal function for deconvoluting multiple spectra.
- (private) `deconvolute_spectrum`: Internal function for deconvoluting a single spectrum.
- (private) `smooth_signals`: Smooths signal intensities using a moving average.
- (private) `find_peaks`: Detects peaks in the spectrum.
- (private) `filter_peaks`: Filters peaks with low scores outside the signal-free region.
- (private) `fit_lorentz_curves`: Fits Lorentz curves to the detected peaks.

### depr.R

Deprecated functions and classes.

- (exported) `MetaboDecon1D`: Deprecated function for deconvoluting 1D NMR spectra.
- (exported) `calculate_lorentz_curves`: Calculates Lorentz curves for analyzed spectra.
- (exported) `plot_triplets`: Plots peak triplets (deprecated).
- (exported) `plot_lorentz_curves_save_as_png`: Plots Lorentz curves and saves as PNG (deprecated).
- (exported) `plot_spectrum_superposition_save_as_png`: Plots spectrum superposition and saves as PNG (deprecated).
- (private) `deconvolution`
- (private) `plot_si_mat`
- (private) `plot_sim_spec`
- (private) `plot_noise_methods`
- (private) `read_decon_params_original`
- (private) `speaq_align_original`
- (private) `simulate_from_decon`
- (private) `count_stretches`
- (private) `analyze_noise_methods`

### mdrb.R

Functions for installing and checking mdrb (metabodecon rust backend).

- (exported) `check_mdrb`: Checks if the Rust backend is installed.
- (exported) `check_mdrb_deps`: Checks dependencies for the Rust backend.
- (exported) `install_mdrb`: Installs the Rust backend.
- (prviate) `get_mdrb_version`: Retrieves the version of the Rust backend.

### paper.R

Functions for creating the figures for the metabodecon 2025 paper.

- (private) `mkfig_nmr_challenges`: Creates a figure illustrating typical NMR challenges.
- (private) `test_plot_nmr_challenges`: Tests the plotting of NMR challenges interactively or non-interactively.
- (private) `plot_nmr_challenges`: Plots a series of subplots for NMR challenges.
- (private) `plot_1_nmr_experiment` to `plot_10_annotated_spectra`: Helper functions for individual subplots.
- (private) `draw_vial` and `draw_nmr_spectrometer`: Draws specific components like vials and spectrometers.
- Various helpers: Functions like `init_dev`, `fill_dev`, and `marbox` assist in plotting.

### plot.R

Functions for plotting single and multiple spectra before and after deconvolution/alignment.

- (exported) `plot_spectra`: Plots a set of deconvoluted spectra.
- (exported) `plot_spectrum`: Plots a single spectrum with zoomed regions.
- (exported) `draw_spectrum`: Draws a single spectrum, used internally by `plot_spectrum`.
- (private) `plot_sfr`: Plots the signal-free region.
- (private) `plot_ws`: Plots the water signal region.
- (private) `plot_align`: Plots aligned and unaligned spectra for comparison.
- (private) `plot_empty`: Creates an empty plot canvas.
- (private) `plot_dummy`: Creates a dummy plot for testing.
- Various helpers: Functions like `draw_legend`, `draw_con_lines`, and `draw_lc_line` assist in drawing specific plot elements.

### spectrum.R

Functions for creating, reading and writing spectra from/to disk.

- (exported) `read_spectrum`: Reads a single spectrum from disk.
- (exported) `read_spectra`: Reads multiple spectra from disk.
- (exported) `make_spectrum`: Creates a spectrum object.
- (exported) `simulate_spectrum`: Simulates a 1D NMR spectrum.
- (private) `read_bruker`
- (private) `read_jcampdx`
- (private) `parse_metadata`
- (private) `read_acqus`
- (private) `read_procs_file`
- (private) `read_simpar`
- (private) `read_one_r`
- (private) `save_spectrum`
- (private) `save_spectra`

### test.R

Utility functions to help with testing

Exported Functions:

- (exported) `evalwith`: Evaluates an expression with predefined global state, including options for capturing output, mocking directories, and caching results.
- (unexported) `get_readline_mock`: Creates a mock `readline` function for testing.
- (unexported) `get_datadir_mock`: Returns a mock for the `datadir` functions.
- (unexported) `run_tests`: Runs tests with options to skip slow tests or focus on specific functions.
- (unexported) `skip_if_slow_tests_disabled`: Skips tests if slow tests are disabled.
- (unexported) `skip_if_not_in_globenv`: Skips tests if not in the global environment.
- (unexported) `expect_file_size`: Checks if file sizes in a directory are within a certain range.
- (unexported) `expect_str`: Tests if the structure of an object matches an expected string.
- (unexported) `vcomp`: Compares two vectors and prints differences.
- (unexported) `compare_spectra`: Compares spectra deconvoluted with different methods.
- (unexported) `calc_prarp`, `calc_prarpx`, `plot_prarp`: Functions for calculating and visualizing PRARP scores.
- (unexported) `MetaboDecon1D_silent`, `MetaboDecon1D_silent_sim`: Silent wrappers for `MetaboDecon1D`.
- (unexported) `get_MetaboDecon1D_answers`: Generates answers for `MetaboDecon1D`.

### util.R

Utility functions for unit conversion, file handling, user input, type checking, and other miscellaneous tasks. It includes both exported and unexported functions to support various operations within the package.

- (Exported) Unit Conversion:
  - `convert_pos`: Converts positions from one unit to another.
  - `convert_width`: Converts widths from one unit to another.
  - `width`: Calculates the width of a numeric vector.

- (Exported) File Handling:
  - `checksum`: Calculates a checksum for files or directories.
  - `tree`: Prints the structure of a directory tree.

- (Exported) Visualization:
  - `transp`: Makes a color transparent by adding an alpha channel.

- (Private) Unit Conversion:
  - `in_hz`: Converts chemical shifts from ppm to Hz.
  - `sfr_in_ppm_bwc`: Converts signal-free region borders from SDP to PPM.
  - `sfr_in_sdp_bwc`: Converts signal-free region borders from PPM to SDP.

- (Private) File Handling:
  - `mkdirs`: Recursively creates directories.
  - `clear`: Clears a directory and recreates it.
  - `norm_path`: Normalizes file paths.
  - `pkg_file`: Returns the path to a file within the package.
  - `store`: Stores an object in a file.

- (Private) User Input:
  - `readline`: Mockable version of `readline`.
  - `get_num_input`: Prompts the user for numeric input.
  - `get_int_input`: Prompts the user for integer input.
  - `get_str_input`: Prompts the user for string input.
  - `get_yn_input`: Prompts the user for yes/no input.

- (Private) Operators:
  - `%||%`, `%&&%`, `%==%`, `%!=%`, `%===%`, `%!==%`, `%notin`: Custom operators for logical and equality checks.

- (Private) Type Checking:
  - Functions like `is_num`, `is_int`, `is_char`, `is_bool`, and their variations check the type and structure of objects.

- (Private) Miscellaneous:
  - `logf`, `stopf`, `human_readable`: Logging and formatting utilities.
  - `du`: Prints the size of an object and its subcomponents.
  - `set`, `pop`: Modifies lists in-place.
  - `timestamp`: Returns the current timestamp.
  - `mcmapply`: Multi-core version of `mapply`.

### zzz.R

Package stuff like .onLoad, imports and a `get_started` functions.

## Tests

Test live inside ./tests/testthat. The followings tests exist:

1. `test-align.R`
2. `test-as_decon.R`
3. `test-cache_example_datasets.R`
4. `test-combine_scores.R`
5. `test-convert_sfr.R`
6. `test-convert_spectrum.R`
7. `test-convert_wsr.R`
8. `test-datadir.R`
9.  `test-deconvolute_spectra.R`
10. `test-deconvolute_spectrum.R`
11. `test-deconvolute.R`
12. `test-download_example_datasets.R`
13. `test-draw_spectrum.R`
14. `test-evalwith.R`
15. `test-generate_lorentz_curves.R`
16. `test-get_decon_params.R`
17. `test-get_names.R`
18. `test-get_sfr.R`
19. `test-get_smopts.R`
20. `test-get_wshw.R`
21. `test-init_lorentz_curves.R`
22. `test-is_float_str.R`
23. `test-is_int_str.R`
24. `test-lorentz_int.R`
25. `test-lorentz_sup.R`
26. `test-mcmapply.R`
27. `test-metabodecon_file.R`
28. `test-MetaboDecon1d.R`
29. `test-pkg_file.R`
30. `test-plot_sfr.R`
31. `test-plot_spectrum.R`
32. `test-plot_ws.R`
33. `test-read_acqus_file.R`
34. `test-read_bruker.R`
35. `test-read_decon_params.R`
36. `test-read_one_r_file.R`
37. `test-read_procs_file.R`
38. `test-read_spectrum.R`
39. `test-smooth_signals.R`
40. `test-speaq_align.R`
41. `test-vcomp.R`
