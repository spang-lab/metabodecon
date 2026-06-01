# Instructions for MetaboDecon

## IMPORTANT

- Never change any files outside the project or the global temp
  directory without asking for permission first.
- Never run `git commit` without an explicit instruction from
  the user to do so. After implementing changes, leave them in
  the working tree (or staged at most) so the user can review
  and commit themselves. The same rule applies to `git push`,
  `git tag`, `git reset --hard` and any other history-touching
  operation: explicit user instruction is required.

## Versioning

- The v2.0.0 release of metabodecon has **not** yet shipped. The
  `tobi` branch prepares the v2.0.1 release.
- During this preparation phase the version in `DESCRIPTION` may be
  bumped in increments of `2.0.0.x` (e.g. `2.0.0.1`, `2.0.0.2`, ...)
  for intermediate working revisions. Bump the version whenever a
  user-visible change ships, and add a one-line entry to `NEWS.md`.
- Because v2 has not yet been released, **no backwards-compatibility
  shims, deprecation warnings, or migration paths are required** while
  working on any `2.0.0.x` revision or the `tobi` branch. Break
  signatures, defaults, and behavior as needed and just update the
  call sites.

## Workspace Conventions

- Use `./tmp` for temporary files (it is gitignored and Rbuildignored).
  Never write temp files outside the project directory.
- If code changes must be tested with multiple workers, remind the user to
  reinstall the package before they run those tests.
- **Installing the package**: always use
  `R CMD INSTALL --no-lock --no-staged-install .`. The user runs radian
  sessions that mmap the package's `.so` files, which makes plain
  `R CMD INSTALL` fail with stale `00LOCK-metabodecon` directories
  (NFS silly-renames the in-use shared library to a `.nfs*` file that
  can't be removed). `--no-lock` skips the lock-dir step and
  `--no-staged-install` writes directly to the final library path.
  Both flags together avoid the stale-lock failure mode without
  touching the user's running R sessions.
- Never run `scripts/kill-nfs-lock.sh` (or `kill`/`pkill` against R
  processes) just to clear a stale install lock. Doing so terminates
  the user's radian sessions. The `--no-lock --no-staged-install`
  approach above sidesteps the lock entirely.

## Coding Guidelines

- All roxygen comments should start with a tag, in particular title and
  description should be formatted as `#' @title ...` and `#' @description ...`.
- Soft Character limit is 80.
  Hard limit is 100.
  Prefer short variable names like `x` and `y` to achieve that.
	If it's not clear from the function docs what a variable means,
	use a comment to describe it upon first use.
- Fill lines up to ~80 chars to minimize vertical space.
  Prefer single-line calls over multi-line when they fit.
	Going slightly over 80 (up to ~100) is tolerable if it avoids splitting
	a call across multiple lines lines.
- Prefer the use of helper variables instead of function nesting to reduce line
  length and improve readability. E.g. `x <- f(a); y <- g(x)` instead of `y <-
  g(f(a))`. Function nesting is ok if everything still fits in 80 chars and the
  function names are short and readable.
- Avoid defining functions inside other functions, except for temporary
  callbacks passed directly to `lapply()`, `sapply()`, `vapply()`, etc.
  Prefer top-level private helper functions instead to keep function bodies
  short and easy to scan.
- Do not use pipe operators `%>%` or `|>`.
- Always use fully qualified names for functions from other packages, e.g.
  `ggplot2::ggplot()`. Exceptions are functions from R's standard library like
  `sum()`, `mean()`, etc.
- Do NOT use spaces around `=` when passing arguments to functions.
  Good: `foo(bar=2)`. Bad: `foo(bar = 2)`.
- Always use fully qualified function names inside roxygen2 docs,
  even when referring to package internal functions (i.e., write
  `[metabodecon::deconvolute()]` instead of just
  `[deconvolute()]`)

### Multi-line formatting (function calls, `if/else`, loops): READ THIS

This rule is **non-negotiable**. The user has flagged violations of it
many times. Re-read this section before writing any R code.

**Rule 1 — Prefer single lines.** If a call, `if`/`else` branch, or
loop body fits within the project's character limit on one line, put
it on one line. Only split if it doesn't fit.

**Rule 2 — When you DO split, ALWAYS use this exact shape:**

```r
funcName(
    arg1, arg2,
    arg3, arg4
)
```

- Opening `(` is the **last character** on the function-name line.
  **Nothing else** follows it — no args, no comments, nothing.
- Args sit on subsequent lines, **indented exactly 4 spaces** beyond
  the start column of the call.
- Closing `)` is **alone on its own line**, dedented to the call's
  start column.

**Rule 3 — Never use "trailing args after the open paren".** This
style is **strictly forbidden**:

```r
# FORBIDDEN — args trailing after open paren, aligned to opening column:
funcName(arg1, arg2,
         arg3, arg4)

# FORBIDDEN — first arg on same line as funcName, closing ) inline:
funcName(arg1,
    arg2, arg3)

# FORBIDDEN — even if only one arg trails:
img_cached("path/to/file.rds",
           expensive_call(a, b, c))
```

If you find yourself aligning arguments under the open paren of the
function name, **stop and reformat using Rule 2**.

**Rule 4 — `if`/`else` follows the same shape.** Don't tail an `else`
branch onto the closing `}` of the `if` body — give it its own block:

```r
# Good:
v <- if (cond) {
    do_a()
    do_b()
} else {
    do_c()
}

# Forbidden — naked else expression after }:
v <- if (cond) {
    do_a()
    do_b()
} else do_c()
```

**Rule 5 — These rules apply recursively to nested calls.** If an
inner call has to wrap, the inner call itself follows Rule 2 — open
paren at end of line, args indented +4, closing paren on its own line.

This rule overrides any aesthetic preference for vertical alignment.
The user does not want column-aligned multi-line calls in this project.
Ever.

## Pseudocode mode for vignettes

- Definition: "pseudocode mode" means code examples are written primarily for
  readability and teaching, while still running under normal/ideal conditions.
- Scope: Use pseudocode mode for all newly added or modified code in
  `vignettes/*.Rmd`.
- Keep code compact and direct: prefer fewer lines, fewer helper variables, and
  simple control flow.
- Prefer clear intent over defensive robustness in vignettes. Avoid extra
  safeguards that distract from the main idea unless they are essential to
  understand the method.
- Keep naming short when context is obvious (e.g. `X`, `y`, `te`, `Xtr`,
  `Xte`).
- Prefer single-line function calls in vignettes. Avoid multiline calls unless
  required for readability of `*apply()` loops or unavoidable long literals.
- Avoid defensive checks and fallback branches in vignette chunks. In pseudocode
  mode, prefer short, direct, didactic code that assumes normal conditions.

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

Functions for aligning deconvoluted spectra. The public alignment pipeline
chains two stages — **CluPA** (continuous shifts) and **RefPA** (discrete
snap to reference) — and [metabodecon::align()] runs both in one call.

- (exported) `align(x, maxShift, maxCombine, ref=NULL, ...)`: chains
  `clupa()` then `snap_to_ref()`. Returns an `aligns` object whose
  per-spectrum `lcpar` has been collapsed onto the reference's peak grid.
- (exported) `clupa(x, maxShift, ref=NULL, ...)`: **CluPA** —
  hierarchical-clustering FFT segment shifts (Beirnaert et al. 2018,
  Vu et al. 2011). FFT input is the Lorentz superposition `$sit$sup`
  already attached at deconvolution time (= speaq-equivalent shape,
  matches v1.7.0's `get_sup_mat(decons2)` → `dohCluster` input).
  Requires every spectrum in `x` to share the same `$cs` grid — call
  `harmonize_grid(x)` upstream if your corpus is from different
  acquisitions. Adds `x0al`, `pcial` to `lcpar`; keeps original peak
  count.
- (exported) `snap_to_ref(x, maxCombine, ref=NULL)`: **RefPA** — for
  every peak, records the nearest reference column on the shared `cs`
  grid (within `maxCombine`) as `pcisn` / `x0sn`. Peaks farther than
  `maxCombine` get `pcisn = NA` / `x0sn = NA`. Original `x0`, `x0al`,
  `A`, `lambda`, `pcide`, `pcial` are all preserved — `snap_to_ref`
  only *adds* fields. No peaks are dropped here and amplitudes are
  not summed; [metabodecon::si_mat()] / [metabodecon::peak_mat()]
  skip `pcisn = NA` peaks and sum collisions on the same `pcisn`
  column at rasterisation time. Clears `sit$supal` (the post-snap
  peak list is no longer Lorentz-compatible).
- (exported) `identity_align(x, ...)`: no-op; returns `x`.
- (private) `ensure_shared_cs`: assertion-only helper used by every
  alignment / snap entry point. Stops with an actionable message if
  inputs don't share a grid; remedy is to call `harmonize_grid()`.
- (private) `ensure_align_aux`: per-spectrum kernel that backfills
  `lcpar$pcide` and `sit$sup` from `x$cs` if missing.
- (private) `noshift_align`, `noshift_one`: CluPA's `maxShift = 0`
  fast-path (sets `x0al = x0`, `pcial = nearest cs column`).
- (private) `align_decon`: CluPA per-spectrum kernel (FFT shift +
  speaq-equivalent hclust). Reads `x$cs` and `x$sit$sup`; writes
  `x0al = cs[pcial]` and `pcial` as integer indices into `x$cs`.
- (private) `snap_lcpar`: RefPA per-spectrum kernel.
- (private) `find_ref`, `find_ref_ind`: pick the reference spectrum by
  minimising the sum, over every target peak in every other spectrum,
  of the ppm distance to the nearest peak in the candidate reference.
  Grid-free (compares `x0` values in ppm directly). **Bias:** the sum
  is over target peaks only, so candidates with dense peak lists
  (incl. noise peaks) are favoured. Mirrors `speaq::findRef` semantics.
- (private) `pci_on_cs`: integer column index for a vector of ppm
  values via `round(convert_pos(...))`, clamped to `[1, length(cs)]`.
- (private) `fft_shift`, `do_shift`, `hclust_align`, `pad_peaks`:
  bundled CluPA implementation that mirrors `speaq::hClustAlign` and is
  byte-equivalent to it; the speaq backend remains available via
  `use_speaq = TRUE`.

VOPA and GloPA were removed in 2.0.0 — the only built-in CluPA-stage
backends are `clupa` and `identity_align`.

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
- (private) `smooth_signals2`: Smooths signal intensities using a moving average.
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
4. `test-convert_sfr.R`
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
19. `test-get_sfr.R`
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
39. `test-smooth_signals2.R`
40. `test-snap_to_ref.R`
41. `test-speaq_align.R`
42. `test-vcomp.R`
