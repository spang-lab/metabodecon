# DONE WITH 1.4.2 (PR 14)

## [x] Remove Remotes field from DESCRIPTION

When you submit your package to CRAN, all of its dependencies must also be available on CRAN. For this reason, release() will warn you if you try to release a package with a Remotes field.

## [x] Check missing pkgs in align

Installation via `install.packages("metabodecon")` does not install `MassSpecWavelet` and `impute`. So if a user doesn't copy paste the installation instructions but installed via `install.packages("metabodecon")`, these dependencies will be missing. In such scenarios, we should print an error message with the required install commands and abort.

## [x] Change verbose defaults

Change verbose argument for deconvolute and align to TRUE.

## [x] Shrink test-install workflow

Currently the test-install workflow is split over three jobs, with huge amounts of code copy pasted. Extract the R commands into a single `test-install.R` R script, that can be called as `Rscript -e test-install.R <method>` and,

1. Deletes all previously available dependencies incl. Rtools on Windows (i.e. it must be removed from the PATH)
2. Does the installation according to the `method` commandline argument, e.g. "CRAN-Modern", "CRAN-Old" or "Github".

The corresponding commands for "CRAN-Modern", "CRAN-Old" or "Github" can be take from the current version `test-install.yaml`.

After ou created `test-install.R`, update the workflow to use it.

# PLANNED

## Improve questions

1. Question `Signal free region correctly selected? (y/n)` should be replaced by `Borders to Signal Free Regions (green) correctly selected? (y/n)`
2. Question `Water artefact fully inside red vertical lines? (y/n)` should be replaced by `Water artefact fully inside blue area? (y/n)`

## Update expno/procno defaults

Expno and procno should use NULL as default.
If NULL, the function should look for the first available expno/procno folder.
If there is only one, it should use that one.
If there are multiple, it should use expno=10 and procno=10.
If there are multiple but no 10, it should throw an error.

## Add function get_si_mat

Add a function `get_si_mat()` for extracting a matrix of signal intensities (SI) from a metabodecon object. The type of returned SI should be `raw` for `spectra`, `sup` for `decons` and `al` for `aligns`.

## Add authors to functions

## Add lifecycle badges to functions

## Improve Sap Dataset

Sap spectra should be simulated as follows:

1. Find three metabolites `mets` related to diabetes with only 1-3 peaks each.
2. Look up their signal centers `x0_` and halfwidths `lambda_`.
3. Define two groups: `healthy` and `diabetes`.
4. Draw signal areas for each sample from a distribution. The distribution should be different for `healthy` and `diabetes`.
5. Draw metabolite shifts for each sample from a normal distribution.
6. Apply the shifts to the signal centers.
7. Simulate data as usual using `simulate_spectra`
8. Make sure, all of the above information is stored inside `$meta$simpar`

## Improve Get-Started vignette

Add nicer plots.
Improve alignment part.

## Write Deconvolution-Details Vignette

## Write Alignment Details Vignette

## Analyze runtime

Do a benchmark about which parts take up the most time in `deconvolute` and `align` and create issues for improving them. Some of the slow parts should already be mentioned in issues. If this is the case, increase the priority of these issue.

## Write paper

Reformat the vignettes as paper and send to Wolfram for proofreading.

## Warn user if peaks in SFR

If delta is small (e.g. 1), peaks in SFR might not be filtered out. Either implement this and warn user about it (this is a strong indication that delta was chosen too small).

Note 23.1.2025: in Rust implementation peaks found in IGNORE-REGION are already filtered out automatically before parameter approximation AND points in IGNORE-REGIONS do not contribute to MSE and/or PRARP.

## Improve mse_normed calculation

In function `add_return_list`:

1. Make the following part faster (or remove completely):

   `s <- sapply(x, function(x_i) sum(abs(A * (lambda / (lambda^2 + (x_i - w)^2))))) # takes approx. 2.2 seconds for urine_1`

2. Return `mse_normed_raw` in addition to `mse_normed` (which is calculated based on `y <- spec$y_smooth`). `mse_normed_raw` should be based on `y <- spec$y_raw`.

## Check negative A values

Check why there are negative values for the estimated lorentz curve area A.

## Improve SFR and WS defaults

Replace the default values `wshw = 0.1527692` and `sfr = c(11.44494, -1.8828)` in `generate_lorentz_curves()` with `wshw = "auto"` and `sfr = "auto"`, which should be calculated as follows:

If `c(11.44494, -1.8828)` is part of the ppm range, use these values, otherwise calculate them as

1. `wshw = 0.01 * width(cs)` (where `0.01` is `round(0.007629452, 2)` and `0.007629452` equals `0.1527692 / 20.0236144338963` which is the width of the default WSHW dividided by the width of the `urine_1` spectrum. I.e., the new calculation would give approximately the same proportion of the spectrum width as the default value.)
2. `sfr = max(cs) - c(1/6, 5/6) * width(cs)`

## Refactor integral calculations

We should use `A * pi` everywhere unless `bwc = 0`.

## Calculate mandatory SITs in as_decon2

Elements `wsrm` and `nvrm` are mandatory in `decon2` and should therefor be calculated during `as_decon2`. If additional data is needed to do this from `decon0` and/or `decon1` objects, add an argument that allows users to provide this information.

## Refactor mse calculations

We should use `mse()` everywhere.

## Refactor parameter approximation

Implement [Max' parameter approximation algorithms](https://gitlab.spang-lab.de/bachelorthesis/ws2425_msombke_metabodecon-v2/-/blob/main/new_method_docs/main.R
) in `calc_A`, `calc_lambda` and `calc_w`.

## Test deconvolute with huge nfit

## Test deconvolute with peak at zero ppm

It's unclear why <https://github.com/spang-lab/metabodecon/blob/v1.2.0/R/21_decon_helpers.R#L315> checks for w[i] == 0.

Assumption: it's a bug and leads to this peak being missed. If this is the case, remove the check from the code as well.

## Remove unneeded checks

Check special handling for cases with A[i] == 0 and lambda[i] == 0 in parameter approximaton. Max analyzed it and concluded that checks are not necessary. My thought: copy paste artifacts from the the w[i] == 0 check (which is wrong).

# BACKLOG

## Show prarp in plot_spectrum

## Show peak scores in plot_spectrum

## Check spectrum type

When reading spectra, the spectrum type should be checked (1D, 2D, etc. and if ncessary an error message should be printed)

## Integrate jcampdx functions

Integrate jcampdx reader and writer functions from Max.

## Remove ispec and idecon

Subtasks:

1. Modify below functions in a way that they return a single value instead
   returning the full input spectrum. The assignement to the spectrum object
   should be done by the caller.
2. Modify below functions in a way, that they don't accept idecon objects, but
   only spectrum objects and/or other v1.2+ elements
3. Remove the `idecon` type from the package.

Functions:

1.  `enrich_sfr(sfr, ispec)`
2.  `enrich_wshw(wshw, ispec)`
3.  `init_lc(ispec)`
4.  `refine_lc_v14(ispec, lc$Z)`
5.  `rm_water_signal(ispec, wshw, bwc)`
6.  `rm_negative_signals(ispec)`
7.  `smooth_signals(ispec, smopts[1], smopts[2], bwc)`
8.  `find_peaks(ispec)`
9.  `filter_peaks(ispec, sfr, delta, force, bwc)`
10. `fit_lorentz_curves(ispec, nfit, bwc)`

Links:

- https://spang-lab.github.io/metabodecon/articles/Classes.html#spectrum
- https://spang-lab.github.io/metabodecon/articles/Classes.html#elements-v1-2

## Turn articles into vignettes

Make sure that all suitable articles are included as vignettes and built from scratch so they don't take up too much space.

# DONE

## DONE WITH PR 12 (v1.4.0)

### Test install on clean OS

Add a workflow for testing installation on a clean Windows/Linux/Mac OS with R pre-installed, but without R-tools and any packages.

*Done: Mar 10. Branch: mdrb. PR: https://github.com/spang-lab/metabodecon/pull/12.*

### Add Rust Backend Installer

The Rust backend should only be built for R versions >= 4.2.0 and if RTools is available. If any of these conditions is not fulfilled, compilation should be skipped.

If skipping compilation is not possible, because we load the dynamic lib in NAMESPACE, it should be generated in a way that does not use any features required by R 4.2.0 or greater. If we choose backend == "rust", we need to check whether the required functions are available.

Update 10.3.2025: Instead of making compilation optional we should provide a seperate R package [mdrb](https://github.com/spang-lab/mdrb) (Metabodecon Rust Backend) which we can install upon request. I.e. we need to implement:

1. A function `check_mdrb()` that checks whether mdrb is already installed.
2. A function `check_mdrb_deps()` that checks whether mdrb can be installed. Required dependencies are
   - R version 4.2 or higher
   - Build tools for R (RTools on Windows, build-essentials on Linux, XCode on Mac)
   - Cargo and rustc version 1.78 or higher
3. A function `install_mdrb()`, that
   1. Does nothing if mdrb is already installed
   2. Prints installation instructions for mdrb dependencies if `check_mdrb_deps()` lists missing dependencies
   3. Calls `pak::pkg_install("spang-lab/mdrb")` if all requirements are satisfied

*Done: Mar 14. Branch: mdrb. PR: https://github.com/spang-lab/metabodecon/pull/12.*

### Rename deconvolute_ispec

This is a preparation for issue *Add Rust Backend Argument*. Rename `deconvolute_ispec()` to `deconvolute_spectrum()` and `deconvolute_ispecs()` to `deconvolute_spectra()`.

*Done: Mar 18-22. Branch: mdrb. PR: https://github.com/spang-lab/metabodecon/pull/12.*

### Add Rust Backend Argument

Add an experimental argument `use_rust` in `deconvolute()` causing the following behaviour:

| use_rust | check_mdrb | call                        |
| -------- | ---------- | --------------------------- |
| FALSE    | anything   | deconvolute_spectrum_r()    |
| NULL     | FALSE      | deconvolute_spectrum_r()    |
| NULL     | TRUE       | deconvolute_spectrum_rust() |
| TRUE     | TRUE       | deconvolute_spectrum_rust() |
| TRUE     | FALSE      | stop(MESSAGE)               |

With MESSAGE being something like "Rust backend not installed yet. Please call install_mdrb() first."

Sub-Tasks

- [x] Add `use_rust` argument to `deconvolute`, `deconvolute_spectra()` and `deconvolute_spectrum()` and assert that it is passed on correctly
- [x] Implement usage of rust in `deconvolute_spectrum()`
- [x] Write testcases for `deconvolute_spectrum(use_rust=TRUE, rtyp="idecon')`
- [x] Write testcases for `deconvolute_spectrum(use_rust=TRUE, rtyp="decon2')`
- [x] Write testcases for `deconvolute_spectra(use_rust=TRUE, rtyp="idecon')`
- [x] Write testcases for `deconvolute_spectra(use_rust=TRUE, rtyp="decon2')`
- [x] Write testcases for `deconvolute(use_rust=TRUE)` with mdrb installed
- [x] Write testcases for `deconvolute(use_rust=TRUE)` with mdrb missing

*Done: Mar 18-25. Branch: mdrb. PR: https://github.com/spang-lab/metabodecon/pull/12.*

## DONE WITH PR 13 (v1.4.1)

### Check mse calculation in Rust

`x$mdrb_decon$mse()` deviates from `mse(si, sup, norm=FALSE)` in `as_decon2.rdecon()`, which is why we have to calculate the MSE ourselves in R. This should be fixed in mdrb.

MSE calculation in R is implemented in `decon.R` as:

```R
mse <- function(y, yhat, normed = FALSE) {
    if (normed) {
        mean(((y / sum(y)) - (yhat / sum(yhat)))^2)
    } else {
        mean((y - yhat)^2)
    }
}
```

*Update 3. April 2025: clarified with Maximilian Sombke. MSE calculation is done using only points in the signal region. This is the cause for the discrepancy. We should add a test for this.*

### Show SFR and WSHW as rects

`plot_ws()` and `plot_sfr()` should show the SFR and WSHW as rectangles instead of border lines. This is more intuitive and allows to see the width of the SFR and WSHW.

*Done: Tue Apr 1 2025. Branch: feat14x. PR: #13.*

### Add Getting Started to Reference

Add a function `Getting_Started()` or `get_started()` to the package that contains a link to the online documentation. This should be the first function in the reference manual (if possible).

*Done: Mon Mar 31 2025. Branch: feat14x. PR: #13.*

### Fix unsafe calls

With the 1.4.0 Release we get the following R CMD check notes:

Found the following possibly unsafe calls:
- In `test.R`: `unlockBinding("assert", ns)`
- In `util.R:699:is_list_of_nums`: no visible global function definition for `returns`
- In `test.R:283:not_cran`: no visible binding for global variable `x`

https://github.com/spang-lab/metabodecon/actions/runs/14069618152/job/39400480535

*Done: Fri Mar 28 2025. Branch: feat14x. PR: #13.*
