# 1. DONE 2025-04-XX (PR 16)

## 1.1. Refactor integral calculations

We should use `A * pi` everywhere unless `bwc = 0`.

# 2. PLANNED

## 2.1. Get codecov to 100%

Make sure all functions are tested or excluded from testing with a corresponding comment giving a reason for the exclusion.

## 2.2. Refactor mse calculations

We should use `mse()` everywhere.

## 2.3. Warn user if peaks in SFR

If delta is small (e.g. 1), peaks in SFR might not be filtered out. Implement this and maybe warn user about it (this is a strong indication that delta was chosen too small).

Note 23.1.2025: in Rust implementation peaks found in IGNORE-REGION are already filtered out automatically before parameter approximation AND points in IGNORE-REGIONS do not contribute to MSE and/or PRARP.

## 2.4. Use R-Universe for mdrb_install

Currently, `mdrb_install()` always installs from source from Github. Installing from R-Univserse should be much faster, as we can install pre-compiled binary packages or at least pre-processed source bundles.

If the user chooses type = "binary", we can also skip the checks for `cargo` and `rustc`.

## 2.5. Refactor parameter approximation

Implement [Max' parameter approximation algorithms](https://gitlab.spang-lab.de/bachelorthesis/ws2425_msombke_metabodecon-v2/-/blob/main/new_method_docs/main.R
) in `calc_A`, `calc_lambda` and `calc_w`.

Probably also solves
[Remove unneeded checks](#remove-unneeded-checks) and
[Check negative A values](#check-negative-a-values).

## 2.6. Improve Sap Dataset

Sap spectra should be simulated as follows:

1. Find three metabolites `mets` related to diabetes with only 1-3 peaks each.
2. Look up their signal centers `x0_` and halfwidths `lambda_`.
3. Define two groups: `healthy` and `diabetes`.
4. Draw signal areas for each sample from a distribution. The distribution should be different for `healthy` and `diabetes`.
5. Draw metabolite shifts for each sample from a normal distribution.
6. Apply the shifts to the signal centers.
7. Simulate data as usual using `simulate_spectra`
8. Make sure, all of the above information is stored inside `$meta$simpar`

## 2.7. Remove generate_lorentz_curves from docs

Replace all occurences of `generate_lorentz_curves()` with `deconvolute()` in the vignettes and articles.

## 2.8. Emit deprecatian warnings

1. Exported functions that will be removed from the package, should emit a deprecation warning when called.
2. Exported functions that will be made private, should check the environment of the calling function:
   - If the environment of the caller is the package namespace, nothing should be done.
   - Else, a deprecation warning should be emitted.

For details see https://lifecycle.r-lib.org/articles/communicate.html

This should be done with a new minor release, e.g. `1.5.0`.

## 2.9. Analyze runtime

Do a benchmark about which parts take up the most time in `deconvolute` and `align` and create issues for improving them. Some of the slow parts should already be mentioned in issues. If this is the case, increase the priority of these issue.

## 2.10. Write paper

Reformat the vignettes as paper and send to Wolfram for proofreading.

# 3. CONDITIONAL

## 3.1. Improve mse_normed calculation

Do after [Analyze Runtime](#analyze-runtime).

In function `add_return_list`:

1. Make the following part faster (or remove completely):

   `s <- sapply(x, function(x_i) sum(abs(A * (lambda / (lambda^2 + (x_i - w)^2))))) # takes approx. 2.2 seconds for urine_1`

2. Return `mse_normed_raw` in addition to `mse_normed` (which is calculated based on `y <- spec$y_smooth`). `mse_normed_raw` should be based on `y <- spec$y_raw`.

## 3.2. Improve Get-Started vignette

Do after [Add function get_si_mat](#add-function-get_si_mat).

In alignment part: add a code snippet showcasing usage of `get_si_mat()` .
In deconvolution part: replace `generate_lorentz_curves()` with `deconvolute()`.

## 3.3. Check negative A values

Do after [Refactor parameter approximation](#refactor-parameter-approximation).

Check why there are negative values for the estimated lorentz curve area A.

## 3.4. Remove unneeded checks

Do after [Refactor parameter approximation](#refactor-parameter-approximation).

Check special handling for cases with A[i] == 0 and lambda[i] == 0 in parameter approximaton. Max analyzed it and concluded that checks are not necessary. My thought: copy paste artifacts from the the w[i] == 0 check (which is wrong).

## 3.5. Write Deconvolution-Details Vignette

Do after [Write paper](#write-paper).

## 3.6. Write Alignment Details Vignette

Do after [Write paper](#write-paper).

# 4. BACKLOG

## 4.1. Test deconvolute with huge nfit

## 4.2. Test deconvolute with peak at zero ppm

It's unclear why <https://github.com/spang-lab/metabodecon/blob/v1.2.0/R/21_decon_helpers.R#L315> checks for w[i] == 0.

Assumption: it's a bug and leads to this peak being missed. If this is the case, remove the check from the code as well.

## 4.3. Check spectrum type

When reading spectra, the spectrum type should be checked (1D, 2D, etc. and if ncessary an error message should be printed)

## 4.4. Integrate jcampdx functions

Integrate jcampdx reader and writer functions from Max.

## 4.5. Remove ispec and idecon

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

## 4.6. Turn articles into vignettes

Make sure that all suitable articles are included as vignettes and built from scratch so they don't take up too much space.

## 4.7. Improve SFR and WS defaults

Replace the default values `wshw = 0.1527692` and `sfr = c(11.44494, -1.8828)` in `generate_lorentz_curves()` with `wshw = "auto"` and `sfr = "auto"`, which should be calculated as follows:

If `c(11.44494, -1.8828)` is part of the ppm range, use these values, otherwise calculate them as

1. `wshw = 0.01 * width(cs)` (where `0.01` is `round(0.007629452, 2)` and `0.007629452` equals `0.1527692 / 20.0236144338963` which is the width of the default WSHW dividided by the width of the `urine_1` spectrum. I.e., the new calculation would give approximately the same proportion of the spectrum width as the default value.)
2. `sfr = max(cs) - c(1/6, 5/6) * width(cs)`

## 4.8. Update expno/procno defaults

Expno and procno should use NULL as default.
If NULL, the function should look for the first available expno/procno folder.
If there is only one, it should use that one.
If there are multiple, it should use expno=10 and procno=10.
If there are multiple but no 10, it should throw an error.

## 4.9. Calculate mandatory SITs in as_decon2

Elements `wsrm` and `nvrm` are mandatory in `decon2` and should therefor be calculated during `as_decon2`. If additional data is needed to do this from `decon0` and/or `decon1` objects, add an argument that allows users to provide this information.

## 4.10. Plot wsr and sfr in plot_spectrum

WSR and SFR should be plotted by plot_spectrum. As soon as this works, we can create a new issue for replacing the `plot_sfr` and `plot_wsr` functions with calls to `plot_spectrum`.

## 4.11. Show prarp in plot_spectrum

## 4.12. Show peak scores in plot_spectrum
