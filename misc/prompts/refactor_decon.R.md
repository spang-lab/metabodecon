
I want to make the following improvements to `decon.R`:

1. In `find_peaks2()` we currently offer a C version and a pure-R version. The
   pure R version is slower, harder to understand and requires 5 helper
   functions. Therefore I want to completely get rid of the R implementation,
   the helper functions and the corresponding tests. Before we do that we just
   need to ensure that the unit tests really cover all edge cases and then there
   is no real reason to keep the R implementation around any longer.

2. In `deconvolute_spectrum()` we currently have a section starting with the
   comment `# Pick best params from attached grid (when npmax > 0)`. I think
   this should be moved into a seperate helper function, so we can keep
   `deconvolute_spectrum()` simple and linear.

3. In `deconvolute_spectrum()` we moved the Rust-Deconvolution into a helper
   function `decon2_from_rust()`. For conistency, we should move the
   R-Deconvolution into a helper function as well. The helpers should be called
   `deconvolute_spectrum_rust` and `deconvolute_spectrum_r` (not
   `decon2_from_rust()`).

4. `decon2` objects currently have an element `mse` (which is a list containing)
   sub-elements `raw`, `norm`, `sm` and `smnorm` (as defined `class.R`). These
   MSEs are useless and we should get rid of them. So my suggestions is search
   all occurences of `mse$raw`, `mse$norm`, `mse$sm`, `mse$smnorm` and
   `mse[["raw"]]`, `mse[["norm"]]`, `mse[["sm"]]`, `mse[["smnorm"]]` in the
   code-base, delete them, make sure the surrounding code still makes sense and
   remove all mentions of these metrics from the documentation and the
   unit-tests.

5. In `grid_deconvolute_spectrum()` we currently use the function
   `deg_to_grid()` to build a `grid` from the `deg` input. BUT, `deg` should
   aready be the grid of deconvolution parameters. In fact 'DeG' is short for
   'Deconvolution Parameter Grid'. I.e., this function should be removed. It is
   completely sufficient to check whether the provided `deg` has the required
   columns. No further transformations needed. Furthermore, we should store the
   "enriched" `deg`, i.e. the `deg` after columns `ar` and `np` have been added
   under a better name than `grid` in the returned spectrum object. I would
   either keep it named `deg` or something like `degar`, `degenr`, `degplus`,
   `degpp` or `degx`. What's your recommendation regarding the name? To make the
   removal of `deg_to_grid()` and the renaming of `x$grid` to something like
   `x$deg` or `x$degar` possible, we need to (1) find all places where
   `grid_deconvolute_spectra`/`grid_deconvolute_spectrum` is called and make
   sure that a proper `deg` is passed and (2) find all places where `grid`
   element is accessed and patch it to the new name.

6. We should remove `write_parameters_txt()`. It was needed to test deprecated
   functions which have been removed by now.

7. We should check whether `mse` and/or `store_as_rds` are still required after
   above changes. If not, they should be removed.