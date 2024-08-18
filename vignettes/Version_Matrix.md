---
title: "Compatibility
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Compatibility}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Version Matrix

| Topic | Function                                | 0.2 | 1.0 | 1.3 | 2.0 |
| ----- | --------------------------------------- | --- | --- | --- | --- |
| decon | MetaboDecon1D                           | s   | s   | d   | -   |
| decon | calculate_lorentz_curves                | s   | s   | d   | -   |
| decon | plot_lorentz_curves_save_as_png         | s   | s   | d   | -   |
| decon | plot_triplets                           | s   | s   | s   | d   |
| decon | generate_lorentz_curves                 | -   | s   | s   | d   |
| decon | deconvolute_spectra                     | -   | -   | x   | s   |
| ----- | --------------------------------------- | --- | --- | --- | --- |
| align | combine_peaks                           | -   | s   | s   | ?   |
| align | dohCluster                              | -   | s   | s   | ?   |
| align | gen_feat_mat                            | -   | s   | s   | ?   |
| align | get_ppm_range                           | -   | s   | s   | ?   |
| align | speaq_align                             | -   | s   | s   | ?   |
| ----- | --------------------------------------- | --- | --- | --- | --- |
| data  | get_data_dir                            | -   | s   | d   | -   |
| data  | datadir                                 | -   | -   | s   | s   |
| data  | datadir_persistent                      | -   | -   | s   | s   |
| data  | datadir_temp                            | -   | -   | s   | s   |
| data  | download_example_datasets               | -   | -   | s   | s   |

* -: internal or not available in package
* ?: not yet decided
* d: deprecated
* s: stable
* x: experimental

# Feature Matrix

| Feature                                   | BWCᵃ  | F1ᵇ | F2ᵇ | F3ᵇ | Issue      |
| ----------------------------------------- | ----- | --- | --- | --- | ---------- |
| Doesn't write to disk by default          | semiᶜ | xᵉ  | x   | x   | CRAN-8     |
| Doesn't change wd or global opts          | semiᶜ | xᵉ  | x   | x   | CRAN-9     |
| Uses faster peak selection implementation | yesᵈ  |     | x   | x   | CHECK-7    |
| Batch Mode                                | yes   |     | x   | x   | FEATURE-3  |
| Parallelized                              | yes   |     | x   | x   | FEATURE-4  |
| Improved plotting speed                   | yes   |     | x   | x   | REFACTOR-4 |
| Uses micro functions                      | yes   |     | x   | x   | REFACTOR-7 |
| Doesn't show License                      | semiᶜ |     | x   | x   | REFACTOR-2 |
| Prints timestamps                         | semiᶜ |     | x   | x   | REFACTOR-2 |
| Uses correctly scaled raw y values        | no    |     |     | x   | CHECK-1    |
| Uses correct sfr calculation              | no    |     |     | x   | CHECK-2    |
| Uses correct ws calculation               | no    |     |     | x   | CHECK-3    |
| Uses dynamic signal removal               | no    |     |     | x   | CHECK-5    |
| Uses improved return list                 | no    |     |     | x   | FEATURE-7  |
| Uses faster smoothing implementation      | no    |     |     | x   | REFACTOR-5 |
| Does all calculations in the same unit    | no    |     |     | x   | REFACTOR-6 |

- ᵃ BWC == backwards compatible
- ᵇ F1 ==MetaboDecon1D, F2 == generate_lorentz_curves, F3 == deconvolute_spectra
- ᶜ The change is backwards compatible if you ignore global state, such as files written to disk or output printed to STDOUT. However, if a script expects such side effects, it will fail. Therefore we set this to "semi".
- ᵈ The faster implementation also fixes an indexing bug, which in most cases shouldn't have any effect, but in some rare cases it might cause one peak to be missed.
- ᵉ Only true for MetaboDecon1D versions >= v1.0.0