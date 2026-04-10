# Cluster Based Peak Alignment

Helper function of
[`align()`](https://spang-lab.github.io/metabodecon/reference/align.md).
Should not be called directly by the user.

Rewrite of
[`speaq::dohCluster()`](https://rdrr.io/pkg/speaq/man/dohCluster.html),
compatible with the data format returned by 'generate_lorentz_curves()'
and 'gen_feat_mat()'. The function name "dohCluster" comes from "Do
Hierarchical Clustering" which is part of the Alignment algorithm
proposed by Vu et al. (2011) in <doi:10.1186/1471-2105-12-405>.

Direct usage of this function has been deprecated with metabodecon
version 1.4.3 and will be removed with metabodecon version 2.0.0.

**\[deprecated\]**

## Usage

``` r
dohCluster(X, peakList, refInd = 0, maxShift = 100, verbose = TRUE, method = 2)
```

## Arguments

- X:

  Dataframe of signal intensities from all spectra as returned by
  [`gen_feat_mat()`](https://spang-lab.github.io/metabodecon/reference/gen_feat_mat.md).

- peakList:

  List of peak indices as returned
  [`gen_feat_mat()`](https://spang-lab.github.io/metabodecon/reference/gen_feat_mat.md).

- refInd:

  Number of the reference spectrum i.e. the spectrum to which all
  signals will be aligned to.

- maxShift:

  Maximum number of points a value can be moved.

- verbose:

  Whether to print additional information during the alignment process.

- method:

  Alignment backend. `1` uses
  [`speaq::hClustAlign()`](https://rdrr.io/pkg/speaq/man/hClustAlign.html),
  `2` (default) uses metabodecon's built-in implementation.

## Value

A list containing two data frames `Y` and `new_peakList`. The first one
contains the aligned spectra, the second one contains the aligned
signals of each spectrum.

## Author

2021-2024 Wolfram Gronwald: initial version.  
2024-2025 Tobias Schmidt: refactored initial version.

## Examples

``` r
deps <- c("MassSpecWavelet", "impute")
deps_installed <- sapply(deps, requireNamespace, quietly = TRUE)
if (all(deps_installed)) {
    # 'speaq' requires 'MassSpecWavelet' and 'impute' to be installed
    sim_subset <- metabodecon_file("bruker/sim_subset")
    decons <- generate_lorentz_curves_sim(sim_subset)
    feat <- gen_feat_mat(decons)
    refObj <- speaq::findRef(feat$peakList)
    hclObj <- dohCluster(
        X = feat$data_matrix,
        peakList = feat$peakList,
        refInd = refObj$refInd,
        maxShift = 100,
        verbose = TRUE
    )
    str(hclObj, 1)
}
#> 2026-04-10 07:15:47.27 Running dohCluster with maxShift = 100 on 2 spectra
#> 2026-04-10 07:15:47.27 Aligning spectrum 1/2
#> 2026-04-10 07:15:47.27 Finished dohCluster in 0.0 s
#> List of 2
#>  $ Y           : num [1:2, 1:2048] 0.0122 0.0105 0.0122 0.0106 0.0122 ...
#>  $ new_peakList:List of 2
```
