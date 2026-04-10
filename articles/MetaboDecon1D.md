# MetaboDecon1D

## MetaboDecon1D Usage Example

The `MetaboDecon1D` package (v0.2.2) is the predecessor of this package
and can be downloaded from
[uni-regensburg.de/medicine/functional-genomics/staff/prof-wolfram-gronwald/software](https://www.uni-regensburg.de/medicine/functional-genomics/staff/prof-wolfram-gronwald/software/index.html).
We include this usage example in here so we can easily reference old
behaviour and point out added features in `metabodecon` (v1.x).

### Install the package

1.  Open the following link in your browser:
    [uni-regensburg.de/medicine/functional-genomics/staff/prof-wolfram-gronwald/software](https://www.uni-regensburg.de/medicine/functional-genomics/staff/prof-wolfram-gronwald/software/index.html)
2.  Click [MetaboDecon1D: An R-package for the Deconvolution and
    Integration of 1D NMR
    data](https://genomics.ur.de/software/NMR/MetaboDecon1D/)
3.  Download
    [MetaboDecon1D_0.2.2.tar.gz](https://genomics.ur.de/software/NMR/MetaboDecon1D/MetaboDecon1D_0.2.2.tar.gz)
4.  Start an R session and enter command
    `install.packages("C:/Users/tobi/Downloads/MetaboDecon1D_0.2.2.tar.gz", repos = NULL, type = "source")`
    (replace `C:/Users/tobi/Downloads/MetaboDecon1D_0.2.2.tar.gz` with
    the path to the downloaded file on your computer)

![1_website_gronwald_software.png](MetaboDecon1D/install/1_website_gronwald_software.png)![2_index_gronwald_software.png](MetaboDecon1D/install/2_index_gronwald_software.png)![3_downloaded_archive.png](MetaboDecon1D/install/3_downloaded_archive.png)![4_commands.png](MetaboDecon1D/install/4_commands.png)

### Load package

``` r
library(MetaboDecon1D)
```

### Deconvolute one spectrum in Bruker format

``` r
result <- MetaboDecon1D(
    filepath = "load_example_path",
    filename = "example_human_urine_spectrum",
    file_format = "bruker"
)
```

![1_signal_borders.png](MetaboDecon1D/bruker/1_signal_borders.png)![2_water_artefact.png](MetaboDecon1D/bruker/2_water_artefact.png)![3_output.png](MetaboDecon1D/bruker/3_output.png)

### Visualize results and store plots as png

``` r
str(result)
plot_triplets(result)
plot_lorentz_curves_save_as_png(result)
plot_spectrum_superposition_save_as_png(result)
```

![4_plots.png](MetaboDecon1D/bruker/4_plots.png)

4_plots.png

### Deconvolute one spectrum in jcampdx format

``` r
result <- MetaboDecon1D(
    filepath = "load_example_path",
    filename = "example_human_urine_spectrum.dx",
    file_format = "jcampdx"
)
str(result)
```

![1_one_jcamp_dx.png](MetaboDecon1D/jcamp/1_one_jcamp_dx.png)

### Deconvolute multiple spectra in Bruker format

``` r
result <- MetaboDecon1D(
    filepath = "load_example_path",
    file_format = "bruker"
)
str(result)
```

![5_n_spectra_1.png](MetaboDecon1D/bruker/5_n_spectra_1.png)![6_n_spectra_2.png](MetaboDecon1D/bruker/6_n_spectra_2.png)

### Deconvolute multiple spectra in jcampdx format

``` r
jcamp_results <- MetaboDecon1D(
    filepath = "load_example_path",
    file_format = "jcampdx"
)
plot_triplets(jcamp_results)
plot_lorentz_curves_save_as_png(jcamp_results)
plot_spectrum_superposition_save_as_png(jcamp_results)
```

![2_many_jcamp_dxs_1.png](MetaboDecon1D/jcamp/2_many_jcamp_dxs_1.png)![3_many_jcamp_dxs_2.png](MetaboDecon1D/jcamp/3_many_jcamp_dxs_2.png)
