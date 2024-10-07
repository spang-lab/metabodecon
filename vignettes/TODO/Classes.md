---
title: "Metabodecon Classes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Metabodecon Classes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Metabodecon introduces a set of classes to highlight the presence of certain elements in corresponding objects.  The order of elements may vary between different versions of Metabodecon, thus elements should always be accessed by name, for example, using `x$si` or `x[["cs"]]`.

A classes itself are described in section [Class Overview].  The elements associated with a class are described in section [Elements].  A grafical visualization of the classes and their corresponding functions is given in figure [Metabodecon Workflow].

# Class Overview

## spectrum

Represents one NMR spectrum.  Objects of class 'spectrum' have at least elements 1-3 from section [Elements v1.2+].

## decon0

Represents one deconvoluted NMR spectrum stored in the old [MetaboDecon1D()] format.  Objects of class 'decon0' have at least elements 1-18 from section [Elements v0.2+].

## decon1

Represents one deconvoluted NMR spectrum stored in the backwards compatible [generate_lorentz_curves()] format.  Objects of class 'decon1' have all elements from section [Elements v0.2+].

## decon2

Represents one deconvoluted NMR spectrum stored in the new [deconvolute()] format.  Objects of class 'decon2' have all elements from section [Elements v1.2+] with elements `sit$al` and `lcpar$x0al` set to `NULL`.

## almnt

Represents one deconvoluted NMR spectrum for which the individual peaks have been aligned using [align()]. Objects of class 'almnt' have all elements from section [Elements v1.2+].

## Collections

The classes mentioned above represent individual objects, such as a single spectrum, deconvolution, or alignment. However, it is often useful to describe collections of these objects, such as a list of spectra or deconvolutions. Therefore, for each individual class, a corresponding "collection" class is provided. These collection classes are named: [spectra], [decons0], [decons1], [decons2], and [almnts].


# Elements

## Elements v1.2+

1.  `cs`: Vector of chemical shifts (CS) in parts per pillion (ppm). Must be of the same length as `si`.
2.  `si`: Vector of signal intensities (SI) in arbitrary units (au). Element `si[i]` is the signal intensity measured at chemical shift `cs[i]`, i.e. `si` must be of the same length as `cs`.
3.  `meta`: Additional metadata about the spectrum, e.g.:
    - `name`: The name of the spectrum, e.g. `"Blood 1"` or `"Urine Mouse X23D"`.
    - `path`: The path of the file/folder containing the spectrum data. E.g. `"example_datasets/jcampdx/urine/urine_1.dx"` or `"example_datasets/bruker/urine/urine"`.
    - `type`: The type of experiment, e.g. `"H1 CPMG"` or `"H1 NOESY"`.
    - `fq`: Vector of signal frequencies in Hertz (Hz). Must be of the same length as `si` and `cs`.
    - `mfs`: Magnetic field strength in Tesla, e.g. `14.1`.
    - `simpar`: True lorentz curve parameters. List with elements `A`, `lambda` and `x0`. For details see element `lcpar`. Only available if a spectrum has been simulated.
4.  `args`: Deconvolution parameters:
    - `nfit`: The number of fitting iterations.
    - `smopts`: The smoothing parameters used for the deconvolution.
    - `delta`: The threshold used for peak filtering.
    - `sfr`: Borders of the signal free region in ppm.
    - `wsr`: Borders of the water signal region in ppm.
5.  `sit`: Signal Intensities (SI) after applying various transformations:
    1.  `wsrm`: SIs after Water Signal Removal (WSRM),
    2.  `nvrm`: SIs after WSRM and Negative Value Removal (NVRM).
    3.  `sm`: SIs after WSRM, NVRM and smoothing.
    4.  `sup`: SIs as superposition of Lorentz curves.
    5.  `al`: SIs after alignment.
6.  `peak`: Peak triplets found during peak selection:
    - `center`: Indices of peak centers.
    - `left`: Indices of left borders.
    - `right`: Indices of right peak borders.
7.  `lcpar`: Lorentz curve parameters after parameter approximation:
    - `A`: Amplitude parameter.
    - `lambda`: Halfwidth parameter.
    - `x0`: Center parameter.
    - `x0al`: Center parameter after the spectrum has been aligned.
8.  `mse`: Mean squared error (MSE):
    1.  `raw`: MSE between `si` and `sit$sup`
    2.  `norm`: MSE between `si` and `sit$sup`, divided by `sum(sit$sup)`
    3.  `sm`: MSE between `sit$sm` and `sit$sup`
    4.  `smnorm`: MSE between `sit$sm` and `sit$sup`, divided by `sum(sit$sup)`. Equals `mse_normed` in [Elements v0.2+].

## Elements v0.2+

1.  `filename`: Name of the analyzed spectrum.
2.  `x_values`: Scaled datapoint numbers (SDP). Datapoints are numbered in descending order going from N to 0, where N equals the total amount of data points. Scaled data point numbers are obtained by dividing these numbers by the scale factor of the x-axis. I.e., for a spectrum with 131072 datapoints and a scale factor of 1000, the first scale datapoint has value 131.071 and the last one has value 0.
3.  `x_values_ppm`: The chemical shift of each datapoint in ppm (parts per million).
4.  `y_values`: The scaled signal intensity (SSI) of each datapoint. Obtained by reading the raw intensity values from the provided `data_path` as integers and dividing them scale factor of the y-axis.
5.  `spectrum_superposition`: Scaled signal intensity obtained by calculating the sum of all estimated Lorentz curves for each data point.
6.  `mse_normed`: Normalized mean squared error. Calculated as \mjeqn{\frac{1}{n} \sum_{i=1}^{n} (z_i - \hat{z}_i)^2}{1/n * sum((z_i - zhat_i)^2)} where \mjeqn{z_i}{z_i} is the normalized, smoothed signal intensity of data point i and \mjeqn{\hat{z}_i}{zhat_i} is the normalized superposition of Lorentz curves at data point i. Normalized in this context means that the vectors were scaled so the sum over all data points equals 1.
7.  `peak_triplets_middle`: Chemical shift of peak centers in ppm.
8.  `peak_triplets_left`: Chemical shift of left peak borders in ppm.
9.  `peak_triplets_right`: Chemical shift of right peak borders in ppm.
10. `index_peak_triplets_middle`: Datapoint numbers of peak centers.
11. `index_peak_triplets_left`: Datapoint numbers of left peak borders.
12. `index_peak_triplets_right`: Datapoint numbers of right peak borders.
13. `integrals`: Integrals of the Lorentz curves.
14. `signal_free_region`: Borders of the signal free region of the spectrum in scaled datapoint numbers. Left of the first element and right of the second element no signals are expected.
15. `range_water_signal_ppm`: Half width of the water signal in ppm. Potential signals in this region are ignored.
16. `A`: Amplitude parameter of the Lorentz curves. Provided as negative number to maintain backwards compatibility with MetaboDecon1D. The area under the Lorentz curve is calculated as \mjeqn{A \cdot \pi}{A * pi}.
17. `lambda`: Half width of the Lorentz curves in scaled data points. Provided as negative value to maintain backwards compatibility with MetaboDecon1D. Example: a value of -0.00525 corresponds to 5.25 data points. With a spectral width of 12019 Hz and 131072 data points this corresponds to a halfwidth at half height of approx. 0.48 Hz. The corresponding calculation is: (12019 Hz / 131071 dp) * 5.25 dp.
18. `x_0`: Center of the Lorentz curves in scaled data points.
19. `y_values_raw`: The raw signal intensity of each datapoint
20. `x_values_hz`: The frequency of each datapoint in Hz
21. `mse_normed_raw`: Normalized mean squared error when comparing the raw signal intensities with the superposition of Lorentz curves.
22. `x_0_hz`: Center of the Lorentz curves in Hz.
23. `x_0_dp`: Center of the Lorentz curves in data points.
24. `x_0_ppm`: Center of the Lorentz curves in ppm.
25. `A_hz`: Amplitude parameter of the Lorentz curves in Hz.
26. `A_dp`: Amplitude parameter of the Lorentz curves in data points.
27. `A_ppm`: Amplitude parameter of the Lorentz curves in ppm.
28. `lambda_hz`: Half width of the Lorentz curves in Hz.
29. `lambda_dp`: Half width of the Lorentz curves in data points.
30. `lambda_ppm`: Half width of the Lorentz curves in ppm.

<!-- Reference Links -->

[alignment]: #alignment
[alignments]: #alignments
[almnt]: #almnt
[almnts]: #almnts
[Class Overview]: #class-overview
[decon0]: #decon0
[decon1]: #decon1
[decon2]: #decon2
[decon3]: #decon3
[decons0]: #decons0
[decons1]: #decons1
[decons2]: #decons2
[decons3]: #decons3
[Elements]: #elements
[Elements v1.2+]: #elements-v1.2+
[Elements v0.2+]: #elements-v0.2+
[Metabodecon Workflow]: #metabodecon-workflow
[Multiple Object Classes]: #multiple-object-classes
[Single Object Classes]: #single-object-classes
[spectra]: #spectra
[spectrum]: #spectrum
[MetaboDecon1D()]: https://spang-lab.github.io/metabodecon/reference/MetaboDecon1D.html
[generate_lorentz_curves()]: https://spang-lab.github.io/metabodecon/reference/generate_lorentz_curves.html
[deconvolute()]: https://spang-lab.github.io/metabodecon/reference/deconvolute()