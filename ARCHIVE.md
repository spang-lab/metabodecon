# DONE 2024

## CHECK-1: Use of DTYPP in load_spectrum

In function `deconvolution` from `MetaboDecon1D_deconvolution.R:121` or the extracted version `read_1r_file`, the y values are read as follows:

```R
int_type <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$DTYPP=")]))
int_size <- if (int_type == 0) 4 else 8
path_1r <- file.path(path, expno, "pdata", procno, "1r")
spec_stream <- file(path_1r, "rb")
on.exit(close(spec_stream), add = TRUE)
y <- readBin(
   spec_stream,
   what = "int",
   size = int_size,
   n = n,
   signed = TRUE,
   endian = "little"
)
```

But in [Bruker_NMR_Data_Formats.pdf](http://www.nmragenda-leiden.nl/static/Bruker_manual_4.1.4/topspin_pdf/Bruker_NMR_Data_Formats.pdf) it says:

> The processing status parameter `DTYPP` defines how the data values are stored. If
> the `DTYPP` is 0 ("int"), the stored value represents a mantissa of the data
> point value, the processing status parameter `NC_proc` is the exponent. In this
> case all data points share the same exponent.
>
> [...]
>
> Their format is given by the parameter `DTYPP`, the byte ordering is given by
> the parameter `BYTORDP`, both may be read from the processing status parameter
> file `procs`.

e.i., we should read y as

```R
path_1r <- file.path(path, expno, "pdata", procno, "1r")
spec_stream <- file(path_1r, "rb")
on.exit(close(spec_stream), add = TRUE)
dtypp <- as.numeric(sub("\\D+", "", procs[startsWith(procs, "##$DTYPP=")]))
ncproc <- as.numeric(strsplit(procs[startsWith(procs, "##$NC_proc=")], "=")[[1]][2])
bytordp <- as.numeric(strsplit(procs[startsWith(procs, "##$BYTORDP")], "=")[[1]][2])
if (dtypp == 0) {
   mantissa <- readBin(
      spec_stream,
      what = "int",
      size = 4,
      n = n,
      endian = if (bytordp == 0) "little" else "big",
      signed = TRUE # this needs to be verified as well
   )
   y <- mantissa * (2 ^ ncproc)
} else {
   stop("Not yet implemented")
}
```

Note: I've checked the [urine/urine_1/10/pdata/10/1r](misc/example_datasets/bruker/urine/urine_1/10/pdata/10/1r) example file. There we have `ncproc == -2`, i.e. using the above described implementation, we would get `y <- mantissa * (2 ^ (-2))` == `y <- mantissa * 0.25`, i.e. the right now the values we use are 4 times larger as they should be. I guess that's not the end of the world, as signals get scaled anyway.

Further useful info from [Bruker_NMR_Data_Formats.pdf](http://www.nmragenda-leiden.nl/static/Bruker_manual_4.1.4/topspin_pdf/Bruker_NMR_Data_Formats.pdf):

> The raw data files `fid` and `ser` contain one dimensional or multi-dimensional
> acquired data, respectively. They consist of a sequence of acquired data point
> values in binary format. The acquisition status parameter `DTYPA` defines, how
> the data values are stored. If the `DTYPA` is "int" the stored value represents
> a mantissa of the data point value, the acquisition parameter NC is the
> exponent. All data points share in this case the same exponent. If `DTYPA` is
> "double", the data points are stored as a double precision 64 bit floating
> number, parameter NC is not used.

Note: there are also images in *Bruker_NMR_Data_Formats.pdf* explaining how to
interpret 64 bit float numbers.

---

2024/06/28: Checked and now traced by issue [FEATURE-9](#feature-9-implement-and-export-read_spectra).

## CHECK-2: SFR calculation

In function `deconvolution` of file [MetaboDecon1D.R](R/MetaboDecon1D.R#1372), the signal free region border in data points is calculated as follows:

```R
# Calculate signal free region
signal_free_region_left  <- (spectrum_length+1)-((ppm_highest_value-signal_free_region[1])/(ppm_range/spectrum_length))
signal_free_region_right <- (spectrum_length+1)-((ppm_highest_value-signal_free_region[2])/(ppm_range/spectrum_length))
```

These versions contain the following two errors:

* Error 1: `spectrum$ppm_nstep` is used instead of `spectrum$ppm_step`. That's not really important, because it only causes a slight shift of the border towards the max value, but the number of data points to the left (or right) of the border stays the same (except for the edge case where the border falls exactly on a datapoint).
* Error 2: The `+ 1` in `(spectrum$n + 1)` is wrong. It should be `- 1`. This causes the signal free border to be shifted two points to the left.

Example code to demonstrate the error:

```R
dp <-  c(7,   6,   5,   4,   3,    2,    1,    0   )
ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1, -4.4)
sfrl_ppm <- 2.0 # i.e. 3 points left
n <- length(dp) # 8 (n) data points, i.e. 7 (n - 1) intervals of size 1.3 in between
ppm_step <- (max(ppm) - min(ppm)) / (n - 1) # 1.3
ppm_nstep <- (max(ppm) - min(ppm)) / n # 1.11
sfrl_dp <- (n + 1) - (max(ppm) - sfrl_ppm) / ppm_nstep
# 8 - (4.7 - 2.0) / 1.11  ==  8 - 2.37  ==  6.63  ==> 1 point left
sfrl_dp_without_error_1 <- sfrl_dp <- (n + 1) - (max(ppm) - sfrl_ppm) / ppm_step
# 8 - (4.7 - 2.0) / 1.3  ==  8 - 2.08  ==  6.92  ==> 1 point left
sfrl_dp_correct <- sfrl_dp <- (n - 1) - (max(ppm) - sfrl_ppm) / ppm_step
# 6 - (4.7 - 2.0) / 1.3  ==  6 - 2.08  ==  4.92  ==> 3 point left
```

It's not a big deal as we have over 100000 data points in our examples and whether we exclude e.g. 20500 or 20502 is not super important, as the user cannot make such a precise estimate anyways based on the shown plot, but it makes the code really confusing to understand.

The correct method of converting ppm to dp is `(sfrl_ppm - min(ppm)) / ppm_step` as done in function `ppm_to_dp` in `00_util.R`.

---

Closed with v1.2 because we have a correct implementation in `ppm_to_dp` in `00_util.R`. However, in `generate_lorentz_curves()` we will continue to use the wrong calculation to stay backwards compatible. The new method will be used in the successor function `deconvolute_spectra()`, tracked by [FEATURE-10](#feature-10-implement-deconvolute_spectra)

## CHECK-3: Water signal calculation

In function `deconvolution` of file [MetaboDecon1D_deconvolution.R](R/MetaboDecon1D_deconvolution.R#443), the water signal border in data points is calculated as follows:

```R
# Remove water signal
water_signal_position <- length(spectrum_x)/2
water_signal_position_ppm <- spectrum_x_ppm[length(spectrum_x_ppm)/2]
# Recalculate ppm into data points
range_water_signal<- range_water_signal_ppm/(ppm_range/spectrum_length)
water_signal_left <- water_signal_position - range_water_signal
water_signal_right <- water_signal_position + range_water_signal
```

This is wrong, because of the following two reasons:

1. We count datapoints from 0 to (n - 1), i.e. taking `length(spectrum_x)/2` is incorrect. It should be `(length(spectrum_x) - 1)/2 ` instead. The reason can be easily seen in the example code below.
2. Taking `spectrum_x_ppm[length(spectrum_x_ppm)/2]` as middle ppm value only works, if `length(spectrum_x_ppm)/2` is an integer, i.e. in case of an even number of data points. As soon as we have an odd number of data points `length(spectrum_x_ppm)/2` will be a float. E.g. for 5 data points, it will be 2.5 and calling `spectrum_x_ppm[2.5]` will give the wrong result.

Example code

```R
dp <-  c(6,   5,   4,   3,    2,    1,    0  ) # median and mean is 3
ppm <- c(4.7, 3.4, 2.1, 0.8, -0.5, -1.8, -3.1) # median and mean is 0.8
n <- length(dp) # 7

# Current wrong code
ws_dp <- n/2 # 3.5
ws_ppm <- ppm[ws_dp] # error

# Correct code
ws_dp <- (n - 1)/2 # 3
ws_ppm <- (max(ppm) + min(ppm)) / 2 # 0.8 and also works for odd ((n - 1) / 2)
```

> [!IMPORTANT] Check result:
> Will be fixed with function `deconvolute_spectra` in [FEATURE-10](#check-10-negative-values-for-estimated-a).

## CHECK-4: Data point format

1. Why do we count data points starting from 0 instead of 1? That makes programming in R complicated, as R starts counting at 1.
2. Why do we use "scaled data points" as x values? That's unintuitive and (as far as I can see) unnecessary.
3. Why do we show large ppm values left and low values right?

> [!IMPORTANT] Check result:
> Discussion with Wolfram lead to the following conclusions:
>
> 1. Reason why MetaboDecon starts counting at 0 is not clear. Probably because Martina has C or python background and was more comfortable with 0-based indexing. For `generate_lorentz_curves` we will keep the 0-based indexing to stay backwards compatible. For `deconvolute_spectra` we will switch to 1-based indexing. Tracked by [FEATURE-10](#feature-10-implement-deconvolute_spectra).
> 2. Intention was to prevent rounding errors, but as far as I can see that's not necessary. For `deconvolute_spectra` we completely remove the scale factor and scaled data point numbers. Tracked by [FEATURE-10](#feature-10-implement-deconvolute_spectra).
> 3. That's the general convention in NMR spectroscopy and stems from the fact that in the early days of NMR the frequency was kept contant and the magnetic field was changed instead. To keep the way the spectra looked like, the frequency had to be shown from high to low values. We will keep this convention.

## CHECK-5: Signal preprocessing

```R
# Remove water signal
for (i in ws$right_dp:ws$left_dp) {
   spectrum$y[i] <- 0.00000001 # <--- why not zero or min(spectrum$y)?
}

# Remove negative values of spectrum by Saving the absolut values
for (i in 1:length(spectrum$y)) {
   spectrum$y[i] <- abs(spectrum$y[i]) # # <--- why abs instead of zero or min(spectrum$y)?
}
```

ToSc:

- WS: check if other zeros. If no, use `min(spectrum$y)`.
- Negatives: leave as is

## CHECK-6: Negative value removal

Why do we "remove negative values" using `abs` instead of setting them to zero or the minimum value of the spectrum? See [R/MetaboDecon1D.R#1781](R/MetaboDecon1D.R#1781)

```R
# Remove negative values of spectrum by Saving the absolut values
for(i in 1:length(spectrum_y)){
   spectrum_y[i] <- abs(spectrum_y[i])
}
```

> [!IMPORTANT] Check result:
> Discussed with Wolfram and decided that it doesn't matter whether we use `abs` for negative removal or something else, as these signals should be filtered out anyways by the algorithm.

## CHECK-7: Peak selection procedure

Currently the [peak selection procedure]((R/MetaboDecon1D.R#1825))

1. Requires two loops in R, which is low
2. Iterates from 1:(m-2) instead of 2:(m-1), which is wrong and might cause the miss of one peak if it's at the very end of the spectrum (which is almost certainly in the signal free region so it's not too bad).

We should use vectorized operations and correct indexing to fix this.

```R
# Peak selection procedure

# Calculate second derivative of spectrum
second_derivative <- matrix(nrow = 2, ncol = length(spectrum_x)-2)
for(i in 2:length(spectrum_x)-1){
   second_derivative[1,i-1] <- spectrum_x[i]
   second_derivative[2,i-1] <- spectrum_y[i-1] + spectrum_y[i+1] -2*spectrum_y[i]
}

# Find all local minima of second derivative
peaks_x <- c()
peaks_index <- c()
second_derivative_border <- ncol(second_derivative)-1
for(i in 2:second_derivative_border){
   if(second_derivative[2,i] < 0){
   if((second_derivative[2,i] <= second_derivative[2,i-1]) & (second_derivative[2,i] < second_derivative[2,i+1])){
      #if(((spectrum_y[i+1] >= spectrum_y[i]) & (spectrum_y[i+1] > spectrum_y[i+2])) | ((spectrum_y[i+1] > spectrum_y[i]) & (spectrum_y[i+1] >= spectrum_y[i+2]))){
      # Add local minima to peak list
      peaks_x <- c(peaks_x, second_derivative[1,i])
      peaks_index <- c(peaks_index, i)
   }
   }
}
```

> [!IMPORTANT] Check result:
> Implemented in function [select_peaks_v2](R/generate_lorentz_curves.R#1325).

## CHECK-8: Fix names in blood dataset

It should be `Blood` not `Bood`.

> [!IMPORTANT] Check result:
> Fixed in branch `test-glc`

## CHECK-9: Discuss peak selection

Ask Wolfram whether it's ok that the peak selection sometimes selects two peaks for one local maximum.

> [!IMPORTANT] Check result:
> That's ok and by intention.

## FEATURE-1: Use tempdir for example data

Function `download_example_data` should allow users to specify a temp dir instead the usual XDG directory. This is useful to pass CRAN checks as CRAN doesn't allow writing to th user' home directory.

_Done with [1.1.0+d65098c](https://github.com/spang-lab/metabodecon/tree/d65098cf869e8959055b16570b5d41b5a5c9b46b)._


## FEATURE-2: Add minimal example dataset

Add a minimal dataset to the package, so that the user can run the examples without having to download the full example data. The minimal dataset should be smaller than 1MB. Idea: remove every second or third datapoint from two example spectra, this be enough to get below 1MB.

_Canceled because with [1.1.0+d65098c](https://github.com/spang-lab/metabodecon/tree/d65098cf869e8959055b16570b5d41b5a5c9b46b) we introduced caching for `download_example_data`, so there is no need for including the minimal dataset anymore._

## FEATURE-3: Batch Mode

We should have a batch mode, that does all the above steps truly automatically and creates a pdf containing all quality control images. The pdf can be inspected later on and based on the findings the function call can be adjusted.

1. [x] __PR test-glc__: Add test cases for `generate_lorentz_curves`. Just copy from test cases of `MetetaboDecon1D` and adjust a little bit.
1. [x] __PR batch-mode__: Add test cases for `generate_lorentz_curves_v2`. Just copy from test cases of `generate_lorentz_curves` and fix `generate_lorentz_curves_v2` to pass these tests. (Now we have backwards compatibility)
2. [x] __PR batch-mode__: Replace `generate_lorentz_curves` with `generate_lorentz_curves_v2`
3. [x] __PR batch-mode__: Add remaining arguments of `MetaboDecon1D` to `generate_lorentz_curves` and make sure the are passed on correctly
4. [x] __PR batch-mode__: Add argument `ask` to `generate_lorentz_curves`
5. [x] __PR batch-mode__: Adjust `mock_readline` so that it throws an exception if called more often than there are answers
6. [x] __PR batch-mode__: Implement the ask parameter in `generate_lorentz_curves` until all tests pass
7. [x] __PR clear-helpers__: Remove helper functions that are not used anymore

## FEATURE-4: Parallelize

Batch mode can also run in parallel to speed up calculations. Instead of waiting 1h we need to wait 3 or 6 minutes then.

## FEATURE-5: Add test suite

Write test cases for every function to ensure that future updates don't break any existing behaviour. Tests should be run automatically upon pull requests and pushes to main.

## FEATURE-6: Return lambda in hertz

Original Teams Messages:

<!-- /* cSpell:disable */ -->
From Wolfram at `2023/11/08 3:29 PM`

*Hi Tobi, ich würde MetaboDecon gerne noch um ca zwei-drei Zeilen Code erweitern. wir geben für jedes Signal einen lambda Wert an, das ist die halbe Signalbreite auf halber Höhe. Dieser Wert wird in Datenpunkten angegeben. Für viele Anwendungen braucht man diesen Wert aber auch in Hertz d.h. ich würde auch die lamda Werte in Hertz angeben die Umwandlung ist ganz einfach.*

From Wolfram at `2023/11/16 10:50 AM`

Hi Tobi hier ist der Code zur Linienbreitenberechnung in Hz, Achtung ich hab hier immer die gesamte Linienbreite auf halber Höhe berechnet

```R
# Function to compute the linewidh for all detected features in Hz based on the original
# values given in points
# lw_hz = line width in Hz
# sw_hz = spectral width in hz
# dp datapoints
return_path=c("C:/Users/Gronwald/Metabolomics/Statistics/Deconvolution/Rechnungen_wolfram/Alex_Ref_metabolites_in_Water/")
lw_hz <- function(spectrum_data, sw_hz) {
   for (entry in spectrum_data)
   {
      dp <- entry$x_values[1] * 1000 + 1 # multipy  with 1000 (scale factor) and add 1 as last datapoint has index 0
      cat("dp=", dp, "\n")
      # multiply with 1000 (sale factor)
      # multiply with 2 to get full linewidth instead of half-linewidth
      # take abs as original lambda is given in negative values
      entry$lambda_hz<-abs((sw_hz/dp)*entry$lambda*1000*2)
      cat("lambda= ",entry$lambda_hz[5],"\n") # random signal to show that everything works
      entry$rl1<-data.frame(entry$peak_triplets_middle,entry$lambda_hz)
      entry$rl2<-data.frame(entry$rl1,t(entry$integrals))
      ret_nam<-paste(return_path,entry$filename,"_lw.csv", sep="")
      write.csv2(entry$rl2,file=ret_nam)
   }
}
```

2024/07/09: Done in branch `v1.2.0` with commit `5a9ed6ab00d8e641a2aa82d209de6604d86bf9be`.

## FEATURE-7: Improve return value

The complete history of transformations, e.g. "removal of water signal", "removal of negative values" or "smoothing" should be traceable from the return value of `generate_lorentz_curves`. Therefore I propose the following changes to the return value of `generate_lorentz_curves`:

- `urine_1`
   - `number_of_files` --> length of list
   - `filename` --> `path`
   - `x_values` --> remove (Currently this is given in "scaled data points", which is a strange unit)
   - `x_values_ppm` --> `ppm`
   - `y_values` --> `si`
   - `spectrum_superposition` --> ??
   - `mse_normed` --> `mse`
   - `index_peak_triplets_middle` --> `peaks$middle`
   - `index_peak_triplets_left` --> `peaks$left`
   - `index_peak_triplets_right` --> `peaks$right`
   - `peak_triplets_middle` --> remove (see REFACTOR-6 "Single source of Truth")
   - `peak_triplets_left` --> remove (see REFACTOR-6 "Single source of Truth")
   - `peak_triplets_right` --> remove  (see REFACTOR-6 "Single source of Truth")
   - `integrals` --> ??
   - `signal_free_region` --> `sfr`
   - `range_water_signal_ppm` --> `ws`
   - `A` --> ??
   - `lambda` --> ??
   - `x_0` --> ??

For the "??"" values I haven't yet checked how they're calculated, so I cannot really say whether or how they could be named. But for the others, the improved return value could look like this:

- `urine_1`
   - `ppm` (parts per million for each data point)
   - `si` (signal intensity)
      - `raw` (raw)
      - `wows` (without water signal)
      - `wows_woneg` (without water signal and negative values)
      - `wows_woneg_smoothed` (without water signal, negative values and smoothed)
   - `mse` (mean squared error)
   - `peak` (peak triplet indices)
      - `left` (index of left data point)
      - `center` (index of middle data point)
      - `right` (index of right data point)
   - `sfr` (signal free region in ppm)
      - `left` (left border)
      - `right` (right border)
   - `wshw` (water signal half width in ppm)
   - `path` (absolute path to urine_1)
- `urine_2`
   - ...

2024/06/28: Closed without implementation, as we will keep the return value of `generate_lorentz_curves` backwards compatible with `MetabDecon1D` and instead fix it in `deconvolute_spectra`.

## FEATURE-9 Implement read_spectra

Implement and export `read_spectra` and `read_spectrum` in a way that

1. DTYPP is interpreted correctly (see [CHECK-1](#check-1-use-of-dtypp-in-load_spectrum)) and
2. The spectrum width (SW) in Hz as well as the Magnetic Field Strength is returned (see [FEATURE-6](#feature-6-return-lambda-in-hertz)).

This function can then be used to read spectra if a character string is provided as argument to `deconvolute_spectra()` or `generate_lorentz_curves` (see [FEATURE-11](#feature-11-accept-dataframes-in-glc)).

2024/07/03: done in branch `v1.2.0` with commit `e3c35dce9965cf9a3a44383be818a8f5ab1b0c6e`. Note: the Magnetic Field Strength is not returned directly, but can be calculated via function `calc_B()`.

## FEATURE-11: Accept dataframes in GLC

Let `generate_lorentz_curves` accept dataframes as input. This is useful for Maximilians Bachelorthesis and also makes testing easier. If necessary, implement a private wrapper around `read_spectra`, called `read_spectra_glc` that converts the return value of `read_spectra` to a format that can be used by `deconvolute_spectra`.

2024/07/15: done in branch `v1.2.0` with commit `fa3c427cc1680925c6a12a0eab17f14673c6ee0f`.

## FEATURE-14: Provide simulated datasets

Provide simulated datasets from blood spectra

2024/17/15: Done in branch `v1.2.0` with commit `d01706c1f5885b6e965b55c6db53c041c866ed47`

## FEATURE-15: Add lifecycle badges

Add lifecycle badges to all non-stable exported functions. Functions which are exported but should not be used should be marked as "experimental" or (if possible) as "internal" (idea: check where other badges are stored and whether they can be modified).

2024/07/15: Closed. Will be done with [DOC-1: Document whole package](#doc-1-document-whole-package).

## FEATURE-17: Discard output

By default output of `generate_lorentz_curves` should be discarded during parallel phase. Before this phase, print a message that the remaining task might take up to a few minutes and that live output can be disabled by settings `share_stdout = TRUE`, but that the output might be scrambled depending on configuration of the R installation and the operating system.

2024/07/09: Done in branch `v1.2.0` with commit `f9bf57a5e4c7167dfc3231cfe0ee515b40ad12bf`.

## FEATURE-20: Implement deconvolute_blood

Implement function `deconvolute_blood()` which should roughly do the following

```R
deconvolute_blood <- function(cache = TRUE, force = FALSE, ...) {
   rds <- file.path(cachedir(), "deconvolute_blood.rds")
   new <- old <- if (file.exists(rds)) readRDS(rds) else NULL
   if (cache && is.null(old)) {
      path <- download_example_datasets(...)
      new <- deconvolute(path)
   }
   if (!(identical(new, old) || is.null(old))) {
      warning("Cache and deconvolution differ.", immediate. = TRUE)
   }
   if (cache && (is.null(old) || force)) {
      logf("Writing results to cache file %s", rds)
      saveRDS(new, rds)
   }
   return(new)
}
```

*2024-10-01 20:21: Done in branch v1.2.0 with commit 9b7d65d*

## FEATURE-16: Improve multiprocessing

Right now, output gets scrambled because all procs share one stdout. We can fix this by not using parLapply, but distributing the tasks ourselver over the cores and in a mainloop collecting outputs and results.

*Done in branch v1.2 with commit bea9348 at Thu Sep 12 17:14:20 2024 +0200*

## FEATURE-21: Implement deconvolute_spectra

Implement `deconvolute_spectra()` and `deconvolute_spectrum()` which should be the successors of `deconvolute_ispec()` and `deconvolute_ispecs()`. In particular it should:

1. Accept `spectrum` objects as input (as returned by `read_spectra`). See FEATURE-9
2. Use the correct SFR calculation as described in CHECK-2
3. Uses the correct water signal calculation as described in CHECK-3
4. Use 1-based indexing for data points as described in CHECK-4
5. Remove the scale factor and scaled data point numbers as described in CHECK-4

*Done with v1.2 or v1.3 (not sure) somewhere in between Sep24 and Jan25*

## FEATURE-18: Implement plot_spectrum

Implement `plot_spectrum` which should be the successor of the following functions:

- plot_triplets
- plot_lorentz_curves_save_as_png
- plot_spectrum_superposition_save_as_png

All the functionality, i.e.

- showing triplets
- showing individual lorentz curves and
- showing superposition of all lorentz curves
- storing plots as png on disk

should be controllable via function arguments.

*2024-12-11: Done in branch v1.2.0 with commit 44f8f02*

## FEATURE-13: Merge into main

- [x] Fix all R CMD Check findings.
- [x] Merge branch `v1.2.0` into `main`.

*2025-01-16: done with commit d9e7442*

## CRAN-10: Resubmit to CRAN

*2025-01-16: done with commit ff7d9ea*

## REFACTOR-1: Combine load functions

Combine `load_jcampdx_spectrum` and `load_bruker_spectrum` into one function, which calls `read_jcampdx` or `read_bruker` depending on the `type` argument. The `read_*_spectrum` function should return the measured signal strengths as vector `y_ss` and the corresponding ppm values as vector `x_ppm`. All other elements returned by the `load_*_spectrum` functions can be calculated from those. This makes the code more maintainable and easier to understand.

Useful info for reading bruker files: according to `Bruker_NMR_Data_Formats.pdf` (available through Google), the text files `acqu?` and `proc?` contain acquisition and processing parameters. Files ending with `s` (`acqus`, `proc2s`, ...) describe the status of the dataset. Other files (`acqu`, `proc2`, ...) contain parameter values used in later processing or acquisition steps. Format of all parameter files corresponds to the JCAMP-DX standard, which allows the inclusion of vendor specific parameters by prefixing them with the character sequence `##$`. For this reason, all TopSpin parameters in the file are preceded by this sequence.

## REFACTOR-2: Improve Text Output

The output should be improved (-License, +Timestamps). The License should not be printed after every function execution, unless there is a strong reason to do so. Timestamps should be added to the output, so the user automatically has a rough idea how long the function will take to finish.

## REFACTOR-4: Plotting speed

Function `plot_lorentz_curves_save_as_png` is suuuuper slow. We should try to make this quicker.

2024/06/28: Closed without implementation, as the function was never exported and is now also properly marked with `noRd`.

## REFACTOR-5: Speedup smoothing


Currently, [smoothing of the spectra](R/MetaboDecon1D.R#1797) is done in a for-loop in R, which is slow. We can speed this up by using the `filter` function, which is implemented in C and therefore much faster.

Implemented in function `smooth_signals_v2` in file [generate_lorentz_curves.R](R/generate_lorentz_curves.R). Unfortunately, results are slight different, due to numeric differences between the two methods, ie.

```R
s1 <- smooth_signals_v1(spec)
s2 <- smooth_signals_v2(spec)
identical(s1, s2) == FALSE
all.equal(s1, s2) == TRUE
```

## REFACTOR-6: Use a single unit ony

Function `generate_lorentz_curves` and `MetaboDecon1D` use different units for their calculations and in their returned list, e.g.
 `ppm` (parts per million), `dp` (data points) and `sdp` (scaled data points) as x values and `si` (signal intensity), `ssi` (scaled signal intensity) as y values. Thats not good, as each conversion introduces numeric rounding errors and whenever we do a transformation, e.g. "removal of water signal", "removal of negative values" or "smoothing", we need to do the transformation for all units. Instead, we should use only one unit as single source of truth and provide conversion functions for the user, e.g. `ppm_to_hz` or `ppm_to_dp`. In the final returned list, it's ok to have all units, but at least during the calculation we should stick to `ppm` and `si` I think.

This also makes input for the user much easier, because something like the following wouldn't occur: in function [deconvolution](R/MetaboDecon1D.R#1245) the argument `signal_free_region` is interpreted as `ppm` if `isFALSE(same_parameter) || current_filenumber == 1`. Else, the argument `signal_free_region` is interpreted as `sdp`. This is confusing and error prone.

2024/06/30: Tracked by [FEATURE-10](#feature-10-implement-deconvolute_spectra) instead.

## REFACTOR-7: Split up big functions

The original `MetaboDecon1D` function has approx. 350 lines of code and the original `deconvolution` function has approx 1000 lines of code. This is way too much for testing and modifying. We should extract copy-pasted parts into indivual functions and test these functions separately.

2024/06/30: Done in function `generate_lorentz_curves` in branch `v1.2.0`.

## REFACTOR-8: Improve Metabodecon1D docs

Original Teams Message from Wolfram from `2023/11/08 3:29 PM`

<!-- /* cSpell:disable */ -->
*Hi Tobi, es kommen ja in letzter Zeit eine ganze Menge Fragen zu MetaboDecon auch nach den Variablen im Output. Ich hab hier eine Zusammenstellung dieser Variablen gemacht, vielleicht können wir das auch noch einbauen? Viele Grüße Wolfram*
<!-- /* cSpell:enable */ -->

Output variables of MetaboDecon1D (These variables will be obtained for each analyzed spectrum):

- __Number_of_files:__ number of analyzed spectra (In case folder contains more than one spectrum value >1).
- __Filename:__ name of analyzed spectrum
- __X_values:__ all data points of original spectrum are numbered in descending order. The first data point has the maximum value and the last one the value 0. If we have for example 131072 (128k) datapoints, the first has the value 131.071. Note that numbers were divided by the scale factor which has a default value of 1000 for the x_axis.
- __X_values_ppm:__ as above but here the corresponding ppm values are provided.
- __Y_values:__ intensities of original datapoints
- __Spectrum_superposition:__ y-values of the superposition of all estimated Lorentzcurves, values are provided in a point-wise manner.
- __Mse_normed:__ final mean-squared-error after optimization of all Lorentz curves.
- __Index_peak_triplets_middle:__ each identified signal is defined by three data points. Here the middle data point is given. Note, numbering starts in ascending order from left to right.
- __Index_peak_triplets_left:__ as above for the left data point.
- __Index_peak_triplets_right:__ as above for the right data point.
- __Peak_triplets_middle:__ for each identified signal the position of the middle data point is given in ppm, order from left to right.
- __Peak_triplets_right:__ as above for the right data point
- __Peak_triplets_left:__ as above for the left datapoint.
- __Integrals:__ integrals of deconvoluted signals from left to right as above.
- __Signal_free_region:__ e.g., 109.03619, 21.79452, left and right of these borders no signals are expected. Values are in points like x_values but here with 5 instead of 3 decimals.
- __Range_water_signal_ppm:__ half width of water signal region i.e., where no signals should be identified in ppm, for example 0.15 ppm.
- __A:__ -A*p area under Lorentz-curve, see also integrals, A is always provided as negative number.
- __Lambda:__ for all identified signals half width at half height. Is provided as negative value in data points divided by scale-factor i.e., 1000. For example, a value of von 0.00525 corresponds to 5.25 data points. With a spectral width 12019 Hz (this example) and 31072 data points this corresponds to a half linewidth at half height of 0.48 Hz.
- __X_0:__ center of all estimated Lorentz curves. Provided in data points divided by scale factor (see also x_values).
- __Scale_factor:__ scale factor for x- and y-axis to reduce numerical instabilities, default 1000 and 1000000. Toy example for TSP signal (note numbers will differ for each spectrum): TSP is signal 979, `index_peak_triplets_middle[979]=96955`, `x_0[979]=34.11715` (Note `(131072-96955)/1000=34.117)`, `peak_triplets_middle=0.000` ppm, `lambda[979]=-0.00525` corresponds to 0.48 Hz. `A[979]=-1.218312`; `A* Pi=-3.82`, `integrals[979]=3.82`)

*2024-06-30: Done with commit ab20d64*

## REFACTOR-9: Replace glc calls

Replace all `glc()` calls with calls to `generate_lorentz_curves()`.

*2024-10-01 17:43: Done in branch v1.2.0 with commit 0b52023*

## REFACTOR-10: Replace all md1d calls

1. Implement a function `get_MetaboDecon1D_answers` that takes the path to the spectra as well as the required `sfr`, `wshw` values as as input and returns a vecotr with the corresponding answers to the questions asked by `MetaboDecon1D`.

2. Then replace all `md1d` calls with code snippets as shown below:

   ```R
   answers <- get_MetaboDecon1D_answers(path, sfr = c(3.5, 3.4), wshw = 0)
   decons <- evalwith(answers = answers, MetaboDecon1D(...))
   ```

   This makes it directly visible how cumbersome it is to use the old function and also can be applied to any input folder (in contrast to the current `md1d` function).

*2024-10-07 09:21:34: Done in branch v1.2.0 with commits 18db936, 8f01fae and 6bdaa6f*

## REFACTOR-11: Implement calc_prarp

Implement a function `calc_prarp` that takes a `decon` object and optionally a `truepar` object. If `truepar` is not given, it shall be taken from `decon$meta$simpar`. The function then calculates the PRARP (peak ratio area ratio product) from it, by comparing the estimated parameters with the true parameters.

See function `check_decon_quality()` for existing code to reuse.

*Done in 2024/10/08 in branch v1.2.0 with commit 1b7b4c1*

## REFACTOR-12: Write compliance tests

Write testcases for the following functions to check whether they produce results that are compliant with `MetaboDecon1D()`:

- [x] `deconvolute_ispec()`
- [x] `deconvolute_ispecs()`

*Done in 2025/01/13 in branch v1.2.0 with commit 8c4a16b..c6c6e6a*

## REFACTOR-13: Write prarp tests

Write testcases for the following functions that test for a good PRARP as well as for the correct return type:

- [x] `MetaboDecon1D()`
- [x] `generate_lorentz_curves()`
- [x] `deconvolute()`
- [x] `deconvolute_ispec()`
- [x] `deconvolute_ispecs()`

*Done in 2025/01/13 in branch v1.2.0 with commit 8c4a16b..c6c6e6a*

## CRAN-0: Improve package title

Omit the redundant "Functions for" in your title.

## CRAN-1: Improve package description

Do not start the description with "Functions for", "This package", package name, title or similar.

## CRAN-2: Explain acronyms like NMR

Always explain all acronyms in the description text. e.g.: NMR

## CRAN-3: Correct description format

Write references in the description of the DESCRIPTION file in the form `authors (year) <doi:...>` with no space after 'doi:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

## CRAN-4: Always Explain return value

Please add `\value` to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. `\value{No return value, called for side effects}` or similar). Missing Rd-tags in up to 11 .Rd files, e.g.: `combine_peaks.Rd: \value`, `dohCluster.Rd: \value`, ...

## CRAN-5: Remove unexported examples

You have examples for unexported functions. Please either omit these examples or export these functions. Examples for unexported function with example: plot_spectrum_superposition_save_as_png().

## CRAN-6: Fix vignettes

In addition, we see: "Unexecutable code in vignettes/metabodecon.Rmd": the `#` should be before `"2)"` instead of afterwards, I guess.

## CRAN-7: Check dontrun examples

Remove `dontrun` from examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing in < 5 sec. Reason: `\dontrun{}` should only be used if the example really cannot be executed by the user, e.g. because of missing additional software, missing API keys, etc. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Alternative: You could also replace `\dontrun{}` with `\donttest`, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.

*Done in 2025/01/13 in branch v1.2.0 with commit 8c4a16b..c6c6e6a*

## CRAN-8: Ask before writing to disk

Please ensure that your functions do not write by default or in your `examples/vignettes/tests` in the user's home filespace (including the package directory and `getwd()`). This is not allowed by CRAN policies. Please omit any default path in writing functions. In your examples/vignettes/tests you can write to `tempdir()`.

ToSc: affected functions are:

- `generate_lorentz_curves`
- `plot_triplets`

2024/06/28: Implemented with f5d63f7, 6491b42 and 76b8c4d.

## CRAN-9: Do not change global state

 Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. E.g. `R/MetaboDecon1D.R`. If you're not familiar with the function, please check `?on.exit`. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

## FIX-1: Fix crashes when high smoothing

> From: Maximilian Sombke <Maximilian.Sombke@stud.uni-regensburg.de>
> Sent: Friday, July 12, 2024 2:59 PM
> To: Tobias2 Schmidt <Tobias2.Schmidt@klinik.uni-regensburg.de>
> Subject: Re: Benchmark Metabodecon>
>
> Hey,
>
> sorry für die ganzen mails. Ich wollte es nur so früh wie möglich
> kommunizieren, falls ich keine Lösung finde. Ich hab glaube ich sogar etwas
> sehr Interessantes gefunden worüber wir nächste Woche nochmal reden sollten.
> Als kurze Zusammenfassung:
>
> Meine Vermutung, dass es daran liegt, dass keine Peaks gefunden werden ist
> nicht richtig. Das tatsächliche Problem ist, dass keine Peaks herausgefiltered
> werden. Hier mal ein output log:
>
> ```
> 2024-07-12 14:54:40.96 Starting deconvolution of 1 spectra with 1 core
> 2024-07-12 14:54:40.96 Starting deconvolution of sim_6
> 2024-07-12 14:54:40.96 Removing water signal
> 2024-07-12 14:54:40.96 Removing negative signals
> 2024-07-12 14:54:40.96 Smoothing signals
> 2024-07-12 14:54:41.06 Starting peak selection
> 2024-07-12 14:54:41.06 Detected 4 peaks
> 2024-07-12 14:54:41.06 Removing peaks with low scores
> 2024-07-12 14:54:41.06 Removed NA peaks
> 2024-07-12 14:54:41.06 Initializing Lorentz curves
> Error in xl[ds] <- 2 * xc[ds] - xr[ds] (deconvolution.R#63): NAs are not allowed in subscripted assignments
> ```
>
> Es findet also 4 peaks aber keiner davon scheint einen score zu haben der
> niedrig genug ist um herausgefiltered zu werden, wodurch das dann einfach NA
> ist und das Program crashed. Das passiert hauptsächlich wenn smoothing einen
> gewissen Threshold überschreitet, z.B. smoothing iterations 20 und smoothing
> window size 29, bzw. 15 und 39, oder 25 und 25. Delta scheint zu mindest
> keinen Einfluss hierauf zu haben, es war nur delta = 1 die erste Zeile im Grid
> in der sollche hohen smoothing kombinationen erreicht werden. Die zu hohen
> Parameter nicht zu verwenden hat das Problem (zu mindest bis jetzt) gelöst.
>
> - Sombke Maximilian

2024/07/09: Done in branch `v1.2.0` with commit `f5c204cab44b838afe5d5e8c7ace8c74f11b293c`.

## DOC-1: Document whole package

Document the whole package in vignettes, including chapters about:

- [x] Classes.Rmd
- [x] Compatibility.Rmd
- [x] Contributing.Rmd
- [x] Datasets.Rmd
- [x] FAQ.Rmd
- [x] MetaboDecon1D.Rmd
- [x] Get_Started.Rmd

*Done on 2025/01/14 in branch v1.2.0 with commit 14e0326..71fdc2c*

# DONE 2025-03-26 (PR 12)

## Test install on clean OS

Add a workflow for testing installation on a clean Windows/Linux/Mac OS with R pre-installed, but without R-tools and any packages.

*Done: Mar 10. Branch: mdrb. PR: https://github.com/spang-lab/metabodecon/pull/12.*

## Add Rust Backend Installer

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

## Rename deconvolute_ispec

This is a preparation for issue *Add Rust Backend Argument*. Rename `deconvolute_ispec()` to `deconvolute_spectrum()` and `deconvolute_ispecs()` to `deconvolute_spectra()`.

*Done: Mar 18-22. Branch: mdrb. PR: https://github.com/spang-lab/metabodecon/pull/12.*

## Add Rust Backend Argument

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

# DONE 2025-04-03 (PR 13)

## Check mse calculation in Rust

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

## Show SFR and WSHW as rects

`plot_ws()` and `plot_sfr()` should show the SFR and WSHW as rectangles instead of border lines. This is more intuitive and allows to see the width of the SFR and WSHW.

*Done: Tue Apr 1 2025. Branch: feat14x. PR: #13.*

## Add Getting Started to Reference

Add a function `Getting_Started()` or `get_started()` to the package that contains a link to the online documentation. This should be the first function in the reference manual (if possible).

*Done: Mon Mar 31 2025. Branch: feat14x. PR: #13.*

## Fix unsafe calls

With the 1.4.0 Release we get the following R CMD check notes:

Found the following possibly unsafe calls:
- In `test.R`: `unlockBinding("assert", ns)`
- In `util.R:699:is_list_of_nums`: no visible global function definition for `returns`
- In `test.R:283:not_cran`: no visible binding for global variable `x`

https://github.com/spang-lab/metabodecon/actions/runs/14069618152/job/39400480535

*Done: Fri Mar 28 2025. Branch: feat14x. PR: #13.*

# DONE 2025-04-13 (PR 14)

## Remove Remotes field from DESCRIPTION

When you submit your package to CRAN, all of its dependencies must also be available on CRAN. For this reason, release() will warn you if you try to release a package with a Remotes field.

## Check missing pkgs in align

Installation via `install.packages("metabodecon")` does not install `MassSpecWavelet` and `impute`. So if a user doesn't copy paste the installation instructions but installed via `install.packages("metabodecon")`, these dependencies will be missing. In such scenarios, we should print an error message with the required install commands and abort.

## Change verbose defaults

Change verbose argument for deconvolute and align to TRUE.

## Shrink test-install workflow

Currently the test-install workflow is split over three jobs, with huge amounts of code copy pasted. Extract the R commands into a single `test-install.R` R script, that can be called as `Rscript -e test-install.R <method>` and,

1. Deletes all previously available dependencies incl. Rtools on Windows (i.e. it must be removed from the PATH)
2. Does the installation according to the `method` commandline argument, e.g. "CRAN-Modern", "CRAN-Old" or "Github".

The corresponding commands for "CRAN-Modern", "CRAN-Old" or "Github" can be take from the current version `test-install.yaml`.

After ou created `test-install.R`, update the workflow to use it.

# DONE 2025-04-17 (PR 15)

## Improve questions

1. Question `Signal free region correctly selected? (y/n)` should be replaced by `Borders to Signal Free Regions (green) correctly selected? (y/n)`
2. Question `Water artefact fully inside red vertical lines? (y/n)` should be replaced by `Water artefact fully inside blue area? (y/n)`

## Add function get_si_mat

Add a function `get_si_mat()` for extracting a matrix of signal intensities (SI) from a metabodecon object. The type of returned SI should be `raw` for `spectra`, `sup` for `decons` and `al` for `aligns`.

## Add authors to functions

Add author information all functions (exported and unexported).

## Add lifecycle badges to functions

Add lifecycle badges to all exported, non-stable functions. I.e., add one of the following code blocks at the end of the function description:

```R
#' Superseded by [FUNCTION()] since metabodecon version X.X.X.
#' Will be replaced with metabodecon version 2.0.0.
#'
#' `r lifecycle::badge("deprecated")`
#'
```

or

```R
#'
#' `r lifecycle::badge("experimental")`
#'
```
