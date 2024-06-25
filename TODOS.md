# Version Overview

## Affected Functions

- F1: `MetaboDecon1D::MetaboDecon1D`
- F2: `metabodecon::MetaboDecon1D`
- F3: `metabodecon::generate_lorentz_curves`
- F4: `metabodecon::deconvolute_spectra`

## Version Matrix

| Package       | Version     | Release | F1  | F2  | F3  | F4  |
| ------------- | ----------- | ------- | --- | --- | --- | --- |
| MetaboDecon1D | == 0.2.2    | 2021    | e   |     |     |     |
| metabodecon   | >= 1.1, < 2 | 2024    |     | ed  | e   | i   |
| metabodecon   | >= 2.0, < 3 | 2025    |     | id  | ed  | e   |

- i: internal
- e: exported
- d: deprecated

## Feature Matrix

| Feature                                   | F1  | F2  | F3  | F4  | BWC[^3]  | Issue      |
| ----------------------------------------- | --- | --- | --- | --- | -------- | ---------- |
| Doesn't write to disk by default          |     | x   | x   | x   | semi[^2] | CRAN-8     |
| Doesn't change wd or global opts          |     | x   | x   | x   | semi[^2] | CRAN-9     |
| Uses faster peak selection implementation |     |     | x   | x   | yes[^1]  | CHECK-7    |
| Batch Mode                                |     |     | x   | x   | yes      | FEATURE-3  |
| Parallelized                              |     |     | x   | x   | yes      | FEATURE-4  |
| Improved plotting speed                   |     |     | x   | x   | yes      | REFACTOR-4 |
| Uses micro functions                      |     |     | x   | x   | yes      | REFACTOR-7 |
| Doesn't show License                      |     |     | x   | x   | semi[^2] | REFACTOR-2 |
| Prints timestamps                         |     |     | x   | x   | semi[^2] | REFACTOR-2 |
| Uses correctly scaled raw y values        |     |     |     | x   | no       | CHECK-1    |
| Uses correct sfr calculation              |     |     |     | x   | no       | CHECK-2    |
| Uses correct ws calculation               |     |     |     | x   | no       | CHECK-3    |
| Uses dynamic signal removal               |     |     |     | x   | no       | CHECK-5    |
| Uses improved return list                 |     |     |     | x   | no       | FEATURE-7  |
| Uses faster smoothing implementation      |     |     |     | x   | no       | REFACTOR-5 |
| Does all calculations in the same unit    |     |     |     | x   | no       | REFACTOR-6 |


[^1] The faster implementation also fixes an indexing bug, which in most cases shouldn't have any effect, but in some rare cases it might cause one peak to be missed.
[^2] The change is backwards compatible if you ignore global state, such as files written to disk or output printed to STDOUT. However, if a script expects such side effects, it will fail. Therefore we set this to "semi".
[^3] BWC: backwards compatible

# Todos

## Features

### FEATURE-6: Return lambda in hertz

Original Teams Messages:

<!-- /* cSpell:disable */ -->
From Wolfram at `2023/11/08 3:29 PM`

*Hi Tobi, ich würde MetaboDecon gerne noch um ca zwei-drei Zeilen Code erweitern. wir geben für jedes Signal einen lambda Wert an, das ist die halbe Signalbreite auf halber Höhe. Dieser Wert wird in Datenpunkten angegeben. Für viele Anwendungen braucht man diesen Wert aber auch in Hertz d.h. ich würde auch die lamda Werte in Hertz angeben die Umwandlung ist ganz einfach.*

From Wolfram at `2023/11/16 10:50 AM`

Hi Tobi hier ist der Code  zur Linienbreitenberechnung in Hz, Achtung ich hab hier immer die gesamte Linienbreite auf halber Höhe berechnet

```R
# Function to compute the linewidh for all detected features in Hz baed on the original
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

### FEATURE-7: Improve return value

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

Alternative

- `urine_1`
   - `x`
      - `ppm` (parts per million for each data point)
      - `in_sfr` (is in signal free region)
      - `in_wsr` (is in water signal region)
      - `is_plb` (is peak border left)
      - `is_pc` (is peak center)
      - `is_prb` (is peak border right)
   - `y`
      - `raw` (raw signal intensity)
      - `wows` (without water signal)
      - `wows_woneg` (without water signal and negative values)
      - `wows_woneg_smooth1` (without water signal and negative values after first round of smoothing)
      - `wows_woneg_smooth2` (without water signal and negative values after second round of smoothing)
      - ...
      - `wows_woneg_smoothed` (without water signal and negative values after last round of smoothing)
   - `mse` (mean squared error)
   - `sfr` (signal free region)
      - `left` (left border in ppm)
      - `right` (right border in ppm)
   - `ws` (water signal)
      - `hw` (half width)
   - `path` (absolute path to urine_1)
- `urine_2`
   - ...

### FEATURE-8: Warn user if peaks are found in SFR

If delta is small (e.g. 1), peaks in SFR might not be filtered out. Either implement this and warn user about it (this is a strong indication that delta was chosen too small).

## Refactorings (4)

### REFACTOR-4: Plotting speed

Function `plot_lorentz_curves_save_as_png` is suuuuper slow. We should try to make this quicker.

### REFACTOR-6: Use a single unit as source of truth

Function `generate_lorentz_curves` and `MetaboDecon1D` use different units for their calculations and in their returned list, e.g.
 `ppm` (parts per million), `dp` (data points) and `sdp` (scaled data points) as x values and `si` (signal intensity), `ssi` (scaled signal intensity) as y values. Thats not good, as each conversion introduces numeric rounding errors and whenever we do a transformation, e.g. "removal of water signal", "removal of negative values" or "smoothing", we need to do the transformation for all units. Instead, we should use only one unit as single source of truth and provide conversion functions for the user, e.g. `ppm_to_hz` or `ppm_to_dp`. In the final returned list, it's ok to have all units, but at least during the calculation we should stick to `ppm` and `si` I think.

This also makes input for the user much easier, because something like the following wouldn't occur: in function [deconvolution](R/MetaboDecon1D.R#1245) the argument `signal_free_region` is interpreted as `ppm` if `isFALSE(same_parameter) || current_filenumber == 1`. Else, the argument `signal_free_region` is interpreted as `sdp`. This is confusing and error prone.


### REFACTOR-7: Split monolithic functions into smaller parts

The original `MetaboDecon1D` function has approx. 350 lines of code and the original `deconvolution` function has approx 1000 lines of code. This is way too much for testing and modifying. We should extract copy-pasted parts into indivual functions and test these functions separately.

### REFACTOR-8: Improve docs for Metabodecon1D return value

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

### REFACTOR-9: Improve mse_normed calculation

In function `add_return_list_v13`:

1. Make the following part faster (or remove completely):

   `s <- sapply(x, function(x_i) sum(abs(A * (lambda / (lambda^2 + (x_i - w)^2))))) # takes approx. 2.2 seconds for urine_1`

2. Check whether it makes sense to calculate y as `y <- spec$y_smooth` or whether it's better to use `y <- spec$y_raw`

## CRAN (10)

### CRAN-3: Use correct reference format in DESCRIPTION

Write references in the description of the DESCRIPTION file in the form `authors (year) <doi:...>` with no space after 'doi:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

### CRAN-4: Explain return value in function docs

Please add `\value` to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. `\value{No return value, called for side effects}` or similar). Missing Rd-tags in up to 11 .Rd files, e.g.: `combine_peaks.Rd: \value`, `dohCluster.Rd: \value`, ...

### CRAN-5: Remove examples from unexported functions

You have examples for unexported functions. Please either omit these examples or export these functions. Examples for unexported function with example: plot_spectrum_superposition_save_as_png().

### CRAN-6: Fix vignettes

In addition, we see: "Unexecutable code in vignettes/metabodecon.Rmd": the `#` should be before `"2)"` instead of afterwards, I guess.

### CRAN-7: Check dontrun examples

Remove `dontrun` from examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing in < 5 sec. Reason: `\dontrun{}` should only be used if the example really cannot be executed by the user, e.g. because of missing additional software, missing API keys, etc. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Alternative: You could also replace `\dontrun{}` with `\donttest`, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.

### CRAN-8: Functions should not write to disk by default

Please ensure that your functions do not write by default or in your `examples/vignettes/tests` in the user's home filespace (including the package directory and `getwd()`). This is not allowed by CRAN policies. Please omit any default path in writing functions. In your examples/vignettes/tests you can write to `tempdir()`.

ToSc: affected functions are:

- `generate_lorentz_curves`
- `plot_triplets`

### CRAN-9: Functions should not change working dir or global options

 Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. E.g. `R/MetaboDecon1D.R`. If you're not familiar with the function, please check `?on.exit`. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

## Checks (6)

### CHECK-1: Use of DTYPP in load_spectrum

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

### CHECK-2: Signal free region calculation

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

The correct method of converting ppm to dp is `(sfrl_ppm - min(ppm)) / ppm_step` as done in function `ppm_to_dp` in `util.R`.

### CHECK-3: Water signal calculation

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

### CHECK-4: Data point format

- Why do we count data points starting from 0 instead of 1? That makes programming in R complicated, as R starts counting at 1.
- Why do we use "scaled data points" as x values? That's unintuitive and (as far as I can see) unnecessary.
- Why do we show large ppm values left and low values right? ToSc: Convention: keep it.

### CHECK-5: Signal preprocessing

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

- WS: check if other zeros. If no, use min(spectrum$y).
- Negatives: leave as is


### CHECK-6: Negative value removal

Why do we "remove negative values" using `abs` instead of setting them to zero or the minimum value of the spectrum? See [](R/MetaboDecon1D.R#1781)

```R
# Remove negative values of spectrum by Saving the absolut values
for(i in 1:length(spectrum_y)){
   spectrum_y[i] <- abs(spectrum_y[i])
}
```

### CHECK-9: Ask Wolfram about peak selection

Ask Wolfram whether it's ok that the peak selection sometimes selects two peaks for one local maximum.

# Done (6)

## Features

### FEATURE-1: Use temp dirs for full example data

Function `download_example_data` should allow users to specify a temp dir instead the usual XDG directory. This is useful to pass CRAN checks as CRAN doesn't allow writing to th user' home directory.

_Done with [1.1.0+d65098c](https://github.com/spang-lab/metabodecon/tree/d65098cf869e8959055b16570b5d41b5a5c9b46b)._


### FEATURE-2: Add minimal example dataset

Add a minimal dataset to the package, so that the user can run the examples without having to download the full example data. The minimal dataset should be smaller than 1MB. Idea: remove every second or third datapoint from two example spectra, this be enough to get below 1MB.

_Canceled because with [1.1.0+d65098c](https://github.com/spang-lab/metabodecon/tree/d65098cf869e8959055b16570b5d41b5a5c9b46b) we introduced caching for `download_example_data`, so there is no need for including the minimal dataset anymore._

### FEATURE-3: Batch Mode

We should have a batch mode, that does all the above steps truly automatically and creates a pdf containing all quality control images. The pdf can be inspected later on and based on the findings the function call can be adjusted.

1. [x] __PR test-glc__: Add test cases for `generate_lorentz_curves`. Just copy from test cases of `MetetaboDecon1D` and adjust a little bit.
1. [x] __PR batch-mode__: Add test cases for `generate_lorentz_curves_v2`. Just copy from test cases of `generate_lorentz_curves` and fix `generate_lorentz_curves_v2` to pass these tests. (Now we have backwards compatibility)
2. [x] __PR batch-mode__: Replace `generate_lorentz_curves` with `generate_lorentz_curves_v2`
3. [x] __PR batch-mode__: Add remaining arguments of `MetaboDecon1D` to `generate_lorentz_curves` and make sure the are passed on correctly
4. [x] __PR batch-mode__: Add argument `ask` to `generate_lorentz_curves`
5. [x] __PR batch-mode__: Adjust `mock_readline` so that it throws an exception if called more often than there are answers
6. [x] __PR batch-mode__: Implement the ask parameter in `generate_lorentz_curves` until all tests pass
7. [x] __PR clear-helpers__: Remove helper functions that are not used anymore

### FEATURE-4: Parallelize

Batch mode can also run in parallel to speed up calculations. Instead of waiting 1h we need to wait 3 or 6 minutes then.

### FEATURE-5: Add test suite to ensure correct behaviour after updates

Write test cases for every function to ensure that future updates don't break any existing behaviour. Tests should be run automatically upon pull requests and pushes to main.

## Refactorings

### REFACTOR-1: Combine load_xxx_spectrum functions

Combine `load_jcampdx_spectrum` and `load_bruker_spectrum` into one function, which calls `read_jcampdx_spectrum` or `read_bruker_spectrum` depending on the `type` argument. The `read_*_spectrum` function should return the measured signal strengths as vector `y_ss` and the corresponding ppm values as vector `x_ppm`. All other elements returned by the `load_*_spectrum` functions can be calculated from those. This makes the code more maintainable and easier to understand.

Useful info for reading bruker files: according to `Bruker_NMR_Data_Formats.pdf` (available through Google), the text files `acqu?` and `proc?` contain acquisition and processing parameters. Files ending with `s` (`acqus`, `proc2s`, ...) describe the status of the dataset. Other files (`acqu`, `proc2`, ...) contain parameter values used in later processing or acquisition steps. Format of all parameter files corresponds to the JCAMP-DX standard, which allows the inclusion of vendor specific parameters by prefixing them with the character sequence `###$`. For this reason, all TopSpin parameters in the file are preceded by this sequence.

### REFACTOR-2: Text Output (-License, +Timestamps)

The output should be improved. The License should not be printed after every function execution, unless there is a strong reason to do so. Timestamps should be added to the output, so the user automatically has a rough idea how long the function will take to finish.

### REFACTOR-5: Speedup smoothing


Currently, [smoothing of the spectra](R/MetaboDecon1D.R#1797) is done in a for-loop in R, which is slow. We can speed this up by using the `filter` function, which is implemented in C and therefore much faster.

Implemented in function `smooth_signals_v2` in file [generate_lorentz_curves.R](R/generate_lorentz_curves.R). Unfortunately, results are slight different, due to numeric differences between the two methods, ie.

```R
s1 <- smooth_signals_v1(spec)
s2 <- smooth_signals_v2(spec)
identical(s1, s2) == FALSE
all.equal(s1, s2) == TRUE
```

## Checks

### CHECK-7: Peak selection procedure

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

Implemented in function [select_peaks_v2](R/generate_lorentz_curves.R#1325).

### CHECK-8: Fix name of samples in blood dataset

It should be `Blood` not `Bood`.

Fixed in branch `test-glc`

## CRAN

### CRAN-0: Omit "Functions for" in title

Omit the redundant "Functions for" in your title.

### CRAN-1: Omit "Functions for" in DESCRIPTION

Do not start the description with "Functions for", "This package", package name, title or similar.

### CRAN-2: Explain acronyms like NMR

Always explain all acronyms in the description text. e.g.: NMR