# Open

## Feature: Batch Mode

We should have a batch mode, that does all the above steps truly automatically and creates a pdf containing all quality control images. The pdf can be inspected later on and based on the findings the function call can be adjusted.

1. [x] __PR test-glc__: Add test cases for `generate_lorentz_curves`. Just copy from test cases of `MetetaboDecon1D` and adjust a little bit.
1. [x] __PR batch-mode__: Add test cases for `generate_lorentz_curves_v2`. Just copy from test cases of `generate_lorentz_curves` and fix `generate_lorentz_curves_v2` to pass these tests. (Now we have backwards compatibility)
2. [ ] __PR batch-mode__: Replace `generate_lorentz_curves` with `generate_lorentz_curves_v2`
3. [ ] __PR batch-mode__: Add remaining arguments of `MetaboDecon1D` to `generate_lorentz_curves` and make sure the are passed on correctly
4. [ ] __PR batch-mode__: Add argument `ask` to `generate_lorentz_curves`
5. [ ] __PR batch-mode__: Adjust `mock_readline` so that it throws an exception if called more often than there are answers
6. [ ] __PR batch-mode__: Add test cases for `generate_lorentz_curves(..., ask = FALSE)` without answers (which should fail)
7. [ ] __PR batch-mode__: Implement the ask parameter in `generate_lorentz_curves` until all tests pass
1. [ ] __PR clear-helpers__: Remove helper functions that are not used anymore

## Feature: Parallelize

Batch mode can also run in parallel to speed up calculations. Instead of waiting 1h we need to wait 3 or 6 minutes then.

## Feature: add test suite to ensure correct behaviour after updates

Write test cases for every function to ensure that future updates don't break any existing behaviour. Tests should be run automatically upon pull requests and pushes to main.

## Fix: generate_lorentz_curves should not write to input folders by default

Function `generate_lorentz_curves` should not write to input folders by default. Instead, all generated output files should be stored inside folder `${cwd}/metabodecon_output` by default with the option to change this path.

## Refactor: Text Output

The output should be improved. The License should not be printed after every function execution, unless there is a strong reason to do so. Timestamps should be added to the output, so the user automatically has a rough idea how long the function will take to finish.

## Refactor: Plotting defaults

Function `plot_triplets` should not store to file by default. If we offer an option to store to a file for convenience, it shouldn't be png and we should print the file path.

## Refactor: Plotting speed

Function `plot_lorentz_curves_save_as_png` is suuuuper slow. We should try to make this quicker.

## Fix CRAN review finding 0

Omit the redundant "Functions for" in your title.

## Fix CRAN review finding 1

Do not start the description with "Functions for", "This package", package name, title or similar.

## Fix CRAN review finding 2

Always explain all acronyms in the description text. e.g.: NMR

## Fix CRAN review finding 3

Write references in the description of the DESCRIPTION file in the form `authors (year) <doi:...>` with no space after 'doi:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

## Fix CRAN review finding 4

Please add `\value` to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. `\value{No return value, called for side effects}` or similar). Missing Rd-tags in up to 11 .Rd files, e.g.: `combine_peaks.Rd: \value`, `dohCluster.Rd: \value`, ...

## Fix CRAN review finding 5

You have examples for unexported functions. Please either omit these examples or export these functions. Examples for unexported function with example: plot_spectrum_superposition_save_as_png().

## Fix CRAN review finding 6

In addition, we see: "Unexecutable code in vignettes/metabodecon.Rmd": the `#` should be before `"2)"` instead of afterwards, I guess.

## Fix CRAN review finding 7

Remove `dontrun` from examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing in < 5 sec. Reason: `\dontrun{}` should only be used if the example really cannot be executed by the user, e.g. because of missing additional software, missing API keys, etc. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Alternative: You could also replace `\dontrun{}` with `\donttest`, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.

## Fix CRAN review finding 8

Please ensure that your functions do not write by default or in your `examples/vignettes/tests` in the user's home filespace (including the package directory and `getwd()`). This is not allowed by CRAN policies. Please omit any default path in writing functions. In your examples/vignettes/tests you can write to `tempdir()`.

## Fix CRAN review finding 9

 Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. E.g. `R/MetaboDecon1D.R`. If you're not familiar with the function, please check `?on.exit`. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

## Refactor: Improve output description of Metabodecon1D

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
- __Scale_factor:__ scale factor for x- and y-axis to reduce numerical instabilities, default 1000 and 1000000. Toy example for TSP signal (note numbers will differ for each spectrum): TSP is signal 979, `index_peak_triplets_middle[979]=96955`, `x_0[979]=34.11715` (Note `(131072-96955)/1000=34.117)`, `peak_triplets_middle=0.000` ppm, `lambda[979]=-0.00525` corresponds to 0.48 Hz. `A[979]=-1.218312`; `A* Pi=-3.82`, `integrals[979]=3.82`

## Refactor: improve output of metabodecon

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
      dp<-entry$x_values[1]*1000+1 # multipy  with 1000 (scale factor) and add 1 as last datapoint has index 0
      cat("dp=",dp,"\n")
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

## Check: check use of DTYPP in load_spectrum

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

## Check: ppm to dp conversion

In function `deconvolution` of file [MetaboDecon1D_deconvolution.R](R/MetaboDecon1D_deconvolution.R#156), the signal free region border in data points is calculated as follows:

```R
# Calculate signal free region
signal_free_region_left  <- (spectrum_length+1)-((ppm_highest_value-signal_free_region[1])/(ppm_range/spectrum_length))
signal_free_region_right <- (spectrum_length+1)-((ppm_highest_value-signal_free_region[2])/(ppm_range/spectrum_length))
```

The exact equivalent in `deconvolute_spectrum_v2` with a less ambiguous variable naming is:

```R
sfrl_dp <- (spectrum$n + 1) - (spectrum$ppm_max - sfrl_ppm) / spectrum$ppm_nstep
sfrr_dp <- (spectrum$n + 1) - (spectrum$ppm_max - sfrr_ppm) / spectrum$ppm_nstep
```

These versions contain the following two errors:

* Error 1: `spectrum$ppm_nstep` is used instead of `spectrum$ppm_step`. That's not really important, because it only causes a slight shift of the border towards the max value, but the number of data points to the left (or right) of the border stays the same (except for the edge case where the border falls exactly on a datapoint).
* Error 2: The `+ 1` in `(spectrum$n + 1)` is wrong. It should be `- 1`. This causes the signal free broder to be shifted two points to the left.

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

## Check: water signal calculation

In function `deconvolution` of file [MetaboDecon1D_deconvolution.R](R/MetaboDecon1D_deconvolution.R#443), the signal free region border in data points is calculated as follows:

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
2. Taking `spectrum_x_ppm[length(spectrum_x_ppm)/2]` as middle ppm value only works, if `length(spectrum_x_ppm)/2` is an integer, i.e. in case of an even number of data points. As soon as we have an odd number of data points `length(spectrum_x_ppm)/2` will be a float. E.g. for 5 data points, it will be 2.5 and calling `spectrum_x_ppm[2.5]` will throw an exception.

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

## Check: data point format

- Why do we count data points starting from 0 instead of 1? That's really annoying when programming in R, as R starts counting at 1.
- Why do use "scaled data points" as x values? That's super unintuitive.
- Why do we show large ppm values left and low values right?

## Refactor: combine load_xxx_spectrum functions

Combine `load_jcampdx_spectrum` and `load_bruker_spectrum` into one function, which calls `read_jcampdx_spectrum` or `read_bruker_spectrum` depending on the `type` argument. The `read_*_spectrum` function should return the measured signal strengths as vector `y_ss` and the corresponding ppm values as vector `x_ppm`. All other elements returned by the `load_*_spectrum` functions can be calculated from those. This makes the code more maintainable and easier to understand.

Useful info for reading bruker files: according to `Bruker_NMR_Data_Formats.pdf` (available through Google), the text files `acqu?` and `proc?` contain acquisition and processing parameters. Files ending with `s` (`acqus`, `proc2s`, ...) describe the status of the dataset. Other files (`acqu`, `proc2`, ...) contain parameter values used in later processing or acquisition steps. Format of all parameter files corresponds to the JCAMP-DX standard, which allows the inclusion of vendor specific parameters by prefixing them with the character sequence `##$`. For this reason, all TopSpin parameters in the file are preceded by this sequence.

# Done

## Fix: fix name of samples in blood dataset

It should be `Blood` not `Bood`.

Fixed in branch `test-glc`

## Feature: use temp dirs for full example data

Function `download_example_data` should allow users to specify a temp dir instead the usual XDG directory. This is useful to pass CRAN checks as CRAN doesn't allow writing to th user' home directory.

_Done with [1.1.0+d65098c](https://github.com/spang-lab/metabodecon/tree/d65098cf869e8959055b16570b5d41b5a5c9b46b)._


## Feature: add minimal example dataset

Add a minimal dataset to the package, so that the user can run the examples without having to download the full example data. The minimal dataset should be smaller than 1MB. Idea: remove every second or third datapoint from two example spectra, this be enough to get below 1MB.

_Canceled because with [1.1.0+d65098c](https://github.com/spang-lab/metabodecon/tree/d65098cf869e8959055b16570b5d41b5a5c9b46b) we introduced caching for `download_example_data`, so there is no need for including the minimal dataset anymore._
