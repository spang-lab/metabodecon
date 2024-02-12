# 1. Feature: Batch Mode

We should have a batch mode, that does all the above steps truly automatically and creates a pdf containing all quality control images. The pdf can be inspected later on and based on the findings the function call can be adjusted.

# 2. Feature: Parallelize

Batch mode can also run in parallel to speed up calculations. Instead of waiting 1h we need to wait 3 or 6 minutes then.

# 3. Feature: add test suite to ensure correct behaviour after updates

Write test cases for every function to ensure that future updates don't break any existing behaviour. Tests should be run automatically upon pull requests and pushes to main.

# 4. Feature: add minimal example dataset

Add a minimal dataset to the package, so that the user can run the examples without having to download the full example data. The minimal dataset should be smaller than 1MB. Idea: remove every second or third datapoint from two example spectra, this be enough to get below 1MB.

# 5. Feature: use temp dirs for full example data

Function `download_example_data` should allow users to specify a temp dir instead the usual XDG directory. This is useful to pass CRAN checks as CRAN doesn't allow writing to th user' home directory.

# 6. Fix: generate_lorentz_curves should not write to input folders by default

Function `generate_lorentz_curves` should not write to input folders by default. Instead, all generated output files should be stored inside folder `${cwd}/metabodecon_output` by default with the option to change this path.

# 7. Fix: fix name of samples in blood dataset

It should be `Blood` not `Bood`.

# 8. Refactor: Text Output

The output should be improved. The License should not be printed after every function execution, unless there is a strong reason to do so. Timestamps should be added to the output, so the user automatically has a rough idea how long the function will take to finish.

# 9. Refactor: Plotting defaults

Function `plot_triplets` should not store to file by default. If we offer an option to store to a file for convenience, it shouldn't be png and we should print the file path.

# 10. Refactor: Plotting speed

Function `plot_lorentz_curves_save_as_png` is suuuuper slow. We should try to make this quicker.

# 11. Fix CRAN review finding 0

Omit the redundant "Functions for" in your title.

# 12. Fix CRAN review finding 1

Do not start the description with "Functions for", "This package", package name, title or similar.

# 13. Fix CRAN review finding 2

Always explain all acronyms in the description text. e.g.: NMR

# 14. Fix CRAN review finding 3

Write references in the description of the DESCRIPTION file in the form `authors (year) <doi:...>` with no space after 'doi:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

# 15. Fix CRAN review finding 4

Please add `\value` to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. `\value{No return value, called for side effects}` or similar). Missing Rd-tags in up to 11 .Rd files, e.g.: `combine_peaks.Rd: \value`, `dohCluster.Rd: \value`, ...

# 16. Fix CRAN review finding 5

You have examples for unexported functions. Please either omit these examples or export these functions. Examples for unexported function with example: plot_spectrum_superposition_save_as_png().

# 17. Fix CRAN review finding 6

In addition, we see: "Unexecutable code in vignettes/metabodecon.Rmd": the `#` should be before `"2)"` instead of afterwards, I guess.

# 18. Fix CRAN review finding 7

Remove `dontrun` from examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing in < 5 sec. Reason: `\dontrun{}` should only be used if the example really cannot be executed by the user, e.g. because of missing additional software, missing API keys, etc. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Alternative: You could also replace `\dontrun{}` with `\donttest`, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.

# 19. Fix CRAN review finding 8

Please ensure that your functions do not write by default or in your `examples/vignettes/tests` in the user's home filespace (including the package directory and `getwd()`). This is not allowed by CRAN policies. Please omit any default path in writing functions. In your examples/vignettes/tests you can write to `tempdir()`.

# 20. Fix CRAN review finding 9

 Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. E.g. `R/MetaboDecon1D.R`. If you're not familiar with the function, please check `?on.exit`. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

# 21. Refactor: Improve output description of Metabodecon1D

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
- __Scale_factor:__ scale factor for x- and y-axis to reduce numerical instabilities, default 1000 and 1000000.
Toy example for TSP signal (note numbers will differ for each spectrum): TSP is signal 979, index_peak_triplets_middle[979]=96955, x_0[979]=34.11715 (Note (131072-96955)/1000=34.117), peak_triplets_middle=0.000 ppm, lambda[979]=-0.00525 corresponds to 0.48 Hz. A[979]=-1.218312; A* Pi=-3.82, integrals[979]=3.82

# 22. Refactor: improve output of metabodecon

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
<!-- /* cSpell:enable */ -->
