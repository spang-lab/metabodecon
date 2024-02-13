# data_path (string): path to the parent folder of where the spectra are stored
# e.g.
data_path <- c("C:/Users/Gronwald/Metabolomics/Statistics/Deconvolution/Rechnungen_wolfram/AKI_data/mydata")
# si_size_real_spectrum (positive integer): how many points were used to
# process the real spectrum (often called "si" inside NMR software)
si_size_real_spectrum <- 131072
# scale_factor (positive integer): A factor which is used to avoid rounding
# errors due to numbers becoming too low for R to handle 1000 is a good value
scale_factor <- 1000
# maxShift (positive integer): maximum number of points along the "ppm-axis"
# which a value can be moved by speaq package e.g. 50
maxShift <- 50
# range (positive integer): amount of adjacent columns which are permitted to
# be used for improving the alignment e.g. 5
range <- 5
# lower_bound (positive integer): amount of columns that need to be skipped
lower_bound <- 1


# 1: Use MetaboDecon1D to convert spectrum into lorentz curves
spectrum_data <- generate_lorentz_curves(
  data_path = data_path,
  file_format = "bruker",
  make_rds = FALSE
)

# 2: Look up the global max and minimum ppm values
ppm_range <- get_ppm_range(
  spectrum_data = spectrum_data
)

# 3: Generates the matrix of features based on spectrum data saved as txt files
feat <- gen_feat_mat(
  data_path = data_path,
  ppm_range = ppm_range,
  si_size_real_spectrum = si_size_real_spectrum,
  scale_factor_x = scale_factor
)

# 4: Start alignment by using speaq package
after_speaq_mat <- speaq_align(
  feat = feat,
  maxShift = maxShift
)

# 5: As alignment is still not perfect we will further optimize it by calling the function "combine_peaks"
aligned_res <- combine_peaks(
  shifted_mat = after_speaq_mat,
  range = range,
  lower_bound = lower_bound,
  spectrum_data = spectrum_data,
  data_path = data_path
)

# The returned results after step 5 contain two matrices `aligned_res$long` and
# `aligned_res$short` where in the short version all columns containing only
# zeros have been removed Furthermore, results will be written into two .csv
# files in your data_path directory "aligned_res_short.csv" and
# "aligned_res_long.csv".
