library(testthat)

# TEST DESCRIPTION
#
# 1. Get the sap2 spectrum
# 2. Deconvolute the sap2 spectrum using deconvolute_ispec with bwc=2
# 3. Plot and print the spectrum
# 4. Write down the expected values for decon0, decon1 and decon2. It should be
#    easy to deduce the expected values from the plot and printout of step 3.
# 5. Convert sap2_idecon to decon0, decon1 and decon2
# 6. Compare the expected values defined in step 4 with the observed values
#    obtained in step 5.
# 7. Repeat steps 2-6 for bwc=1 and bwc=0 (maybe in a separate test file)

# 1. Get the sap2 spectrum
sap2_spectrum <- get_sap2_ispec()

# 2. Deconvolute the sap2 spectrum using deconvolute_ispec with bwc=2
sap2_idecon <- deconvolute_ispec(
    ispec = sap2_ispec,
    sfr = c(3.2, -3.2),
    smopts = c(0, 3),
    delta = 1
)
if (identical(environment(), globalenv())) {
    str(sap2_idecon, 1)
    plot_spectrum(sap2_idecon, sub_show = FALSE)
}

# Decon1 --> Decon2
expect_equal(names(sap2_decon1_bwc1), decon1_members)
expect_equal(sap2_decon1_bwc1$number_of_files            , 1)
expect_equal(sap2_decon1_bwc1$filename                   , "sap2")
expect_equal(sap2_decon1_bwc1$x_values                   , x$sdp)
expect_equal(sap2_decon1_bwc1$x_values_ppm               , x$cs)


 $ cs   : num [1:128] 6.4 6.3 6.2 6.1 6 5.9 5.8 5.7 5.6 5.5 ...
 $ si   : num [1:128] 162 172 183 198 217 243 282 344 449 610 ...
 $ meta :List of 4
  ..$ name          : chr "sap2"
  ..$ fq            : num [1:128] 6e+08 6e+08 6e+08 6e+08 6e+08 ...
  ..$ simpar        :List of 7
  .. ..$ name  : chr "sap2"
  .. ..$ cs    : num [1:128] 6.4 6.3 6.2 6.1 6 5.9 5.8 5.7 5.6 5.5 ...
  .. ..$ pkr   : num [1:2] 3.2 -3.2
  .. ..$ x0    : num [1:8] 5.4 4.2 2 1.4 0 -1.8 -4 -5.2
  .. ..$ A     : num [1:8] 100 150 2000 500 5000 2500 150 100
  .. ..$ lambda: num [1:8] 0.2 0.2 0.4 0.2 0.6 0.5 0.2 0.2
  .. ..$ noise : num 0
 $ args :List of 10
  ..$ nfit    : num 3
  ..$ smopts  : num [1:2] 0 3
  ..$ delta   : num 1
  ..$ sfr     : num [1:2] 3.2 -3.2
  ..$ wshw    : num 0
  ..$ ask     : logi FALSE
  ..$ force   : logi FALSE
  ..$ verbose : logi FALSE
  ..$ bwc     : num 2
  ..$ nworkers: num 1
 $ sit  :List of 4
  ..$ wsrm: num [1:128] 0.000162 0.000172 0.000183 0.000198 0.000217 0.000243 0.000282 0.000344 0.000449 0.00061 ...
  ..$ nvrm: num [1:128] 0.000162 0.000172 0.000183 0.000198 0.000217 0.000243 0.000282 0.000344 0.000449 0.00061 ...
  ..$ sm  : num [1:128] 0.000162 0.000172 0.000183 0.000198 0.000217 0.000243 0.000282 0.000344 0.000449 0.00061 ...
  ..$ sup : NULL
 $ peak :'data.frame':  8 obs. of  6 variables:
  ..$ left  : num [1:8] 10 22 43 50 62 81 104 116
  ..$ center: int [1:8] 11 23 45 51 65 83 105 117
  ..$ right : num [1:8] 12 24 46 52 68 85 106 118
  ..$ score : num [1:8] 0.000247 0.000368 0.000887 0.00096 0.001079 ...
  ..$ high  : logi [1:8] FALSE FALSE TRUE TRUE TRUE TRUE ...
  ..$ region: chr [1:8] "sfrl" "sfrl" "norm" "norm" ...
 $ lcpar:List of 4
  ..$ A        : num [1:4] -1.80e-05 -6.96e-06 -4.98e-05 -2.52e-05
  ..$ lambda   : num [1:4] -0.0038 -0.00236 -0.00599 -0.00502
  ..$ w        : num [1:4] 0.0832 0.0771 0.063 0.045
  ..$ integrals: num [1:4] 5.42e-05 2.13e-05 1.47e-04 7.49e-05
 $ mse  :List of 4
  ..$ raw   : NULL
  ..$ normed: NULL
  ..$ sm    : NULL
  ..$ smnorm: NULL
 - attr(*, "class")= chr "decon2"