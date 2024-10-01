# GLC v13 #####

# IMPORTANT: we dont test for jcampdx files because after calling `read_spectrum()`, the data is the same as for bruker, which is tested in `test-read_spectrum.R`. Also, the calculations in the old `MetaboDecon1D()` function are slightly different for jcampdx and bruker files (in the jcampdx case it calculates with n-1 instead of n) and our `compare_spectra_v13` function currently only accounts for bruker-type errors of MetaboDecon1D, but not jcampdx errors.

# test_that("GLC works for 1 bruker", {
#     data_path <- metabodecon_file("bruker/sim_subset")
#     x <- generate_lorentz_curves_sim(data_path, verbose = FALSE)
#     # new <- generate_lorentz_curves_sim(dp = "sim_01", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv
#     # old <- md1d(dp = "sim_01", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE)$rv
#     # r <- compare_spectra_v13(new, old, silent = TRUE)
#     # expect_true(sum(r %in% 0:1) >= 60 && sum(r %in% 2:3) == 0) # >=60 identical/equal && no diffs/errors
# })

# test_that("GLC works for bruker folder", {
#     data_path <- metabodecon_file("bruker/sim_subset")
#     x <- generate_lorentz_curves_sim(data_path, verbose = FALSE)
#     evobj <- md1d(dp = "sim_subset", ff = "bruker", nfit = 3, simple = TRUE, cache = FALSE, debug = FALSE)
#     y <- evobj$rv
#     expect_equal(names(x), names(y)) # chr [1:2] "sim_01" "sim_02"
#     expect_equal(names(x$sim_01)[1: length(y$sim_01)], names(y$sim_01))
#     new <- x$sim_01
#     old <- y$sim_01
#     expect_equal(new$number_of_files,            old$number_of_files           )
#     expect_equal(new$filename,                   old$filename                  )
#     expect_equal(new$x_values,                   old$x_values                  )
#     expect_equal(new$x_values_ppm,               old$x_values_ppm              )
#     expect_equal(new$y_values,                   old$y_values                  )
#     expect_equal(new$spectrum_superposition,     old$spectrum_superposition    )
#     expect_equal(new$mse_normed,                 old$mse_normed                )
#     expect_equal(new$index_peak_triplets_middle, old$index_peak_triplets_middle)
#     expect_equal(new$index_peak_triplets_left,   old$index_peak_triplets_left  )
#     expect_equal(new$index_peak_triplets_right,  old$index_peak_triplets_right )
#     expect_equal(new$peak_triplets_middle,       old$peak_triplets_middle      )
#     expect_equal(new$peak_triplets_left,         old$peak_triplets_left        )
#     expect_equal(new$peak_triplets_right,        old$peak_triplets_right       )
#     expect_equal(new$integrals,                  old$integrals                 )
#     expect_equal(new$signal_free_region,         old$signal_free_region        )
#     expect_equal(new$range_water_signal_ppm,     old$range_water_signal_ppm    )
#     expect_equal(new$A,                          old$A                         )
#     expect_equal(new$lambda,                     old$lambda                    )
#     expect_equal(new$x_0,                        old$x_0                       )
# })

# test_that("GLC works when no peaks are filtered out", {
#     x <- simulate_spectrum(ndp = 256, npks = 3)
#     expect_error(generate_lorentz_curves(x, sfr = c(Inf, -Inf), wshw = 0, smopts = c(0, 3), ask = FALSE))
#     decon <- generate_lorentz_curves(x, sfr = c(Inf, -Inf), wshw = 0, smopts = c(0, 3), ask = FALSE, force = TRUE)
#     expect_identical(length(decon), 31L)
# })
