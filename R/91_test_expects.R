# Expect #####

#' @title Check if the size of each file in a directory is within a certain range
#' @description This function checks if the size of each file in a directory is within 90% to 110% of the expected size.
#' If a file size is not within this range, a message is printed and an error is thrown.
#' @param testdir A character string specifying the directory to check.
#' @param size_exp A named numeric vector where the names are filenames and the values are the expected file sizes.
#' @examples
#' \dontrun{
#' testdir <- tempdir()
#' file.create(file.path(testdir, "file1.txt"))
#' file.create(file.path(testdir, "file2.txt"))
#' size_exp <- c(file1.txt = 100, file2.txt = 200)
#' expect_file_size(testdir, size_exp)
#' }
#' @noRd
expect_file_size <- function(testdir, size_exp) {
    paths <- file.path(testdir, names(size_exp))
    size_obs <- file.info(paths)$size
    file_has_correct_size <- isTRUE(size_obs > size_exp * 0.9 & size_obs < size_exp * 1.1)
    lapply(seq_along(size_exp), function(i) {
        if (!isTRUE(file_has_correct_size[i])) {
            message(sprintf("Size of %s is %s which is not between %s and %s", paths[i], size_obs[i], size_exp[i] * 0.9, size_exp[i] * 1.1))
        }
    })
    testthat::expect_true(all(file_has_correct_size))
}

#' @title Expect Structure
#' @description Test if the structure of an object matches the expected string
#' @param obj The object to test
#' @param expected_str The expected structure of the object as a string. Can be obtained by calling `dput(capture.output(str(obj)))`.
#' @return A logical value indicating whether the structure of the object matches the expected string
#' @examples
#' expect_str(list(a = 1, b = 2), c("List of 2", " $ a: num 1", " $ b: num 2"))
#' @noRd
expect_str <- function(obj, expected_str) {
    testthat::expect_identical(capture.output(str(obj)), expected_str)
}

# Structures #####

str_urine_1_deconvoluted <- function(cf = 1, nf = 1, dx = FALSE, nested = TRUE, ni = 10) {
    fn <- if (dx) "urine_1.dx" else "urine_1"
    elemstr <- c(
        sprintf("number_of_files           : int %d", nf),
        sprintf('filename                  : chr "%s"', fn),
        sprintf("x_values                  : num [1:131072] 131 131 131 131 131 ..."),
        sprintf("x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ..."),
        sprintf("y_values                  : num [1:131072] 0.000831 0.000783 0.000743 0.000717 0.00065 ..."),
        sprintf("spectrum_superposition    : num [1, 1:131072] 3.51e-05 3.51e-05 3.51e-05 3.51e-05 3.52e-05 ..."),
        sprintf("mse_normed                : num 3.92e-11"),
        sprintf("index_peak_triplets_middle: num [1:1227] 36159 37149 37419 37435 38943 ..."),
        sprintf("index_peak_triplets_left  : num [1:1227] 36161 37160 37423 37438 38949 ..."),
        sprintf("index_peak_triplets_right : num [1:1227] 36156 37140 37415 37432 38938 ..."),
        sprintf("peak_triplets_middle      : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("peak_triplets_left        : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("peak_triplets_right       : num [1:1227] 9.28 9.13 9.09 9.08 8.85 ..."),
        sprintf("integrals                 : num [1, 1:1227] 0.000501 0.026496 0.000402 0.000375 0.008274 ..."),
        # sprintf("signal_free_region        : num [1:2] 109.1 21.9"),
        # sprintf("range_water_signal_ppm    : num 0.153"),
        sprintf("A                         : num [1:1227] -0.00016 -0.008436 -0.000128 -0.000119 -0.002634 ..."),
        sprintf("lambda                    : num [1:1227] -0.00775 -0.02188 -0.00675 -0.00562 -0.01343 ..."),
        sprintf("x_0                       : num [1:1227] 94.9 93.9 93.7 93.6 92.1 ...")
    )
    ne <- length(elemstr)
    plainstr <- c(
        paste0("List of ", ne),
        paste0(" $ ", elemstr)
    )
    nestedstr <- c(
        paste0("List of ", nf),
        paste0(" $ ", fn, ":List of ", ne),
        paste0("  ..$ ", elemstr)
    )
    if (nested) nestedstr else plainstr
}

str_urine_2_deconvoluted <- function(cf = 2, nf = 2, dx = FALSE, nested = TRUE, ni = 10) {
    fn <- if (dx) "urine_2.dx" else "urine_2"
    elemstr <- c(
        sprintf("number_of_files           : int %d", nf),
        sprintf('filename                  : chr "%s"', fn),
        sprintf("x_values                  : num [1:131072] 131 131 131 131 131 ..."),
        sprintf("x_values_ppm              : num [1:131072] 14.8 14.8 14.8 14.8 14.8 ..."),
        sprintf("y_values                  : num [1:131072] 0.00586 0.00578 0.00569 0.00557 0.00548 ..."),
        sprintf("spectrum_superposition    : num [1, 1:131072] 4.21e-05 4.21e-05 4.21e-05 4.21e-05 4.21e-05 ..."),
        sprintf("mse_normed                : num 2.86e-11"),
        sprintf("index_peak_triplets_middle: num [1:1393] 36290 37241 38346 38826 39025 ..."),
        sprintf("index_peak_triplets_left  : num [1:1393] 36297 37244 38349 38835 39028 ..."),
        sprintf("index_peak_triplets_right : num [1:1393] 36285 37234 38343 38823 39019 ..."),
        sprintf("peak_triplets_middle      : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ..."),
        sprintf("peak_triplets_left        : num [1:1393] 9.26 9.12 8.95 8.87 8.85 ..."),
        sprintf("peak_triplets_right       : num [1:1393] 9.26 9.12 8.95 8.88 8.85 ..."),
        sprintf("integrals                 : num [1, 1:1393] 0.00683 0.00504 0.00322 0.0174 0.00274 ..."),
        # sprintf("signal_free_region        : num [1:2] 109.1 21.9"),
        # sprintf("range_water_signal_ppm    : num 0.153"),
        sprintf("A                         : num [1:1393] -0.002176 -0.001604 -0.001025 -0.005541 -0.000872 ..."),
        sprintf("lambda                    : num [1:1393] -0.0189 -0.0168 -0.014 -0.0291 -0.0146 ..."),
        sprintf("x_0                       : num [1:1393] 94.8 93.8 92.7 92.2 92 ...")
    )
    ne <- length(elemstr)
    plainstr <- c(
        paste0("List of ", ne),
        paste0(" $ ", elemstr)
    )
    nestedstr <- c(
        paste0("List of ", nf),
        paste0(" $ ", fn, ":List of ", ne),
        paste0("  ..$ ", elemstr)
    )
    if (nested) nestedstr else plainstr
}

str_urine_deconvoluted <- function(nf = 2, dx = FALSE, nested = TRUE, ni = 10) {
    u1 <- str_urine_1_deconvoluted(nf = 2, dx = dx, nested = TRUE, ni = ni)
    u2 <- str_urine_2_deconvoluted(nf = 2, dx = dx, nested = TRUE, ni = ni)
    c("List of 2", u1[2:length(u1)], u2[2:length(u2)])
}

# MetaboDecon1D #####

cache_MetaboDecon1D_results <- function(overwrite = FALSE) {
    ns <- asNamespace("metabodecon")
    fsall <- ls(ns)
    fs <- fsall[grep("^MetaboDecon1D_", fsall)]
    for (f in fs) {
        callstr <- sprintf("%s(overwrite = %s)", f, overwrite)
        cat(callstr, "...", sep="")
        tryCatch(
            expr = {
                eval(parse(text = callstr))
                cat2(GREEN, "ok", RESET, sep = "")
            },
            error = function(e) {
                cat2(RED, e$message, RESET, sep = "")
            }
        )
    }
}

MetaboDecon1D_urine1_1010yy_ni1 <- function(cache = TRUE, overwrite = FALSE) {
    # Fastest possible way of running MetaboDecon1D
    x <- evalwith(
        testdir = "MetaboDecon1D_urine1_1010yy_ni1", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 1)
    )
}

MetaboDecon1D_dot_urine1_1010y1yy_ni3 <- function(cache = TRUE, overwrite = FALSE) {
    x <- evalwith(
        testdir = "MetaboDecon1D_dot_urine1_1010y1yy_ni3", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SameParam="y", AdjNo=1, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", file_format = "bruker")
    )
}

MetaboDecon1D_dot_urine1_1010y1yy_ni10_dbg <- function(cache = TRUE, overwrite = FALSE) {
    x <- evalwith(
        testdir = "MetaboDecon1D_dot_urine1_1010y1yy_ni10_dbg", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SameParam="y", AdjNo=1, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", file_format = "bruker", debug = TRUE)
    )
}

MetaboDecon1D_urine1_1010yy_ni1_run2 <- function(cache = TRUE, overwrite = FALSE) {
    # To test independence of seed. Should give same result as MetaboDecon1D_urine1_1010yy_ni1.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine1_1010yy_ni1_run2", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 1)
    )
}

MetaboDecon1D_urine1_1010yy_ni1_run3 <- function(cache = TRUE, overwrite = FALSE) {
    # To test independence of seed. Should give same result as MetaboDecon1D_urine1_1010yy_ni1.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine1_1010yy_ni1_run3", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 1)
    )
}

MetaboDecon1D_urine1_1010yy_ni1_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Test that debug doesnt change result. Should give same result as MetaboDecon1D_urine1_1010yy_ni1 except for additional debuglist element in the returned list.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine1_1010yy_ni1_dbg", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", # output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 1, debug = TRUE)
    )
}

MetaboDecon1D_urine1_1010yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of a single bruker file with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine1_1010yy_ni3_dbg", inputs = "bruker/urine/urine_1",
        answers = c(ExpNo=10, ProcNo=10, SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", filename = "urine_1", file_format = "bruker", number_iterations = 3, debug = TRUE)
    )
}

MetaboDecon1D_urine1dx_yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of a single jcampdx file with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine1dx_yy_ni3_dbg", inputs = "jcampdx/urine/urine_1.dx",
        answers = c(SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = ".", filename = "urine_1.dx", file_format = "jcampdx", number_iterations = 3, debug = TRUE)
    )
}

MetaboDecon1D_urine_1010y1yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple bruker files with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine_1010y1yy_ni3_dbg", inputs = "bruker/urine",
        answers = c(ExpNo = "10", ProcNo = "10", SameParam = "y", AdjNo="1", SFRok = "y", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = "urine", file_format = "bruker", number_iterations = 3, debug = TRUE)
    )
}

MetaboDecon1D_urinedx_y1yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple jcampdxs file with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "MetaboDecon1D_urinedx_y1yy_ni3_dbg", inputs = "jcampdx/urine",
        answers = c(SameParam = "y", AdjNo="1", SFRok = "y", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = "urine", file_format = "jcampdx", number_iterations = 3, debug = TRUE)
    )
}

MetaboDecon1D_urine_1010y1yy_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple bruker files with default args.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine_1010y1yy_ni3_dbg", inputs = "bruker/urine",
        answers = c(ExpNo = "10", ProcNo = "10", SameParam = "y", AdjNo="1", SFRok = "y", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = "urine", file_format = "bruker", debug = TRUE)
    )
}

MetaboDecon1D_urinedx_y1yy_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple jcampdx files with default args.
    x <- evalwith(
        testdir = "MetaboDecon1D_urinedx_y1yy_dbg", inputs = "jcampdx/urine",
        answers = c(SameParam = "y", AdjNo="1", SFRok = "y", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = "urine", file_format = "jcampdx", debug = TRUE)
    )
}

MetaboDecon1D_urine_1010nyyn11m1yXn013y_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple bruker files with complicated args.
    x <- evalwith(
        testdir = "MetaboDecon1D_urine_1010nyyn11m1yXn013y_ni3_dbg", inputs = "bruker/urine",
        answers = c(ExpNo = "10", ProcNo = "10", SameParam = "n", SFRok = "y", WSok = "y", SFRok = "n", Left = "11", Right = "-1", SFRok = "y", WSok = "asdf", WSok = "n", WSHW = "0.13", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = MetaboDecon1D(filepath = "urine", file_format = "bruker", number_iterations = 3, debug = TRUE)
    )
}

# Generate Lorentz Curves #####

glc_testid <- function(dp, ff, nfit) paste("glc", dp, ff, nfit, sep = "-")

glc <- function(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE, overwrite = FALSE, cout = TRUE, cplot = TRUE, cache = TRUE, debug = TRUE) {
    inputs <- file.path(ff, "urine")
    answers <- c(SFRok="y", WSok="y")
    if (dp == "urine") answers <- c(SameParam = "y", AdjNo="1", answers)
    if (dp != "urine") inputs <- file.path(ff, "urine", dp)
    if (dp %in% c("urine_1", "urine_2") && ff == "jcampdx") dp <- paste0(dp, ".dx")
    x <- evalwith(
        testdir = glc_testid(dp, ff, nfit),
        inputs = inputs,
        answers = answers,
        cache = cache,
        overwrite = overwrite,
        plot = if (cplot) "plots.pdf" else NULL,
        output = if (cout) "captured" else NULL,
        message = if (cout) "captured" else NULL,
        expr = generate_lorentz_curves_v12(data_path = dp, file_format = ff, nfit = nfit, debug = debug)
    )
}

cache_glc_results <- function(overwrite = FALSE) {
    df <- expand.grid(dp = c("urine_1", "urine_2", "urine"), ff = c("bruker", "jcampdx"), nfit = c(1, 3, 10), simple = c(TRUE, FALSE), skip = TRUE, stringsAsFactors = FALSE)
    df$skip[df$ff == "bruker" | (df$ff == "jcampdx" & df$nfit == 3 & df$simple == TRUE)] <- FALSE
    x <- df$dp %in% c("urine_1", "urine_2") & df$ff == "jcampdx"
    df$dp[x] <- paste0(df$dp[x], ".dx")
    cdir <- cachedir()
    process_row <- function(i) {
    row <- as.list(df[i, colnames(df)])
    rds <- file.path(cdir, paste0(glc_testid(row$dp, row$ff, row$nfit), ".rds"))
    status <- if (row$skip) "skip" else if (file.exists(rds)) "cached" else "..."
    callstr <- sprintf("glc(dp='%s', ff='%s', nfit=%d, simple=%s, overwrite=%s)", row$dp, row$ff, row$nfit, row$simple, overwrite)
    col <- if (status == "skip") BLUE else if (status == "cached") GREEN else YELLOW
    cat2(callstr, " ", col, status, RESET, sep = "")
    if (status != "skip") eval(parse(text = callstr))
    }
    mclapply(seq_len(nrow(df)), process_row, mc.cores = 8)
}

glc_urine2_yy_ni1_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of a single bruker file with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "glc_urine2_yy_ni1_dbg", inputs = "bruker/urine/urine_2",
        answers = c(SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = generate_lorentz_curves_v12(data_path = "urine_2", file_format = "bruker", nfit = 1, debug = TRUE)
    )
}

glc_urine1_yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of a single bruker file with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "glc_urine1_yy_ni3_dbg", inputs = "bruker/urine/urine_1",
        answers = c(SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = generate_lorentz_curves_v12(data_path = "urine_1", file_format = "bruker", nfit = 3, debug = TRUE)
    )
}

glc_urine1dx_yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of a single bruker file with only 3 iterations instead of 10 to speed things up.
    x <- evalwith(
        testdir = "glc_urine1dx_yy_ni3_dbg", inputs = "jcampdx/urine/urine_1.dx",
        answers = c(SFRok="y", WSok="y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", output = "captured", message = "captured",
        expr = {
            generate_lorentz_curves_v12(data_path = "urine_1.dx", file_format = "jcampdx", nfit = 3, debug = TRUE)
        }
    )
}

glc_urine_y1yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple bruker files with default args.
    x <- evalwith(
        testdir = "glc_urine_y1yy_ni3_dbg", inputs = "bruker/urine",
        answers = c(SameParam = "y", AdjNo="1", SFRok = "y", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", # output = "captured", message = "captured",
        expr = generate_lorentz_curves_v12(data_path = "urine", file_format = "bruker", nfit = 3, debug = TRUE)
    )
}

glc_urinedx_y1yy_ni3_dbg <- function(cache = TRUE, overwrite = FALSE) {
    # Deconvolution of multiple bruker files with default args.
    x <- evalwith(
        testdir = "glc_urinedx_y1yy_ni3_dbg", inputs = "jcampdx/urine",
        answers = c(SameParam = "y", AdjNo="1", SFRok = "y", WSok = "y"),
        cache = cache, overwrite = overwrite,
        plot = "plots.pdf", # output = "captured", message = "captured",
        expr = generate_lorentz_curves_v12(data_path = "urine", file_format = "jcampdx", nfit = 3, debug = TRUE)
    )
}
