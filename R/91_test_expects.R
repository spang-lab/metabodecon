# Expect #####

#' @title Check if the size of each file in a directory is within a certain range
#' @description Check if the size of each file in a directory is within 90% to 110% of the expected size.
#' If a file size is not within this range, a message is printed and an error is thrown.
#' @param testdir A character string specifying the directory to check.
#' @param size_exp A named numeric vector where the names are filenames and the values are the expected file sizes.
#' @examples
#' \dontrun{
#' testdir <- tmpdir()
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
#' @description Tests if the structure of an object matches the expected string
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

testmatrix <- local({
    df <- expand.grid(dp = c("urine_1", "urine_2", "urine"), ff = c("bruker", "jcampdx"), nfit = c(1, 3, 10), simple = c(TRUE, FALSE), skip = TRUE, stringsAsFactors = FALSE)
    df$skip[df$ff == "bruker" | (df$ff == "jcampdx" & df$nfit == 3 & df$simple == TRUE)] <- FALSE
    x <- df$dp %in% c("urine_1", "urine_2") & df$ff == "jcampdx"
    df$dp[x] <- paste0(df$dp[x], ".dx")
    df
})

#' @description Generates a unique identifier for a test of `generate_lorentz_curves_v12` or `MetaboDecon1D`
#' @noRd
get_tid <- function(func, dp, ff, nfit, simple) {
    paste(func, dp, ff, nfit, simple, sep = "-")
}

#' @description Calls `func` for each row in `testmatrix` and caches the results
#' @param func Either "glc" or "MD1D"
#' @param overwrite Logical indicating whether to overwrite cached results if they already exist
#' @noRd
cache_func_results <- function(func = "glc", overwrite = FALSE) {
    df <- testmatrix
    cdir <- cachedir()
    tid <- get_tid(func, df$dp, df$ff, df$nfit, df$simple)
    rds <- file.path(cdir, paste0(tid, ".rds"))
    status <- ifelse(file.exists(rds), "cached", "todo")
    status[df$skip] <- "skip"
    callstr <- sprintf("%s(dp='%s', ff='%s', nfit=%d, simple=%s, overwrite=%s)", func, df$dp, df$ff, df$nfit, df$simple, overwrite)
    col <- ifelse(status == "cached", GREEN,  YELLOW)
    col[df$skip] <- BLUE
    df[, c("rds", "status", "col", "callstr")] <- list(rds, status, col, callstr)
    idxtodo <- which(status == "todo")
    idxdone <- which(status != "todo")
    process_row <- function(i) {
        row <- as.list(df[i, colnames(df)])
        cat2(row$callstr, " ", row$col, row$status, RESET, sep = "")
        x <- if (row$status == "todo") try(eval(parse(text = row$callstr))) else 0
        return(if (inherits(x, "try-error")) x else 0)
    }
    x <- lapply(idxdone, process_row) # dont spawn new processes for cached or skipped function calls
    y <- parallel::mclapply(idxtodo, process_row, mc.cores = ceiling(parallel::detectCores() / 2))
    return(unlist(c(x, y)))
}

cache_glc_results <- function(overwrite = FALSE) cache_func_results("glc", overwrite = FALSE)

cache_MD1D_results <- function(overwrite = FALSE) cache_func_results("MD1D", overwrite = FALSE)

glc_v13 <- function(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE, overwrite = FALSE, cout = TRUE, cplot = TRUE, cache = TRUE, debug = TRUE, nworkers = "auto") {
    tid <- get_tid("glc_v13", dp, ff, nfit, simple)
    if (dp %in% c("urine", "blood")) {
        answers <- c(SameParam = "y", AdjNo = "1", SFRok = "y", WSok = "y")
        inputs <- file.path(ff, dp)
    } else {
        answers <- c(SFRok = "y", WSok = "y")
        inputs <- file.path(ff, strsplit(dp, "_")[[1]][1], dp)
    }
    invisible(evalwith(
        testdir = tid,
        inputs = inputs,
        answers = answers,
        cache = cache,
        overwrite = overwrite,
        plot = if (cplot) "plots.pdf" else NULL,
        output = if (cout) "captured" else NULL,
        message = if (cout) "captured" else NULL,
        expr = generate_lorentz_curves(data_path = dp, file_format = ff, nfit = nfit, debug = debug, nworkers = nworkers)
    ))
}

glc_v12 <- function(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE, overwrite = FALSE, cout = TRUE, cplot = TRUE, cache = TRUE, debug = TRUE) {
    tid <- get_tid("glc", dp, ff, nfit, simple)
    inputs <- file.path(ff, "urine")
    if (dp != "urine") inputs <- file.path(ff, "urine", dp)
    answers <- c(SFRok = "y", WSok = "y")
    if (dp == "urine") answers <- c(SameParam = "y", AdjNo = "1", answers)
    invisible(evalwith(
        testdir = tid,
        inputs = inputs,
        answers = answers,
        cache = cache,
        overwrite = overwrite,
        plot = if (cplot) "plots.pdf" else NULL,
        output = if (cout) "captured" else NULL,
        message = if (cout) "captured" else NULL,
        expr = generate_lorentz_curves_v12(data_path = dp, file_format = ff, nfit = nfit, debug = debug)
    ))
}

MD1D <- function(dp = "urine_1", ff = "bruker", nfit = 3, simple = TRUE, overwrite = FALSE, cout = TRUE, cplot = TRUE, cache = TRUE, debug = TRUE) { # nolint: object_usage_linter.
    tid <- get_tid("MD1D", dp, ff, nfit, simple)
    if (dp %in% c("urine", "blood")) {
        # DATADIR/example_datasets
        #       {bruker|jcampdx} === ff (file format, input)
        #           {urine|blood} == dp (data path, input)
        # TESTDIR
        #       urine ============== fp (file path, MetaboDecon1D arg)
        inputs <- file.path(ff, dp)
        fp <- dp
        fn <- NA
    } else {
        # DATADIR/example_datasets
        #   {bruker|jcampdx} ====================== ff (file format, input)
        #       {urine|blood} ===================== pp (parent path, new variable)
        #           {blood_n|urine_n|urine_n.dx} == dp (data path, input)
        # TESTDIR
        #   . ===================================== fp (file path, MetaboDecon1D arg)
        #       {blood_n|urine_n|urine_n.dx} ====== fn (file name, MetaboDecon1D arg))
        pp <- strsplit(dp, "_")[[1]][1]
        inputs <- file.path(ff, pp, dp)
        fp <- "."
        fn <- dp
    }
    if (simple && !is.na(fn))  answers <- c(SFRok = "y", WSok = "y")
    if (simple && is.na(fn))   answers <- c(SameParam = "y", AdjNo = "1", SFRok = "y", WSok = "y")
    if (!simple && !is.na(fn)) answers <- c(SFRok = "n", Left = "11", Right = "-1", SFRok = "y", WSok = "asdf", WSok = "n", WSHW = "0.13", WSok = "y")
    if (!simple && is.na(fn))  answers <- c(SameParam = "n", answers, answers)
    if (ff == "bruker") answers <- c(ExpNo = "10", ProcNo = "10", answers)
    invisible(evalwith(
        testdir = tid,
        inputs = inputs,
        answers = answers,
        cache = cache,
        overwrite = overwrite,
        plot = if (cplot) "plots.pdf" else NULL,
        output = if (cout) "captured" else NULL,
        message = if (cout) "captured" else NULL,
        expr = MetaboDecon1D(filepath = fp, filename = fn, file_format = ff, number_iterations = nfit, debug = debug)
    ))
}
