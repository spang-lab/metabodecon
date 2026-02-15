#' @onRd
#' @title AKI Klassifikation
#' @description
#' Compares the performance of metabodecon against three other approaches
#' by performing an deconvolution, alignment and classification of the AKI
#' dataset.
#'
#' Requires execution of `update_aki_dataset()` first, which in turn
#' requires availability of Spang datasets under /data.
#'
#' @details
#' To estimate the end-to-end performance of metabodecon, this function
#' evaluates its deconvolution and alignment workflow against three preprocessing
#' strategies in a supervised classification setting: equidistant binning,
#' wavelet-based peak detection and alignment using speaq, and automated peak
#' picking and alignment as implemented in Bruker TopSpin. Equidistant binning
#' represents a widely adopted baseline in NMR metabolomics (Emwas et al.,
#' 2019), while wavelet-based peak detection improves peak localization and
#' noise robustness compared to fixed binning (Vu et al., 2011; Vu et al.,
#' 2019). Vendor-provided pipelines such as TopSpin reflect common analytical
#' practice and provide a real-world reference. Each preprocessing method yields
#' a feature matrix of aligned chemical shift positions with corresponding peak
#' integrals or binned intensities. To assess discriminative performance between
#' AKI patients and controls, the function trains logistic regression models with
#' an elastic-net penalty, which is well suited for high-dimensional, correlated
#' metabolomics data with limited sample size (Zou & Hastie, 2005; Kirpich et
#' al., 2018). Model hyperparameters (regularization strength and mixing
#' parameter) are tuned using nested cross-validation, and predictive performance
#' is quantified by balanced accuracy and area under the ROC curve (AUROC), with
#' feature scaling performed within each training fold to prevent information
#' leakage. This evaluation framework ensures an unbiased comparison of
#' preprocessing strategies and quantifies their end-to-end performance.
#'
classify_aki_patients <- function(data, labels) {

}

#' @noRd
#' @title Copies the AKI dataset from the Spang server to the package folder
#' @description
#' Copies the AKI dataset from the Spang server to the package folder. This
#' function must be run on a machine with access to the Spang datasets, e.g. r4,
#' r30 or r31 as of January 2026. It is not meant to be run by users of the
#' package, but only by developers during development. It is not exported and
#' not part of the package API.
#' @details
#' The following recommendations and/or hints are taken from Wolfram's email,
#' sent on Thursday, 12. February 2026 at 15:34 to tobias2.schmidt@ukr.de. The
#' exact email content is available at /data/aki/email_wolfram.txt on Spang
#' machines with access to /data. Ask Tobias Schmidt or Christian Kohler for
#' access, in case you need it.
#'
#' 1. The AKI dataset consists of 109 urine spectra from 110 study participants
#'    (sample 42 did not give urine).
#' 2. Samples were collected in Erlangen as part of the "AKI study".
#' 3. Samples include 73 healthy controls, 34 AKI patients and 3 cases having
#'    unclear diagnosis.
#' 4. Phenodata are stored in
#'    '/data/aki/Patientenvorhersagen_Erlangen_30Juli2012.xls'
#' 5. Use spectra folders named 'AKI_8_24_01-110' and pick '10/pdata10/1r'
#' 6. Sample 42 is missing (no probe)
#' 7. Sample 50 was measured twice; either measurement is fine (suffix _2)
#' 8. Samples 94 and 96 were remeasured after a bad first run; use _2 or _B
#' 9. Samples 47, 54, 83 should be excluded from classification due to unclear
#'    diagnosis
#' 10. AKI status is color coded in
#'     'Patientenvorhersagen_Erlangen_30Juli2012.xls' (yellow background
#'     indicates AKI, white background indicates healthy control)
#' 11. Always apply PQN normalization after bucketing/deconvolution to account
#'     for urine dilution
#'
update_aki_dataset <- function() {

    pheno_xls <- "/data/aki/Patientenvorhersagen_Erlangen_30Juli2012.xls"

    # Check Existence
    if (!file.exists(pheno_xls)) {
        msg <- paste(
            "Missing /data/aki/Patientenvorhersagen_Erlangen_30Juli2012.xls.",
            "Run this on a Spang server with /data mounted",
            "(r4, r30, r31; Jan 2026). Contact Tobias Schmidt or Christian Kohler."
        )
        stop(msg, call. = FALSE)
    }

    # Read in raw values
    if (!requireNamespace("xlsx", quietly = TRUE)) {
        stop("Missing package 'xlsx'. Install it before running this.",
             call. = FALSE)
    }
    if (!requireNamespace("rJava", quietly = TRUE)) {
        stop("Missing package 'rJava'. Install it before running this.",
             call. = FALSE)
    }
    raw <- xlsx::read.xlsx(pheno_xls, sheetIndex = 1, header = FALSE, stringsAsFactors = FALSE)
    raw <- as.data.frame(raw, stringsAsFactors = FALSE)

    # Now read in background color of first column as feature.
    # Color code 64 (default) indicates "healthy control".
    # Color code 51 (yellow) indicates "AKI patient".
    wb <- xlsx::loadWorkbook(pheno_xls)
    sheet <- xlsx::getSheets(wb)[[1]]
    rows <- xlsx::getRows(sheet)
    cells <- xlsx::getCells(rows, colIndex = 1)
    cell_rows <- sub("\\..*$", "", names(cells))
    cell_map <- stats::setNames(cells, cell_rows)
    row_idx <- seq_along(cell_rows)
    background <- vapply(row_idx, FUN.VALUE = integer(1), function(idx) {
        cell <- cell_map[[as.character(idx)]]
        if (is.null(cell)) return(NA_integer_)
        style <- xlsx::getCellStyle(cell)
        as.integer(style$getFillForegroundColor())
    })

    # Check and patch column names first
    names_observed <- as.character(raw[2, ]) # row 1 is empty
    names_expected <- c(
        "Pat. ID", "RIFLE", "AKIN 2d", "AKIN 3d", "RRT",
        "radiale SVM 24hrs \n106 urine samples \n5 runs",
        "radiale SVM 04hrs \n106 urine samples \n5 runs",
        "radiale SVM 24hrs \n105 plasma samples \n5 runs"
    )
    names_nice <- c(
        "PatID", "RIFLE", "AKIN_2d", "AKIN_3d", "RRT",
        "SVM_24h_urine", "SVM_4h_urine", "SVM_24h_plasma"
    )
    if (!isTRUE(all.equal(names_observed, names_expected))) {
        msg <- paste(
            "Unexpected column names in the AKI XLS.",
            "Open the file in Excel and update names_expected",
            "to match the current header row."
        )
        stop(msg, call. = FALSE)
    }
    pheno <- raw
    colnames(pheno) <- names_nice

    # Now check and filter rows. Make sure to append AKI status first, so vector
    # and matrix lengths match.
    pheno$HasAKI <- background == 51L
    pheno <- pheno[3:nrow(pheno), , drop = FALSE]
    pheno <- pheno[!is.na(pheno$PatID), , drop = FALSE]
    pheno$PatID <- as.integer(pheno$PatID)
    if (!isTRUE(all.equal(pheno$PatID, 1:110))) {
        msg <- paste(
            "Patient IDs are not 1..110 as expected.",
            "Check the XLS rows and the header row index",
            "and update the extraction accordingly."
        )
        stop(msg, call. = FALSE)
    }

    # As a last check, we compare the automatically extracted AKI status with a
    # vector created manually by opening
    # 'Patientenvorhersagen_Erlangen_30Juli2012.xls' in Excel and writing down
    # the AKI status for each patient, based on the background color.
    akiIDs_manual <- c(
      2, 4, 5, 6, 8, 15, 17, 18, 28, 29, 35, 36, 38, 45, 46, 48, 49,
      62, 65, 70, 71, 74, 75, 76, 80, 82, 89, 90, 92, 94, 98, 99,
      103, 106
    )
    akiIDs_auto <- pheno$PatID[pheno$HasAKI]
    if (!isTRUE(all.equal(akiIDs_auto, akiIDs_manual))) {
        msg <- paste(
            "AKI labels from background colors do not match",
            "the manually recorded IDs. Re-check the",
            "color palette mapping or the manual list."
        )
        stop(msg, call. = FALSE)
    }

    # Filter out patients with unclear diagnosis (47, 54, 83) and/or no urine
    # sample (42)
    pheno <- pheno[!pheno$PatID %in% c(42, 47, 54, 83), , drop = FALSE]

    # Final check, we should have 72 healthy controls and 34 AKI patients left.
    if (!isTRUE(sum(pheno$HasAKI) == 34L && sum(!pheno$HasAKI) == 72L)) {
        msg <- paste(
            "Unexpected class counts after filtering.",
            "Expected 72 controls and 34 AKI cases",
            "after removing 42, 47, 54, 83."
        )
        stop(msg, call. = FALSE)
    }

    # Now find over the actual spectra files
    src_dir <- "/data/aki/AKI-Studie_Urin_Erlangen"
    files <- dir(src_dir, pattern = "AKI_8_24_.*?_110.*")
    pattern <- "^AKI_8_24_([0-9]{2,3})_110[0-9]{3}(?:_([A-Za-z0-9]+))?$"
    parts <- regmatches(files, regexec(pattern, files, perl = TRUE))
    X <- data.frame(
      File = files,
      PatID = as.integer(sapply(parts, function(x) x[2])),
      Replicate = sapply(parts, function(x) x[3])
    )
    X <- X[order(X$PatID, X$Replicate), , drop = FALSE]
    rownames(X) <- NULL
    if (!all(X$Replicate %in% c("", "2", "B"))) {
        msg <- paste(
            "Unexpected Replicate suffix in spectra folders.",
            "Only '', '2', or 'B' are supported in",
            "AKI_8_24_* names."
        )
        stop(msg, call. = FALSE)
    }
    X$PatID <- as.integer(X$PatID)
    X$Replicate[is.na(X$Replicate) | X$Replicate == ""] <- 1
    X$Replicate[X$Replicate == "B"] <- 2
    X$Replicate <- as.integer(X$Replicate)

    # Remove patients with unclear diagnosis (47, 54, 83) and/or no urine sample (42)
    # and/or a second (potentially better) measurement (50, 94, 96)
    rm <- which(
        (X$PatID %in% c(42, 47, 54, 83)) |
        (X$PatID %in% c(50, 94, 96) & X$Replicate == 1)
    )
    X <- X[-rm, , drop = FALSE]

    # Final check, we should have the same PatIDs left as in the phenodata, and
    # each PatID should be represented once.
    if (!isTRUE(all.equal(sort(X$PatID), sort(pheno$PatID)))) {
        msg <- paste(
            "Phenodata IDs do not match available spectra IDs",
            "after filtering. Check file naming and",
            "the filtering rules for both tables."
        )
        stop(msg, call. = FALSE)
    }

    # Add filenames to phenodata
    phenoX <- merge(pheno, X, by = "PatID", all.x = TRUE, sort = FALSE)

    # Write phenodata to misc/example_datasets/bruker/aki/pheno.csv
    out_dir <- pkg_file("misc/example_datasets/bruker/aki")
    if (out_dir == "") {
        msg <- paste(
            "Output path is missing in the package tree.",
            "Run this from a dev checkout (devtools::load_all),",
            "not from an installed package."
        )
        stop(msg, call. = FALSE)
    }
    out_path <- file.path(out_dir, "pheno.csv")
    utils::write.csv(phenoX, out_path, row.names = FALSE)

    # Copy spectra files to misc/example_datasets/bruker/aki.
    req_files <- c("10/acqus", "10/pdata/10/procs", "10/pdata/10/1r")
    n <- nrow(phenoX)
    for (i in seq_len(n)) {
        logf("[%d/%d] Copying %s", i, n, phenoX$File[i])
        spl_name <- phenoX$File[i]
        src_dir_i <- file.path(src_dir, spl_name)
        dst_spl_i <- file.path(out_dir, spl_name)
        mkdirs(file.path(dst_spl_i, "10/pdata/10"))
        src_files <- file.path(src_dir_i, req_files)
        dst_files <- file.path(dst_spl_i, req_files)
        if (!all(file.exists(src_files))) {
            missing <- src_files[!file.exists(src_files)]
            msg <- paste("Missing files:\n", paste(missing, collapse = "\n"))
            stop(msg, call. = FALSE)
        }
        file.copy(src_files, dst_files, overwrite = TRUE)
        Sys.chmod(dst_files, mode = "0660", use_umask = FALSE)
    }

    invisible(NULL)
}

#' @noRd
#' @title Zips the AKI Dataset and prints its size
#' @description
#' Zips the AKI Dataset and prints the resulting size of the zip file
#' so we can easily update the `xds` variables in `data.R`. Requires
#' execution of `update_aki_dataset()` first.
zip_aki_dataset <- function() {
    out_dir <- pkg_file("misc/example_datasets/bruker/aki")
    if (out_dir == "") {
        msg <- paste(
            "Output path is missing in the package tree.",
            "Run this from a dev checkout (devtools::load_all),",
            "not from an installed package."
        )
        stop(msg, call. = FALSE)
    }
    pheno_csv <- file.path(out_dir, "pheno.csv")
    n_files <- length(dir(out_dir))
    if (!file.exists(pheno_csv) || n_files < 107) {
        msg <- paste(
            "AKI data not found under misc/example_datasets/bruker/aki.",
            "Run update_aki_dataset() to generate pheno.csv",
            "and copy spectra before zipping."
        )
        stop(msg, call. = FALSE)
    }

    zip_path <- file.path(dirname(out_dir), "aki.zip")
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(dirname(out_dir))
    if (file.exists(zip_path)) unlink(zip_path)
    utils::zip(zipfile = zip_path, files = "aki")

    size <- file.info(zip_path)$size
    message(sprintf("Wrote %s (%s)", zip_path, human_readable(size, "B")))
    zip_path
}