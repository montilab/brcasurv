library(GSVA)
library(SummarizedExperiment)
library(survival)
library(dplyr)
library(tibble)
library(assertthat)
library(S4Vectors)

#' @title GSVA scoring on TCGA + METABRIC
#' @description Processes gene set variation analysis (GSVA) data for breast cancer datasets.
#' @param sigs_list A list of character vectors representing gene signatures.
#' @param brca_data Character string specifying the breast cancer dataset to use. Options are "TCGA", "METABRIC", or "SCANB".
#' @param adjust_prolif Logical indicating whether to adjust for proliferation signature.
#' @param adjust_inflam Logical indicating whether to adjust for inflammation signature.
#' @return A GSVA object containing processed data.
#' @import GSVA SummarizedExperiment S4Vectors survival tibble dplyr assertthat
#' @export
gsva_data <- function(sigs_list,
                      brca_data = c("TCGA", "METABRIC", "SCANB"),
                      adjust_prolif = TRUE,
                      adjust_inflam = TRUE) {

  brca_data <- match.arg(brca_data, c("TCGA", "METABRIC", "SCANB"))

  assert_that(is.list(sigs_list), msg = "sigs_list must be a list")
  assert_that(!is.null(names(sigs_list)) && all(names(sigs_list) != ""), msg = "Elements of sigs_list must have names")
  assert_that(all(sapply(sigs_list, is.character)), msg = "Elements of sigs_list must be character vectors")

  # Loading signatures
  if (adjust_prolif) {
    stopifnot(exists("prolif_sig", where = "package:brcasurv"))
    data("prolif_sig", envir=environment())
    sigs_list <- c(sigs_list, prolif_sig)
  }
  if (adjust_inflam) {
    stopifnot(exists("inflam_sig", where = "package:brcasurv"))
    data("inflam_sig", envir=environment())
    sigs_list <- c(sigs_list, inflam_sig)
  }

  if (brca_data == "TCGA") {
    stopifnot(exists("tcga_data", where = "package:brcasurv"))
    data("tcga_data", envir=environment())

    gsva_param <- GSVA::gsvaParam(tcga_data, sigs_list, maxDiff = TRUE)
    gsva_data <- GSVA::gsva(gsva_param)

    # Survival Filtering
    na_filter <- !is.na(gsva_data$vital_status)
    missing_death_day_filter <- !((gsva_data$vital_status == "Dead") & is.na(gsva_data$days_to_death))
    data_filter <- na_filter & missing_death_day_filter

    gsva_data <- gsva_data[, data_filter]
    colData(gsva_data) <- colData(gsva_data) |>
      data.frame() |>
      tibble::rownames_to_column("ID") |>
      dplyr::mutate(time = if_else(!is.na(days_to_death), days_to_death, days_to_last_follow_up)) |>
      dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) |>
      dplyr::mutate(vital_status_1 = if_else(vital_status == "Alive", 1, 2)) |>
      dplyr::mutate(vital_status_5 = if_else(vital_status == "Dead" & (time_5 > 1825.0), 1, vital_status_1)) |>
      tibble::column_to_rownames("ID") |>
      S4Vectors::DataFrame()

  } else if (brca_data == "METABRIC") {
    stopifnot(exists("metabric_data", where = "package:brcasurv"))
    data("metabric_data", envir=environment())

    gsva_param <- GSVA::gsvaParam(metabric_data, sigs_list, maxDiff = TRUE)
    gsva_data <- GSVA::gsva(gsva_param)

    # Survival Filtering
    na_filter <- !is.na(gsva_data$OS_STATUS)
    missing_death_day_filter <- !((gsva_data$OS_STATUS == "DECEASED") & is.na(gsva_data$OS_MONTHS))
    data_filter <- na_filter & missing_death_day_filter
    gsva_data <- gsva_data[, data_filter]

    colData(gsva_data) <- colData(gsva_data) |>
      data.frame() |>
      dplyr::mutate(time = OS_MONTHS * 30.437) |>
      dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) |>
      dplyr::mutate(vital_status_1 = if_else(OS_STATUS == "LIVING", 1, 2)) |>
      dplyr::mutate(vital_status_5 = if_else(OS_STATUS == "DECEASED" & (time_5 > 1825.0), 1, vital_status_1)) |>
      S4Vectors::DataFrame()

  } else if (brca_data == "SCANB") {
    stopifnot(exists("scanb_data", where = "package:brcasurv"))
    data("scanb_data", envir=environment())

    gsva_param <- GSVA::gsvaParam(scanb_data, sigs_list, maxDiff = TRUE)
    gsva_data <- GSVA::gsva(gsva_param)

    # Survival Filtering
    na_filter <- !is.na(gsva_data$OS_event)
    missing_death_day_filter <- !((gsva_data$OS_event == 1) & is.na(gsva_data$OS_days))
    data_filter <- na_filter & missing_death_day_filter
    gsva_data <- gsva_data[, data_filter]

    colData(gsva_data) <- colData(gsva_data) |>
      data.frame() |>
      dplyr::mutate(time = OS_days) |>
      dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) |>
      dplyr::mutate(vital_status_1 = if_else(OS_event == 0, 1, 2)) |>
      dplyr::mutate(vital_status_5 = if_else(OS_event == 1 & (time_5 > 1825.0), 1, vital_status_1)) |>
      S4Vectors::DataFrame()
  }

  return(gsva_data)
}

#' @title Cox Proportional Hazards Model Fitting for GSVA Data with PH Check
#' @description Fits Cox proportional hazards models to GSVA data for breast cancer datasets,
#'   tests the proportional hazards assumption, and optionally filters out signatures that violate the PH assumption.
#' @param gsva_data A GSVA object containing processed data.
#' @param brca_data Character string specifying the breast cancer dataset used. Options are "TCGA", "METABRIC", or "SCANB".
#' @param subtype Character string specifying a PAM50 subtype to subset the data.
#' @param subset_subtype Logical indicating whether to subset data to a specified subtype.
#' @param five_year Logical indicating whether to use 5-year survival data.
#' @param adjust_age Logical indicating whether to adjust for age.
#' @param adjust_prolif Logical indicating whether to adjust for proliferation signature.
#' @param adjust_inflam Logical indicating whether to adjust for inflammation signature.
#' @param ph_filter Logical indicating whether to filter out signatures that violate the proportional hazards (PH) assumption based on \code{ph_alpha}. Default is \code{TRUE}.
#' @param ph_alpha Numeric threshold for the PH assumption test p-value. Signatures with PH p-value greater than this will be filtered out if \code{ph_filter} is \code{TRUE}. Default is 0.05.
#'
#' @return A list with components:
#' \describe{
#'   \item{cox_fits}{A list of fitted Cox proportional hazards models.}
#'   \item{ph_tests}{A list of proportional hazards test results for each model (objects returned by `cox.zph`).}
#' }
#' @import survival
#' @export
gsva_cox_fit <- function(gsva_data,
                         brca_data = c("TCGA", "METABRIC", "SCANB"),
                         subtype = c("LumA", "LumB", "Her2", "Basal", "Normal"),
                         subset_subtype = FALSE,
                         five_year = FALSE,
                         adjust_age = TRUE,
                         adjust_prolif = TRUE,
                         adjust_inflam = TRUE,
                         ph_filter = TRUE,
                         ph_alpha = 0.05) {

  brca_data <- match.arg(brca_data, c("TCGA", "METABRIC", "SCANB"))
  subtype <- match.arg(subtype, c("LumA", "LumB", "Her2", "Basal", "Normal"))

  all_sigs <- rownames(gsva_data)
  exp_sigs <- all_sigs[!(all_sigs %in% c("prolif", "inflam"))]

  if (five_year) {
    surv_response <- "Surv(as.numeric(time_5), vital_status_5)"
  } else {
    surv_response <- "Surv(as.numeric(time), vital_status_1)"
  }

  if (subset_subtype) {
    subtype_col_map <- list(
      TCGA = "subtype_selected",
      METABRIC = "Pam50_SUBTYPE",
      SCANB = "SSP.Subtype"
    )
    subtype_col <- subtype_col_map[[brca_data]]

    if (brca_data == "TCGA") {
      subtype <- paste0("BRCA.", subtype)
    }

    stopifnot("Subtype not present in dataset" = (subtype %in% gsva_data[[subtype_col]]))
    message("Subsetting for ", paste(subtype, "samples"))
    gsva_data <- gsva_data[,gsva_data[[subtype_col]] == subtype]
  }

  gsva_mat <- assays(gsva_data)[["es"]]

  if (adjust_prolif) {
    gsva_data[["prolif"]] <- t(gsva_mat["prolif",])
  }
  if (adjust_inflam) {
    gsva_data[["inflam"]] <- t(gsva_mat["inflam",])
  }

  age_var_map <- list(
    TCGA     = "age_at_index",
    METABRIC = "AGE_AT_DIAGNOSIS",
    SCANB    = "Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc.."
  )
  covariate_flags <- c(adjust_age, adjust_prolif, adjust_inflam)
  covariate_names <- c(age_var_map[[brca_data]], "prolif", "inflam")
  covariates <- covariate_names[covariate_flags]
  message("Adjusting for: ", paste(covariates, collapse = " "))

  cox_fits <- list()
  ph_tests <- list()

  for (sig in exp_sigs) {
    # Adding gsva scores to colData
    gsva_data[[sig]] <- t(gsva_mat[sig,])

    # Fitting Cox model
    fit <- coxph(reformulate(c(covariates, sig), response = surv_response), data = colData(gsva_data))
    # Check proportional hazards assumption
    ph_test <- cox.zph(fit)

    if(ph_filter) {
      ph_pvalue <- ph_test$table[sig, "p"]
      if(ph_pvalue <= ph_alpha) {
        next
      }
    }

    ph_tests[[sig]] <- ph_test
    cox_fits[[sig]] <- fit
  }

  return(list(cox_fits = cox_fits, ph_tests = ph_tests))
}

.onAttach <- function(libname, pkgname) {
  data_dir <- system.file("data", package = pkgname)

  base_url <- "https://js2.jetstream-cloud.org:8001/swift/v1/brcasurv/"
  files_to_download <- c("inflam_sig.rda", "prolif_sig.rda", "scanb_data.rda", "metabric_data.rda", "tcga_data.rda")

  # Check if all files exist
  all_files_exist <- all(sapply(files_to_download, function(f) file.exists(file.path(data_dir, f))))

  if (!all_files_exist) {
    packageStartupMessage("Downloading required data files.")
    tryCatch({
      for (filename in files_to_download) {
        file_path <- file.path(data_dir, filename)
        download_url <- paste0(base_url, filename)

        # Download each file
        packageStartupMessage(paste("Downloading", filename, "..."))
        utils::download.file(download_url, file_path, mode = "wb")
      }
      packageStartupMessage("Download complete.")
    }, error = function(e) {
      # Handle any errors that occur
      stop(paste("Error during data download:", e$message))
    })
  }
}
