library(GSVA)
library(Biobase)
library(survival)

#' @title GSVA scoring on TCGA + METABRIC
#' @description Processes gene set variation analysis (GSVA) data for breast cancer datasets.
#' @param sigs_list A list of character vectors representing gene signatures.
#' @param brca_data Character string specifying the breast cancer dataset to use. Options are "TCGA" or "METABRIC".
#' @param adjust_prolif Logical indicating whether to adjust for proliferation signature.
#' @param adjust_inflam Logical indicating whether to adjust for inflammation signature.
#' @return A GSVA object containing processed data.
#' @import GSVA Biobase survival here
#' @export
gsva_data <- function(sigs_list,
                      brca_data = c("TCGA", "METABRIC"),
                      adjust_prolif = TRUE,
                      adjust_inflam = TRUE) {

  brca_data <- match.arg(brca_data, c("TCGA", "METABRIC"))

  # Assertion statements
  stopifnot("sigs_list must be a list" = is.list(sigs_list))
  stopifnot("elements of sigs_list must be char vectors" = is.character(sigs_list[[1]]))

  # Loading signatures
  if (adjust_prolif) {
    data("prolif_sig")
    sigs_list <- c(sigs_list, prolif_sig)
  }
  if (adjust_inflam) {
    data("inflam_sig")
    sigs_list <- c(sigs_list, inflam_sig)
  }

  if (brca_data == "TCGA") {
    data("tcga_data")
    gsva_data <- GSVA::gsva(tcga_data, sigs_list, mx.diff=TRUE, verbose=FALSE)

    # Survival Filtering
    na_filter <- !is.na(gsva_data$vital_status)
    missing_death_day_filter <- !((gsva_data$vital_status == "Dead") & is.na(gsva_data$days_to_death))
    subtype_filter <- !is.na(gsva_data$subtype_selected)
    data_filter <- na_filter & missing_death_day_filter & subtype_filter

    gsva_data <- gsva_data[, data_filter]
    pData(gsva_data) <- pData(gsva_data) %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::mutate(time = if_else(!is.na(days_to_death), days_to_death, days_to_last_follow_up)) %>%
      dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) %>%
      dplyr::mutate(vital_status_1 = if_else(vital_status == "Alive", 1, 2)) %>%
      dplyr::mutate(vital_status_5 = if_else(vital_status == "Dead" & (time_5 > 1825.0), 1, vital_status_1)) %>%
      tibble::column_to_rownames("ID")

  } else if (brca_data == "METABRIC") {
    data("metabric_data")
    gsva_data <- GSVA::gsva(metabric_data, sigs_list, mx.diff=TRUE, verbose=FALSE)

    pData(gsva_data) <- pData(gsva_data) %>%
      dplyr::mutate(time = OS_MONTHS * 30.437) %>%
      dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) %>%
      dplyr::mutate(vital_status_1 = if_else(OS_STATUS == "LIVING", 1, 2)) %>%
      dplyr::mutate(vital_status_5 = if_else(OS_STATUS == "DECEASED" & (time_5 > 1825.0), 1, vital_status_1))
  }

  return(gsva_data)
}

#' @title Cox Proportional Hazards Model Fitting for GSVA Data
#' @description Fits Cox proportional hazards models to GSVA data for breast cancer datasets.
#' @param gsva_data A GSVA object containing processed data.
#' @param brca_data Character string specifying the breast cancer dataset used. Options are "TCGA" or "METABRIC".
#' @param five_year Logical indicating whether to use 5-year survival data.
#' @param adjust_age Logical indicating whether to adjust for age.
#' @param adjust_prolif Logical indicating whether to adjust for proliferation signature.
#' @param adjust_inflam Logical indicating whether to adjust for inflammation signature.
#' @return A list of fitted Cox proportional hazards models.
#' @import survival
#' @export
gsva_cox_fit <- function(gsva_data,
                         brca_data = c("TCGA", "METABRIC"),
                         five_year = FALSE,
                         adjust_age = TRUE,
                         adjust_prolif = TRUE,
                         adjust_inflam = TRUE) {

  brca_data <- match.arg(brca_data, c("TCGA", "METABRIC"))

  all_sigs <- rownames(gsva_data)
  exp_sigs <- all_sigs[!(all_sigs %in% c("prolif", "inflam"))]

  if(five_year) {
    surv_response <- "Surv(as.numeric(time_5), vital_status_5)"
  } else {
    surv_response <- "Surv(as.numeric(time), vital_status_1)"
  }

  if (adjust_prolif) {
    gsva_data[["prolif"]] <- t(exprs(gsva_data["prolif",]))
  }
  if (adjust_inflam) {
    gsva_data[["inflam"]] <- t(exprs(gsva_data["inflam",]))
  }

  if(brca_data == "TCGA") {
    covariate_flags <- c(adjust_age, adjust_prolif, adjust_inflam)
    covariates <- c("age_at_index", "prolif", "inflam")[covariate_flags]
    print(paste(c("Adjusting for:", covariates), collapse = " "))
  } else if (brca_data == "METABRIC") {
    covariate_flags <- c(adjust_age, adjust_prolif, adjust_inflam)
    covariates <- c("AGE_AT_DIAGNOSIS", "prolif", "inflam")[covariate_flags]
    print(paste(c("Adjusting for:", covariates), collapse = " "))
  }

  cox_fits <- list()
  for (sig in exp_sigs) {
    # Adding gsva scores to pdata
    gsva_data[[sig]] <- t(exprs(gsva_data[sig,]))

    # Fitting Cox model
    cox_fits[[sig]] <- coxph(reformulate(c(covariates, sig), response = surv_response), data = pData(gsva_data))
  }

  return(cox_fits)
}

.onLoad <- function(libname, pkgname) {
  data_dir <- system.file("data", package = pkgname)
  tcga_path <- file.path(data_dir, "tcga_data.rda")
  metabric_path <- file.path(data_dir, "metabric_data.rda")
  
  if(!file.exists(tcga_path) | !file.exists(metabric_path)) {
    message("Downloading TCGA/METABRIC data.")
    tryCatch({
      # Download zip file from dropbox
      temp_zip <- tempfile(fileext = ".zip")
      utils::download.file("https://www.dropbox.com/scl/fo/ncrbv8i7g1fs5xpofeoa8/AGdjjxSgSiG6tDuRb6RZimE?rlkey=a20j6knqys4npqev1fcul7kas&e=5&dl=1",
                           temp_zip,
                           mode = "wb")
      # Unzip contents into package directory
      utils::unzip(temp_zip, exdir = data_dir)
      # Delete zip file
      unlink(temp_zip)
      message("Download complete.")
    }, error = function(e) {
      # Handle any errors that occur
      stop(paste("Error during download or unzip:", e$message))
      
    }, finally = {
      # Ensure the temporary file is deleted even if there's an error
      if (file.exists(temp_zip)) {
        unlink(temp_zip)
      }
    })
  }
  
}