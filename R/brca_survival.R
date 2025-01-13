library(GSVA)
library(Biobase)
library(survival)

tcga_gsva_data <- function(sigs_list, 
                           adjust_prolif = TRUE,
                           adjust_inflam = TRUE) {
  
  # Loading Data
  tcga_data <- load(file.path("data", "tcga_data.rda"))
  
  if (adjust_prolif) {
    prolif_sig <- load(file.path("data", "prolif_sig.rda"))
    sigs_list <- c(sigs_list, prolif_sig)
  }
  if (adjust_inflam) {
    inflam_sig <- load(file.path("data", "inflam_sig.rda"))
    sigs_list <- c(sigs_list, inflam_sig)
  }
  
  # Running GSVA
  tcga_gsva <- GSVA::gsva(tcga_data, sigs_list, mx.diff=TRUE, verbose=FALSE)
  
  # Survival Filtering 
  na_filter <- !is.na(tcga_gsva$vital_status)
  missing_death_day_filter <- !((tcga_gsva$vital_status == "Dead") & is.na(tcga_gsva$days_to_death))
  subtype_filter <- !is.na(tcga_gsva$subtype_selected)
  data_filter <- na_filter & missing_death_day_filter & subtype_filter
  
  tcga_gsva <- tcga_gsva[, data_filter]
  pData(tcga_gsva) <- pData(tcga_gsva) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::mutate(time = if_else(!is.na(days_to_death), days_to_death, days_to_last_follow_up)) %>%
    dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) %>%
    dplyr::mutate(vital_status_1 = if_else(vital_status == "Alive", 1, 2)) %>%
    dplyr::mutate(vital_status_5 = if_else(vital_status == "Dead" & (time_5 > 1825.0), 1, vital_status_1)) %>%
    tibble::column_to_rownames("ID")
  
  return(tcga_gsva)
}

tcga_cox_fit <- function(tcga_gsva,
                         five_year = FALSE,
                         adjust_age = TRUE,
                         adjust_prolif = TRUE,
                         adjust_inflam = TRUE) {
  
  all_sigs <- rownames(tcga_gsva) 
  all_sigs <- all_sigs[!(all_sigs %in% c("prolif", "inflam"))]
  
  if(five_year) {
    surv_response <- "Surv(as.numeric(time_5), vital_status_5)"
  } else {
    surv_response <- "Surv(as.numeric(time), vital_status_1)"
  }
  
  covariate_flags <- c(adjust_age, adjust_prolif, adjust_inflam)
  covariates <- c("age_at_index", "prolif", "inflam")[covariate_flags]
  print(paste(c("Adjusting for:", covariates), collapse = " "))
  
  cox_fits <- list()
  for (sig in all_sigs) {
    # Adding gsva scores to pdata
    tcga_gsva[[sig]] <- t(exprs(tcga_gsva[sig,]))
    
    # Fitting Cox model
    cox_fits[[sig]] <- coxph(reformulate(c(covariates, sig), response = surv_response), data = pData(tcga_gsva))
  }
  
  return(cox_fits)
}

metabric_gsva_data <- function(sigs_list, 
                               adjust_prolif = TRUE,
                               adjust_inflam = TRUE) {
  
  # Loading Data
  metabric_data <- load(file.path("data", "metabric_data.rda"))
  
  if (adjust_prolif) {
    prolif_sig <- load(file.path("data", "prolif_sig.rda"))
    sigs_list <- c(sigs_list, prolif_sig)
  }
  if (adjust_inflam) {
    inflam_sig <- load(file.path("data", "inflam_sig.rda"))
    sigs_list <- c(sigs_list, inflam_sig)
  }
  
  # Running GSVA
  metabric_gsva <- GSVA::gsva(metabric_data, sigs_list, mx.diff=TRUE, verbose=FALSE)

  pData(metabric_gsva) <- pData(metabric_gsva) %>%
    dplyr::mutate(time = OS_MONTHS * 30.437) %>%
    dplyr::mutate(time_5 = if_else(as.numeric(time) < 1825.0, as.numeric(time), 1826.0)) %>%
    dplyr::mutate(vital_status_1 = if_else(OS_STATUS == "LIVING", 1, 2)) %>%
    dplyr::mutate(vital_status_5 = if_else(OS_STATUS == "DECEASED" & (time_5 > 1825.0), 1, vital_status_1))
  
  return(metabric_gsva)
}

metabric_cox_fit <- function(metabric_gsva,
                             five_year = FALSE,
                             adjust_age = TRUE,
                             adjust_prolif = TRUE,
                             adjust_inflam = TRUE) {
  
  all_sigs <- rownames(metabric_gsva) 
  all_sigs <- all_sigs[!(all_sigs %in% c("prolif", "inflam"))]
  
  if(five_year) {
    surv_response <- "Surv(as.numeric(time_5), vital_status_5)"
  } else {
    surv_response <- "Surv(as.numeric(time), vital_status_1)"
  }
  
  covariate_flags <- c(adjust_age, adjust_prolif, adjust_inflam)
  covariates <- c("AGE_AT_DIAGNOSIS", "prolif", "inflam")[covariate_flags]
  print(paste(c("Adjusting for:", covariates), collapse = " "))
  
  cox_fits <- list()
  for (sig in all_sigs) {
    # Adding gsva scores to pdata
    metabric_gsva[[sig]] <- t(exprs(metabric_gsva[sig,]))
    
    # Fitting Cox model
    cox_fits[[sig]] <- coxph(reformulate(c(covariates, sig), response = surv_response), data = pData(metabric_gsva))
  }
  
  return(cox_fits)
}