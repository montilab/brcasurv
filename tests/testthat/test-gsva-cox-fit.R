make_mock_gsva_data <- function(n = 120L) {
  set.seed(1)

  es_mat <- matrix(rnorm(4L * n), nrow = 4L)
  rownames(es_mat) <- c("sigA", "sigB", "prolif", "inflam")
  colnames(es_mat) <- paste0("sample", seq_len(n))

  # Ensure events and censoring are both present for Cox fitting.
  status <- rep(c(1, 2), length.out = n)

  col_data <- S4Vectors::DataFrame(
    age_at_index = rnorm(n, mean = 55, sd = 10),
    AGE_AT_DIAGNOSIS = rnorm(n, mean = 55, sd = 10),
    `Age..5.year.range..e.g...35.31.35...40.36.40...45.41.45..etc..` = rnorm(n, mean = 55, sd = 10),
    time = rexp(n, rate = 1 / 500) + 30,
    time_5 = rexp(n, rate = 1 / 500) + 30,
    vital_status_1 = status,
    vital_status_5 = status,
    subtype_selected = rep(c("BRCA.LumA", "BRCA.Basal"), length.out = n),
    Pam50_SUBTYPE = rep(c("LumA", "Basal"), length.out = n),
    SSP.Subtype = rep(c("LumA", "Basal"), length.out = n),
    row.names = colnames(es_mat)
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(es = es_mat),
    colData = col_data
  )
}

test_that("gsva_cox_fit returns expected structure", {
  gsva_obj <- make_mock_gsva_data()

  fit_res <- gsva_cox_fit(
    gsva_data = gsva_obj,
    brca_data = "TCGA",
    adjust_prolif = FALSE,
    adjust_inflam = FALSE,
    ph_filter = FALSE
  )

  expect_type(fit_res, "list")
  expect_named(fit_res, c("cox_fits", "ph_tests"))
  expect_length(fit_res$cox_fits, 2)
  expect_length(fit_res$ph_tests, 2)
  expect_s3_class(fit_res$cox_fits[[1]], "coxph")
  expect_s3_class(fit_res$ph_tests[[1]], "cox.zph")
})

test_that("gsva_cox_fit PH filtering can drop all signatures", {
  gsva_obj <- make_mock_gsva_data()

  fit_res <- gsva_cox_fit(
    gsva_data = gsva_obj,
    brca_data = "TCGA",
    ph_filter = TRUE,
    ph_alpha = 1
  )

  expect_length(fit_res$cox_fits, 0)
  expect_length(fit_res$ph_tests, 0)
})

test_that("gsva_cox_fit subtype subsetting works for TCGA labels", {
  gsva_obj <- make_mock_gsva_data()

  expect_no_error(
    gsva_cox_fit(
      gsva_data = gsva_obj,
      brca_data = "TCGA",
      subset_subtype = TRUE,
      subtype = "LumA",
      adjust_prolif = FALSE,
      adjust_inflam = FALSE,
      ph_filter = FALSE
    )
  )

  expect_error(
    gsva_cox_fit(
      gsva_data = gsva_obj,
      brca_data = "TCGA",
      subset_subtype = TRUE,
      subtype = "Normal",
      adjust_prolif = FALSE,
      adjust_inflam = FALSE,
      ph_filter = FALSE
    ),
    "Subtype not present"
  )
})
