test_that("gsva_data validates sigs_list input", {
  expect_error(
    gsva_data(
      sigs_list = "not-a-list",
      brca_data = "TCGA",
      adjust_prolif = FALSE,
      adjust_inflam = FALSE
    ),
    "sigs_list must be a list"
  )

  expect_error(
    gsva_data(
      sigs_list = list(c("A", "B")),
      brca_data = "TCGA",
      adjust_prolif = FALSE,
      adjust_inflam = FALSE
    ),
    "must have names"
  )

  expect_error(
    gsva_data(
      sigs_list = list(sig1 = 1:3),
      brca_data = "TCGA",
      adjust_prolif = FALSE,
      adjust_inflam = FALSE
    ),
    "must be character vectors"
  )
})

test_that("gsva_data validates brca_data values", {
  expect_error(
    gsva_data(
      sigs_list = list(sig1 = c("A", "B")),
      brca_data = "INVALID",
      adjust_prolif = FALSE,
      adjust_inflam = FALSE
    ),
    "should be one of"
  )
})
