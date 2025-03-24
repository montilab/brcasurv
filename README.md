<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# brcasurv

Helper functions to estimate the survival associations of user-provided gene signatures in TCGA & METABRIC.

## Installation:
```
devtools::install_github("montilab/brcasurv")
```

## Example Usage:

Given a named list of gene sets e.g. list(GS1 = c("gene1", "gene2", "gene3")) one can score each sample in TCGA/METABRIC with `gsva_data` and use `gsva_cox_fit` to fit coxph models.

```
gsva_data <- gsva_data(sigs_list = sigs_list, 
                       brca_data = "TCGA", # Or "METABRIC"
                       adjust_prolif = TRUE, 
                       adjust_inflam = TRUE)

gsva_cox_fits <- gsva_cox_fit(gsva_data,
                              brca_data = "TCGA", # Or "METABRIC"
                              adjust_age = TRUE,
                              adjust_prolif = TRUE,
                              adjust_inflam = TRUE)
```

## Plotting:

Use `survminer::ggadjustedcurves` to plot results from the cox models.

```
if (!require("survminer", quietly = TRUE)) {
  install.packages("survminer")
}

# Replace parameters with relevant data.
ggadjustedcurves(fit = {cox_fit},
                 data = {pData(gsva_data)},
                 method = "conditional",
                 variable = {signature_name},
                 xlab = "Days",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 title = "TCGA/METABRIC age-adjusted cox model"
)
```
