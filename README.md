<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# brcasurv

Helper functions to estimate the survival associations of user-provided gene signatures in TCGA & METABRIC.

## Example Usage:

Given a named list of gene sets e.g. list(GS1 = c("gene1", "gene2", "gene3")) one can score each sample in TCGA/METABRIC with `gsva_data` and use `gsva_cox_fit` to fit coxph models.

```
gsva_data <- gsva_data(sigs_list = sigs_list, 
                       brca_data = "TCGA", # Or "METABRIC"
                       adjust_prolif = TRUE, 
                       adjust_inflam = TRUE)

gsva_cox_fits <- gsva_cox_fit(gsva_data, 
                              adjust_age = TRUE,
                              adjust_prolif = TRUE,
                              adjust_inflam = TRUE)
```
