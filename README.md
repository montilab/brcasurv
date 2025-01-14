# brcasurv

Helper functions to estimate the survival associations of user-provided gene signatures in TCGA & METABRIC.

## Example Usage:

Given a named list of gene sets e.g. list(GS1 = c("gene1", "gene2", "gene3")) one can score each sample in TCGA/METABRIC with `tcga_gsva` and use `tcga_cox_fit` to fit coxph models.

```
tcga_gsva <- tcga_gsva_data(sigs_list = sigs_list, 
                            adjust_prolif = TRUE, 
                            adjust_inflam = TRUE)

tcga_cox_fits <- tcga_cox_fit(tcga_gsva, 
                              five_year = FALSE,
                              adjust_age = TRUE,
                              adjust_prolif = TRUE,
                              adjust_inflam = TRUE)

tcga_cox_fits_5 <- tcga_cox_fit(tcga_gsva, 
                              five_year = TRUE,
                              adjust_age = TRUE,
                              adjust_prolif = TRUE,
                              adjust_inflam = TRUE)
```
