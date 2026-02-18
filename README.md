<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# brcasurv

Helper functions to estimate the survival associations of user-provided gene signatures in TCGA/METABRIC/SCANB.

## Installation:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(
  c("GSVA", "SummarizedExperiment", "S4Vectors"),
  dependencies = TRUE
)
remotes::install_github("montilab/brcasurv", dependencies = TRUE)
```

If installation fails due to missing transitive dependencies, install them explicitly and retry. For example:
```
install.packages("magick")
BiocManager::install(c("Biostrings", "rhdf5", "h5mread"))
```

## Getting Started Vignette

After installing the package:
```
browseVignettes("brcasurv")
```

The source vignette is available at `vignettes/getting-started.Rmd`.

## Example Usage:

Given a named list of gene sets e.g. list(GS1 = c("gene1", "gene2", "gene3")) one can score each sample in TCGA/METABRIC with `gsva_data` and use `gsva_cox_fit` to fit coxph models.

```
gsva_data <- gsva_data(sigs_list = sigs_list, 
                       brca_data = "TCGA", # Or "METABRIC/SCANB"
                       adjust_prolif = TRUE, 
                       adjust_inflam = TRUE)

gsva_cox_fits <- gsva_cox_fit(gsva_data,
                              brca_data = "TCGA", # Or "METABRIC/SCANB"
                              adjust_age = TRUE,
                              adjust_prolif = TRUE,
                              adjust_inflam = TRUE)
```

## Plotting:

Example of plotting code that uses `survminer::ggadjustedcurves` to plot results from the cox models. 
Since gsva scores of a geneset are continuous, we need to first make the geneset scores into a categorical variable and define a new cox model. 
```
sig_name <- "GS1"
gsva_data$sig <- as.numeric(SummarizedExperiment::assay(gsva_data, "es")[sig_name, ])
gsva_sig_median <- median(gsva_data$sig, na.rm = TRUE)
gsva_data$stat <- ifelse(gsva_data$sig < gsva_sig_median, "Low", "High")
cox_fit <- survival::coxph(
  survival::Surv(as.numeric(time_5), vital_status_5) ~ age_at_index + stat,
  data = SummarizedExperiment::colData(gsva_data)
)
```

Now we can plot the two adjusted survival curves.
```
if (!require("survminer", quietly = TRUE)) {
  install.packages("survminer")
}

# Replace parameters with relevant data.
survminer::ggadjustedcurves(fit = cox_fit,
                 data = SummarizedExperiment::colData(gsva_data),
                 method = "conditional",
                 variable = "stat",
                 xlab = "Days",
                 ylab = "Survival Probability",
                 pval = TRUE,
                 title = "TCGA/METABRIC/SCANB age-adjusted cox model"
)
```
