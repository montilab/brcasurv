#' TCGA BRCA Data
#'
#' ExpressionSet Object containing 37335 features (HGNC) and 1222 samples. Counts are filtered to remove
#' low-expressed genes, DESeq2normalized to normalize for library size, and log transformed (log2(x+1)).
"tcga_data"

#' METABRIC BRCA Data
#'
#' ExpressionSet Object containing 24368 features (HGNC) and 1980 samples. Counts are log transformed
"metabric_data"

#' SCANB Data
#'
#' SummarizedExperiment Object containing 18129 features (HGNC) and 9142 samples. Counts are filtered to remove
#' low-expressed genes, DESeq2normalized to normalize for library size, and log transformed (log2(x+1)).
"scanb_data"

#' Prolif Sig
#'
#' Proliferation signature from Venet et. al 2011 Most random gene expression signatures are
#' significantly associated with breast cancer outcome.
"prolif_sig"

#' Inflam Sig
#'
#' Inflammation signature from Winslow et. al 2015 Prognostic stromal gene signatures in breast cancer.
"inflam_sig"
