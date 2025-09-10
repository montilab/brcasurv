library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

do_save <- FALSE
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas")
DATA_PATH <- file.path(Sys.getenv("CBM"), "SCANB/StringTie_prepDE_gene_count_data_unadjusted")
options(box.path=file.path(Sys.getenv("MLAB"), "personal/andrewdr/MLscripts"))
box::use(R/rm_low_rnaseq_counts)

# Need to download this Rdata from https://data.mendeley.com/datasets/yzxtxn4nmd/4
load(file.path(DATA_PATH, "SCANB.9142.matrixprepDEgenecount.Rdata"))
scanb_metadata <- readxl::read_excel(file.path(PATH, "data/revision/41523_2022_465_MOESM2_ESM.xlsx"), sheet = 1)
stopifnot(all(colnames(matrixprepDEgenecount) %in% scanb_metadata$GEX.assay))
# stopifnot(all(scanb_metadata$GEX.assay %in% colnames(matrixprepDEgenecount)))

rownames(scanb_metadata) <- scanb_metadata$GEX.assay
scanb_metadata_filtered <- scanb_metadata[colnames(matrixprepDEgenecount),]

# Creating the summarized experiment
scanb_se <- SummarizedExperiment(assays=list(counts=matrixprepDEgenecount),
                                 colData=S4Vectors::DataFrame(scanb_metadata_filtered))

# Filtering low count genes. Reduced from 19675 to 18426 genes.
scanb_se <- rm_low_rnaseq_counts$rm_low_rnaseq_counts(scanb_se, min_samples = 3)
saveRDS(scanb_se, file.path(DATA_PATH, "scanb_se.rds"))

# ----- DESeq NORMALIZATION -----
# Convert to DESeqDataSet
dds <- DESeqDataSet(scanb_se, design = ~ 1)  # design can be updated if needed

# Calculate size factors
dds <- estimateSizeFactors(dds)

# Obtain normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# ----- LOG2 TRANSFORMATION -----
log_norm_counts <- log2(norm_counts + 1)

# ----- UPDATE assay in SummarizedExperiment -----
# Replace counts assay directly (overwrites raw counts)
assay(scanb_se) <- log_norm_counts

# ----- UPDATE features for SummarizedExperiment ----
# Suppose 'se' is your SummarizedExperiment object
ensembl_ids <- rownames(scanb_se)

# If Ensembl IDs have version numbers (e.g., ENSG000001234.5), strip them
ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)

# Map Ensembl IDs to HGNC symbols using org.Hs.eg.db
hgnc_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add symbols as rowData or as a new column
rowData(scanb_se)$HGNC <- hgnc_symbols
rowData(scanb_se)$ENSEMBL <- rownames(scanb_se)

# Subset: keep entries with non-NA, non-duplicated symbols
keep <- !is.na(hgnc_symbols) & !duplicated(hgnc_symbols)

# 18426 features to 18129
scanb_se <- scanb_se[keep, ]

# Set rownames to HGNC symbols
rownames(scanb_se) <- rowData(scanb_se)$HGNC

scanb_data <- scanb_se
usethis::use_data(scanb_data, overwrite = TRUE)
