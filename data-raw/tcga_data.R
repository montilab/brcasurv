require(Biobase)
require(DESeq2)
require(edgeR)

eset.symbols <- function(eset) {
  eset <- eset[fData(eset)$hgnc_symbol != "" & !duplicated(fData(eset)$hgnc_symbol) & !is.na(fData(eset)$hgnc_symbol),]
  rownames(eset) <- fData(eset)$hgnc_symbol
  rownames(fData(eset)) <- fData(eset)$hgnc_symbol
  return(eset)
}

eset.filter <- function(id.data) {
  indir <- "/restricted/projectnb/montilab-p/personal/anthony/gdc-data/esets_normalized"
  outdir <- "/restricted/projectnb/montilab-p/personal/anthony/gdc-data/esets_filtered"

  eset.edger <- readRDS(file.path(indir, paste0("TCGA-", id.data, "_2020-03-22_edgeR_log_eset.rds")))
  saveRDS(eset.symbols(eset.edger), file.path(outdir, paste0("TCGA-", id.data, "_2020-03-22_edgeR_log_filtered_eset.rds")))

  eset.deseq <- readRDS(file.path(indir, paste0("TCGA-", id.data, "_2020-03-22_DESeq2_log_eset.rds")))
  saveRDS(eset.symbols(eset.deseq), file.path(outdir, paste0("TCGA-", id.data, "_2020-03-22_DESeq2_log_filtered_eset.rds")))
}

eset.normalize <- function(id.data) {

  cat("Normalizing: ", id.data, "\n")

  # Using eset as a varibale name for convenience, not because this function is specific to TCGA-eset
  indir <- "/restricted/projectnb/montilab-p/personal/anthony/gdc-data/esets_annotated"
  outdir <- "/restricted/projectnb/montilab-p/personal/anthony/gdc-data/esets_normalized"
  eset <- readRDS(file=file.path(indir, paste0("TCGA-", id.data, "_2020-03-22_annotated_eset.rds")))
  print(dim(eset))

  # edgeR Normalization
  eset.tmm <- eset.tmm.log <- eset
  dge.tmm <- calcNormFactors(DGEList(eset), method = "TMM") # Create edgeR specific object
  exprs(eset.tmm) <- cpm(dge.tmm) # Calculate Scaling Factors
  exprs(eset.tmm.log) <- log(exprs(eset.tmm) + 1, 2) # Calculate log2 counts per million

  saveRDS(eset.tmm.log, file.path(outdir, paste0("TCGA-", id.data, "_2020-03-22_edgeR_log_eset.rds")))

  # DESeq2 Normalization
  eset.rle <- eset.rle.log <- eset
  dds <- DESeqDataSetFromMatrix(exprs(eset), pData(eset), formula( ~ 1)) # Create DEseq2 specific object
  dds.rle <- estimateSizeFactors(dds) # Estimate scaling factors
  exprs(eset.rle) <- counts(dds.rle, normalized=TRUE) # Get counts
  exprs(eset.rle.log) <- log(exprs(eset.rle) + 1, 2) # Calculate log2 counts

  saveRDS(eset.rle.log, file.path(outdir, paste0("TCGA-", id.data, "_2020-03-22_DESeq2_log_eset.rds")))
}

if (F) {
  eset.normalize("BRCA")
  eset.filter("BRCA")
}

usethis::use_data(tcga_data, overwrite = TRUE)
