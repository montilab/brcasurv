library(Biobase)

# This script takes a look the some raw/processed METABRIC BRCA Gene Expression data, performs sanity checks, and assembles an ESet pertaining to the CBioPortal data for the same
# For more information on the data formats from CBioPortal, visit https://github.com/cBioPortal/cbioportal/blob/master/docs/File-Formats.md
# Source directories
# Amy's raw data directory
src_dir1 <- "/restricted/projectnb/montilab-p/CBMrepositoryData/CURTIS/raw_data/"
# CBioPortal data directory
src_dir2 <- "/restricted/projectnb/montilab-p/CBMrepositoryData/CURTIS/CBioPortal_data/"

# Read in gene expression file
# Discovery matrix
discov.mat <- read.delim(paste(src_dir1,"discovery_ExpressionMatrix.txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
dim(discov.mat)
discov.samples <- gsub(".","-",colnames(discov.mat),fixed = TRUE)

# Validation matrix
valid.mat <- read.table(paste(src_dir1,"validation_ExpressionMatrix.txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
dim(valid.mat)
valid.samples <- gsub(".","-",colnames(valid.mat),fixed=TRUE)

# CBioportal data matrix (this has almost both combined, looks like)
cbio.d <- read.delim(paste(src_dir2,"data_expression.txt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE)
# Note that this matrix has Hugo gene symbol and Entrez gene id as the first two fields
dim(cbio.d)

# Matrix of counts
cbio.mat <- as.matrix(cbio.d[,3:ncol(cbio.d)])
# Change column names/sample IDs to shorter form (i.e. leaving out the leading "BRCA.METABRIC.S1.", which was found in every single sample ID), replacing some characters for consistency
colnames(cbio.mat) <- gsub(".","-",substr(colnames(cbio.mat),18,24),fixed=TRUE)
# Making sure that there are no redundant (more than 1) gene symbol in the feature data
length(unique(cbio.d$Hugo_Symbol))==nrow(cbio.mat)
# Make these the rownames
rownames(cbio.mat) <- cbio.d$Hugo_Symbol

# Write to a txt file the entire list of sample ID's (in standard format) that are covered by CBioPortal
# Making sure all column names are fixed character size (to use substr)
table(sapply(colnames(cbio.mat),nchar))

# Fetch the actual sample names from the CBioPortal naming convention
cbio.samples <- gsub(".","-",substr(colnames(cbio.mat),18,24),fixed=TRUE)
head(cbio.samples)
write.table(data.frame("Sample_ID"=cbio.samples),paste(src_dir,"GE_sample_IDs.txt",sep=""),quote=FALSE,row.names=FALSE)

# Look at overlap in samples (to see what we're missing if we use CBioPortal)
# All 997 samples are covered in CBioPortal from the Discovery dataset
length(intersect(discov.samples,cbio.samples))
# 12/995 samples are not covered in CBioPortal from the Validation datset
length(intersect(valid.samples,cbio.samples))


# Confirming that the data is already log-transormed
# First two fields are gene symbol and Entrez gene ID, respectively
summary(cbio.mat)[,1]

# Load clinical data (CBioPortal)
clin.mat <- read.delim(paste(src_dir2,"data_clinical.txt",sep=""),sep="\t",header=TRUE,stringsAsFactors=FALSE,comment.char = "#")

# Showing that each patient has only one sample
all.equal(clin.mat$PATIENT_ID,substr(clin.mat$SAMPLE_ID,18,24))
length(unique(clin.mat$PATIENT_ID))

# Make sure that the clinical data Sample IDs are all covered in the expression matrix
length(intersect(clin.mat$PATIENT_ID,colnames(cbio.mat))) == ncol(cbio.mat)


# Keep in mind that all discovery and validation samples are apparently primary tumor data (looking at overlap with CBioPortal confirmed this)
table(clin.mat$CANCER_TYPE,exclude = NULL)
# Also includes PAM50 classification of tumor subtype status
table(clin.mat$Pam50_SUBTYPE,exclude = NULL)
# The same subtype classification, but using the "THREE Gene" classifider (ER/HER2/Proliferation)
table(clin.mat$THREEGENE,exclude = NULL)

# Match order of samples in count matrix and phenotype data
cbio.mat <- cbio.mat[,sort(colnames(cbio.mat))]
clin.mat <- clin.mat[order(clin.mat$PATIENT_ID),]

rownames(clin.mat) <- clin.mat$PATIENT_ID
# Check if they align
all.equal(colnames(cbio.mat),rownames(clin.mat))

# Assemble ESet object now that everything is in place
ES <- ExpressionSet(assayData = cbio.mat,
                    annotation = "Expression log intensity levels (Illumina Human v3 microarray). Targeted sequencing of 1980 primary breast cancer samples. (Pereira et al. Nat Commun 2015)",
                    phenoData = AnnotatedDataFrame(clin.mat),
                    featureData = AnnotatedDataFrame(data.frame(cbio.d[,c(1,2)],row.names=rownames(cbio.mat)))
)

ES

saveRDS(ES,"/restricted/projectnb/montilab-p/CBMrepositoryData/CURTIS/ESets/metabrick_GE_ESet.rds")

usethis::use_data(metabric_data, overwrite = TRUE)
