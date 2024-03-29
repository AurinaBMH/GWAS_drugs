# get BIOMART file using querie
# first, will need to install the corresponding packages (do it only once), uncomment if running on a different computer

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

#BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
#BiocManager::install(c("biomaRt"))

# load libraries
library("biomaRt")
library("dplyr")

setwd("~/Google_drive/PostDoc/projects/GWASdrugs")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

DATAtable = getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'external_gene_name', 'hgnc_symbol', 'hgnc_id', 'entrezgene_id'), 
      mart = ensembl)

write.table(DATAtable, file = sprintf("data/GWASlists/BIOMART_geneIDs.txt"), 
            sep = "\t", row.names = FALSE, quote=FALSE)
