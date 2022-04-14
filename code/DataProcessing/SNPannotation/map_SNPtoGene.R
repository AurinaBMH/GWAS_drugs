library("biomaRt")
library("dplyr")


setwd("~/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists/GWAS")
# load data that will be used for mapping: SNP-level and gene level
grch38.snp = useMart(biomart="ENSEMBL_MART_SNP", host="ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")
grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# this mapping requires rs IDs for each SNP, use the original GWAS summary statistics from PGC
# read GWAS
# For all disorders do clumping using 1000GPhase 3 dataset to keep independent SNPs;
# the list for ADHD sataset contains SNPs in LD with each other, need to clump and keep only independent ones; 

DISORDERS=c('ADHD', 'AUT', 'BIP2', 'IQSavage', 'MDD2', 'SCZ', 'DIABETES')
pThr = 10^-5
for (disorder in DISORDERS){
  
# load the dataset 'AUT', have independent SNPs in the original file
if (disorder=='AUT'){
  fileName = sprintf('pgc%s.txt', disorder)
# for others load clumped data
} else {
  fileName = sprintf('pgc%s.clumped', disorder)}

  
GWASlist <- read.delim(fileName, sep="")

# select SNPs at a selected p threshold, Janette jused 10^-5
GWASlist = GWASlist[GWASlist$P<pThr,]
SNPlist = GWASlist$SNP
# find only SNPs that have rs in their name
keepSNP = grepl("rs", GWASlist$SNP, fixed = TRUE)
SNPlistKeep = SNPlist[keepSNP]

sprintf("%s: unique SNPs p<%e identified n=%d", fileName, pThr, length(SNPlistKeep))

# extract information relevant for SNPs for a selected list
SNPtable <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"), 
                filters = "snp_filter",
                values = SNPlistKeep, 
                mart = grch38.snp,   
                useCache = FALSE, 
                uniqueRows = TRUE)

# extract information relevant for genes for a selected list
GENEtable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "ensembl_gene_id", 
                values =  SNPtable$ensembl_gene_stable_id, 
                mart = grch38, 
                useCache = FALSE, 
                uniqueRows = TRUE)

results <- merge(SNPtable,GENEtable, by.x = "ensembl_gene_stable_id", by.y="ensembl_gene_id", all.x=T)
# keep only unique rows of the data wirh SNPs with all data in each row
FINALresults = unique(results[complete.cases(results), 2:3])


write.table(FINALresults, file = sprintf("pgc%s_genes.txt",disorder), 
            sep = "\t", row.names = FALSE, quote=FALSE)

sprintf("%s: unique SNPs identified n=%d", fileName, nrow(FINALresults))
}
