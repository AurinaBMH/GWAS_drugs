library("biomaRt")
library("dplyr")

setwd("/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/enrichment/")
# load data that will be used for mapping: SNP-level and gene level
grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

geneNames <- read.delim('/Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/enrichment/drugbank_geneNames.txt', sep="")

# extract information relevant for genes for a selected list
GENEtable <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                filters = "external_gene_name", 
                values =  geneNames, 
                mart = grch38, 
                useCache = FALSE, 
                uniqueRows = TRUE)

results <- merge(SNPtable,GENEtable, by.x = "ensembl_gene_stable_id", by.y="ensembl_gene_id", all.x=T)
# keep only unique rows of the data wirh SNPs with all data in each row
FINALresults = unique(results[complete.cases(results), 2:3])


write.table(FINALresults, file = sprintf("pgc%s_genes.txt",disorder), 
            sep = "\t", row.names = FALSE, quote=FALSE)

