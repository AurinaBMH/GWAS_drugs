# get BIOMART file using querie
setwd("~/Google_drive/PostDoc/projects/GWASdrugs/data/GWASlists")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

DATAtable = getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id', 'external_gene_name', 'hgnc_symbol', 'hgnc_id', 'entrezgene_id'), 
      mart = ensembl)

write.table(DATAtable, file = sprintf("BIOMART_geneIDs.txt"), 
            sep = "\t", row.names = FALSE, quote=FALSE)
