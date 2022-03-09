
This repository provides Matlab, R and MAGMA code for reproducing results presented in the manuscript entitled: 

- Arnatkeviciute et al. (2022) [:green_book: 'Linking GWAS to pharmacological treatments for psychiatric disorders'](DOI).

The code was written using MATLAB_R2020b. 

Contact Aurina Arnatkeviciute by [email](mailto:aurina.arnatkeviciute@monash.edu).


## Workflow
### Data processing
First, add all sub-folders to the path using startup() function from the root directory. 

#### :label: Aggregate PPI-based information
Information on the PPI data (file: `9606.protein.links.v11.0.txt`) can be found in `rawData/README_PPI.txt`
1. Replace protein IDs with gene names using `make_PPI_linkfile()`; 
2. Generate a binary PPI network thresholdeds at different evidence thresholds: 0,400,600,900: 
```matlab
% BINARY networks:
PPIthrs = [0,400,600,900];
for t=1:length(PPIthrs)   
    % save PPI Adj and distance matrix
    [AdjPPI,geneNames] = PPINImport(false, PPIthrs(t), 'HGNCmatch');
    distMatrix = ComputePPIDist(PPIthrs(t), false);  
end
% WEIGHTED network:
[AdjPPI,geneNames] = PPINImport(true);
```
These commands will save `PPI_HGNC_Adj_th0.mat/PPI_HGNC_Dist_th0.mat/PPI_HGNC_geneLabels_th0.mat` files.   

#### :label: Aggregate GWAS-based information

1. Map genes based on GWAS summary statistics for each disorder using `HMAGMA_code_2022.sh`. 
First, modify paths in lines 1-5 of `code/DataProcessing/HMAGMA/HMAGMA_code_2022.sh` to indicate the location of code, .annot files and reference genome. 

2. Gene names are in the `ENSG` format. Get gene name to entrezID mapping using `code/DataProcessing/HMAGMA/get_BIOMARTdata.R`.
The output is saved to `BIOMART_geneIDs.txt`; 

3. Update gene IDs for MAGMA outputs and collate all results into a single .mat file: 
```matlab
save_MAGMAHresults()
```

4. Create GWAS-based gene scores and save them for each disorder: 
```matlab
GenerateResultsTables()
```
This will create `geneScores` structure for each disorder (takes several hours to run). 

#### :label: Aggregate drug target information

1. Combine drug target information from `.txt` files into matlab format
```matlab
dataTable = give_drugTargets('all', 'drugbank'); 
```
This will save `drugTargets_2020_all_drugbank.mat` file; 

2. Generate 5000 drug-based null vectors for each disorder. 
For each disorder a corresponding number of random drugs is selected and treatment-based scores are calculated across all 2155 genes;   
For example, there are 14 drugs for ADHD, 22 for bopolar disorder and 45 for diabetes, so for each disorder that number of random treatments is selected. 
```matlab
generate_randomDrug_nulls('drugbank')
```


### Analysis

#### :scroll: Reproduce results presented in the manuscript: 

For psychiatric disorders and diabetes (`Figure 2`, `Figure 3`, `Figure S1`, `Figure S2`, `Figure S4`) generate figures using: 
```matlab
plot_Psych_figures()
```

For non-psychiatric disorders (`Figure S3`): 
```matlab
plot_Body_figures()
```

Save gene scores for enrichment analysis: 
```matlab
save_enrichment_scores()
```

Run the enrichment analysis using ermineJ software and aggregate results using: 
```matlab
save_enrichment_results(); 
```

## Data Information

#### :dna: GWAS summary statistics: 
1. ADHD GWAS summary statistics based on [:green_book: 'Demontis et al (2019)'](https://doi.org/10.1038/s41588-018-0269-7)
2. Bipolar disorder GWAS summary statistics based on [:green_book: 'Mullins et al (2021)'](https://doi.org/10.1038/s41588-021-00857-4)
3. Major depression GWAS summary statistics based on [:green_book: 'Howard et al (2019)'](https://doi.org/10.1038/s41593-018-0326-7)
4. Schizophrenia GWAS summary statistics based on [:green_book: 'The Schizophrenia Working Group of the Psychiatric Genomics Consortium et al (2020)'](https://doi.org/10.1101/2020.09.12.20192922)
5. Diabetes GWAS summary statistics based on [:green_book: 'Xue et al (2018)'](https://doi.org/10.1038/s41467-018-04951-w)
6. Heart failure GWAS summary statistics based on [:green_book: 'Shah et al (2021)'](https://doi.org/10.1038/s41467-019-13690-5)
7. Inflammatory bowel disease GWAS summary statistics based on [:green_book: 'Lange et al (2017)'](https://doi.org/10.1038/ng.3760)
8. Rheumatoid arthritis GWAS summary statistics based on [:green_book: 'Okada et al (2013)'](https://doi.org/10.1038/nature12873)


#### :books: PPI network data

`9606.protein.links.v11.0.txt.gz (71.2 Mb)` and `9606.protein.info.v11.0.txt.gz (1.9 Mb)` - downloaded from [:books: 'STRING database (version 11.0)'](https://string-db.org/cgi/download.pl?sessionId=a1fHJhN5R9Md&species_text=Homo+sapiens) on the 24th of June 2020;

#### :pill: Treatments lists
Treatments for different conditions of interest were selected by searching the [:medical_symbol: 'DrugBank database'](www.drugbank.ca), accessed on September 3, 2020. 
Specifically, drugs for each indication were searched in the DrugBank database using the following search terms: 
1. "attention deficit" (for ADHD);
2. "bipolar" (for bipolar disorder) excluding "bipolar depression";
3. "schizophrenia" (for schizophrenia); 
4. "major depression" (for major depression); - confirm this term!!! 
5. "diabetes" (for type 2 diabetes) excluding "type I diabetes" and "diabetes insipidus"; 
6. "heart failure" (for heart failure); 
7. "Crohn's" and "ulcerative colitis" (for inflammatory bowel disease"); and 
8. "rheumatoid arthritis" (for rheumatoid arthritis).
