
## Workflow
### Data processing
First, add all sub-folders to the path using startup() function from the root directory. 

#### Aggregate PPI-based information
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

#### Aggregate GWAS-based information

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
This will create `geneScores` structure for each disorder. This step takes several hours to run. 

#### Aggregate drug target information

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

#### Reproduce results presented in the manuscript: 

1. For psychiatric disorders and diabetes (Figure 2, Figure 3, Figure S1, Figure S2, Figure S4) generate figures using: 
```matlab
plot_Psych_figures
```

2. For non-psychiatric disorders (Figure S3): 
```matlab
plot_Body_figures
```

3. Save gene scores for enrichment analysis: 
```matlab
save_enrichment_scores()
```

4. Run the enrichment analysis using ermineJ software and aggregate results using: 
```matlab
GOtable = save_enrichment_results(); 
```

## Data Information
