
## Workflow
### Data processing

#### Set default parameters
All default parameters to be used repeatedly are to be set in `SetDefaultParams`.

#### PPI network data
Process PPI network data, save to processed sparse matrix:
```
    % Generate a weighted PPI network with evidence thresholds as weights
    PPINImport(true);
    % Generate a binary PPI network thresholded at evidence threshold of 0.4
    PPINImport(false,400);
```

### Analysis

#### Compute SNP annotation details

Reads in information about SNPs from GWAS hits, computes LD neighbor information
from a mySQL database, saves output to `SNPAnnotationTable_X` for a given LD threshold X:
```
    SNPAnnotationProcessing;
```

#### Characterize drugs used for each disorder with respect to GWAS hits
Run:
```
    GenerateResultsTables
```
Saves `resultsTable_X.mat` for each disorder, X.
These `.mat` files contain a structure with diverse information about how each
drug-target gene relates to the set of GWAS hits for a given disorder X.

### Visualization
#### `DistinguishingCharBar`
Allows visualization of how a set of drug-target genes relates to GWAS hits for a given disorder.


## Data Information
`eQTL_protein_names.csv` contains protein nodes for the entire PPI network.

`eQTL_edges.csv` contains all interactions formed by the nodes. For a given interaction, i code whether the nodes are from first degree partner of disease genes (binary with 1 denoting first degree partner, 0 for disease gene), and which disease is the disease gene from.

`eQTL_identifier.csv` contains which disease is the protein from, the association type (LD or GWAS relationship) with the disease, whether it is a first degree interaction partner or not.
