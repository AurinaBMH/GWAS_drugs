% --------------------------------------------
% Aggregating GWAS-based information
% --------------------------------------------
% updating data using new GWAS for MDD, BIP and SCZ

% 1. run /Users/aurinaa/Documents/PostDoc/projects/GWASdrugs/code/DataProcessing/HMAGMA/HMAGMA_code_2022.sh 
% to create new MAGMA outputs

% 2.0 update default disorder abbreviations and paths for MAGMA outputs in
% SetDefaultParams(); 

% 2.1 update gene IDs for MAGMA outputs and collate all results into a
% single .mat file using save_MAGMAHresults() function; 

% 3.1 creating results tables for each disorder using
% GenerateResultsTables() function; This will create geneScores structure for each disorder. 
% This step takes the longest to run; 

% --------------------------------------------
% Aggregating drug target information
% --------------------------------------------
% 1. Combine drug target information from .txt files into matlab format
% dataTable = give_drugTargets('all', 'drugbank')

% 2. Generate nulls for each disorder - 5000 drug-based vectors selecting X number of random drugs: 
% generate_randomDrug_nulls('drugbank')

% --------------------------------------------
% Reproducing results
% --------------------------------------------
% Generate figures using: plot_Psych_figures.m
% For non-psychiatric disorders: plot_Body_figures.m
