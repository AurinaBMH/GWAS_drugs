% Loop over different versions of the GWAS-gene mapping and select
% combination that gives the strongest match between GWAS and drugs.
clear all; close all; 

load('resultsTable_ADHD_FDR.mat', 'geneScores')

% select all available similarity types
similarityTypes = setdiff(fieldnames(geneScores), {'gene', 'params'});

% for PPI select 1 neighbor and mean and median measures;
whatDiseases_GWAS = {'BIP2'};
whatDiseases_Drug = 'BIP';
whatNull = 'randomDisease';
PPImeasures_names = {'numPPIneighbors1','percPPIneighbors1','weiPPIneighbors1','expWeiPPIneighbors1', 'numPPIneighbors2','percPPIneighbors2'}; 

compareGWASvsDRUGmatches(whatDiseases_GWAS, whatNull, similarityTypes, PPImeasures_names, whatDiseases_Drug); 

