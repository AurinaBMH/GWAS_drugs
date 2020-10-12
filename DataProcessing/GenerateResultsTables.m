% GenerateResultsTables
% Idea is to run through pipeline for each disease, and save results tables
% for each:
clear all; close all; 
% Diseases with SNP info, etc.:
params = SetDefaultParams();
whatDiseases = params.whatGWAS; 

whatThreshold='BF'; 
% load information on genes mapped using different methods
load('GWAS_disordersMAGMA.mat')
params = SetDefaultParams();

%-------------------------------------------------------------------------------
numDiseases = length(whatDiseases);
for k = 1:numDiseases
    whatDisease = whatDiseases{k};
    geneScores = pipeline(DISORDERlist, whatDisease, whatThreshold);
    % Save:
    fileName = sprintf('resultsTable_%s_%s_%s_drugbank.mat',whatDisease, whatThreshold, params.whatDrugTargets);
    fileName = fullfile('DataOutput',fileName);
    save(fileName,'geneScores');
    fprintf(1,'Saved results to %s!!!\n\n\n',fileName);
end
