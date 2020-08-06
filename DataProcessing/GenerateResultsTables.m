% GenerateResultsTables
% Idea is to run through pipeline for each disease, and save results tables
% for each:

% Diseases with SNP info, etc.:
whatDiseases = {'ADHD', 'MDD2', 'SCZ', 'BIP2', 'DIABETES', 'HF', 'AD'};
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
    fileName = sprintf('resultsTable_%s_%s.mat',whatDisease, whatThreshold);
    fileName = fullfile('DataOutput',fileName);
    save(fileName,'geneScores');
    fprintf(1,'Saved results to %s!!!\n\n\n',fileName);
end
