% GenerateResultsTables
% Idea is to run through pipeline for each disease, and save results tables
% for each:
function GenerateResultsTables()

% Diseases with SNP info, etc.:
params = SetDefaultParams();
whatDiseases = params.whatGWAS_2022;
whatThreshold='BF';

% load information on genes mapped using different methods
load('DataOutput_2022/GWAS_disordersMAGMA.mat')

%-------------------------------------------------------------------------------
numDiseases = length(whatDiseases);
for k = 1:numDiseases
    %for k = 1:numDiseases
    whatDisease = whatDiseases{k};
    geneScores = pipeline(DISORDERlist, whatDisease, whatThreshold);
    % Save:
    fileName = sprintf('resultsTable_%s_%s_%s_%s_drugbank.mat',whatDisease, whatThreshold, params.whatDrugTargets, params.whatTargets);
    fileName = fullfile('DataOutput_2022',fileName);
    save(fileName,'geneScores');
    fprintf(1,'Saved results to %s!!!\n\n\n',fileName);
end
end
