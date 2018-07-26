% GenerateResultsTables
% Idea is to run through pipeline for each disease, and save results tables
% for each:

% Diseases with SNP info, etc.:
whatDiseases = {'ADHD','BIP','SZP','MDD','diabetes'};
params = SetDefaultParams();

%-------------------------------------------------------------------------------
numDiseases = length(whatDiseases);
for k = 1:numDiseases
    whatDisease = whatDiseases{k};
    geneScores = pipeline(whatDisease);
    % Save:
    fileName = sprintf('resultsTable_%s.mat',whatDisease);
    fileName = fullfile('DataOutput',fileName);
    save(fileName,'geneScores');
    fprintf(1,'Saved results to %s!!!\n\n\n',fileName);
end
