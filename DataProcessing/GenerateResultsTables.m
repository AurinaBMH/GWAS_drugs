% Idea is to run through pipeline for each disease, and save results tables
% for each:

% Diseases with SNP info, etc.:
whatDiseases = {'ADHD','BIP','SZP','MDD','diabetes'};
PPINevidenceThreshold = 0.4;

%-------------------------------------------------------------------------------
numDiseases = length(whatDiseases);
for k = 1:numDiseases
    whatDisease = whatDiseases{k};
    resultsTable = pipeline(whatDisease,PPINevidenceThreshold);
    % Save:
    fileName = sprintf('resultsTable_%s-%.1f.mat',whatDisease,PPINevidenceThreshold);
    fileName = fullfile('DataOutput',fileName);
    save(fileName,'resultsTable');
    fprintf(1,'Saved results to %s!!!\n\n\n',fileName);
end
