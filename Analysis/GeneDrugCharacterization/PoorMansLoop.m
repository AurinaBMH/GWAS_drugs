% Computes results tables across a range of diseases (iteratively):
%-------------------------------------------------------------------------------
diseaseList = {'all','SZP','ASD','ADHD','BIP','MDD'};
PPINevidenceThreshold = 0.4;
numDiseases = length(diseaseList);

resultsTable = struct();
for i = 1:numDiseases
    resultsTable.(diseaseList{i}) = pipeline(diseaseList{i},PPINevidenceThreshold);
end

% Save results:
fileNameSave = sprintf('AllResultsTables_th_0%s.mat',PPINevidenceThreshold);
fileNameSave = fullfile('Data',fileNameSave);
save(fileNameSave,'resultsTable')
fprintf(1,'Saved all results tables to %s\n',fileNameSave);

%-------------------------------------------------------------------------------
% Display some results:
%-------------------------------------------------------------------------------
whatDisease = 'SZP';
% Display just with custom columns
customColumns = {'gene','numGWASMapped','numLDSNPs',... %,'meanPPIDistance'
                'numPPIneighbors1DiseaseLD','percPPIneighbors1DiseaseLD','matchingDrugsString'};
display(resultsTable.(whatDisease)(1:60,ismember(resultsTable.(whatDisease).Properties.VariableNames,customColumns)));
