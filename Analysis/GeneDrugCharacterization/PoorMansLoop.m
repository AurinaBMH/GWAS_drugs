% Computes results tables across a range of diseases (iteratively):
%-------------------------------------------------------------------------------
diseaseList = {'all','SZP','ASD','ADHD','BIP','MDD'};
PPINevidenceThreshold = 0;
numDiseases = length(diseaseList);

resultsTable = struct();
for i = 1:numDiseases
    resultsTable.(diseaseList{i}) = pipeline(diseaseList{i},PPINevidenceThreshold);
end

%-------------------------------------------------------------------------------
% Display some results:
%-------------------------------------------------------------------------------
whatDisease = 'SZP';
% Display just with custom columns
customColumns = {'gene','numGWASMapped','numLDSNPs','numPPIneighbors1DiseaseMapped',... %,'meanPPIDistance'
                'numPPIneighbors1DiseaseLD','percPPIneighbors1DiseaseMapped',...
                'percPPIneighbors1DiseaseLD','matchingDrugsString'};
display(resultsTable.(whatDisease)(1:60,ismember(resultsTable.(whatDisease).Properties.VariableNames,customColumns)));
