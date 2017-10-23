% Computes results tables across a range of diseases (iteratively):
%-------------------------------------------------------------------------------
diseaseList = {'all','SZP','ASD','ADHD','BIP','MDD'};
PPINevidenceThreshold = 0;
numDiseases = length(diseaseList);

resultsTable = struct();
for i = 1:numDiseases
    resultsTable.(diseaseList{i}) = pipeline(diseaseList{i},PPINevidenceThreshold);
end
