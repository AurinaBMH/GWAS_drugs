% Idea is to get disease-unique genes after comparing interactions
% across all five disorders
%-------------------------------------------------------------------------------

diseases = {'ADHD','ASD','BIP','MDD','SZP'};

%-------------------------------------------------------------------------------
% Load data:
load(fullfile('Data','processedData.mat'),'Adj','eQTLproteinnames');
eQTLidentifier = importIdentifier();


%-------------------------------------------------------------------------------
numDiseases = length(diseases);
diseaseMarker = false(length(eQTLproteinnames),numDiseases);

for i = 1:numDiseases
    diseaseMarker(:,i) = ismember(eQTLproteinnames,...
                unique(eQTLidentifier.Name(eQTLidentifier.(diseases{i}) & ~eQTLidentifier.Partners)));
end

% Filter by genes that mark for disorders (i.e., no partners):
isAny = any(diseaseMarker,2);
diseaseMarker_filt = diseaseMarker(isAny,:);
eQTLproteinnames_filt = eQTLproteinnames(isAny);
Adj_filter = Adj(isAny,isAny);

%-------------------------------------------------------------------------------
% Get unique genes for each disorder:
for i = 1:numDiseases
    isD = diseaseMarker_filt(:,i);
    notOthers = ~any(diseaseMarker_filt(:,setxor(1:numDiseases,i)),2);
    isD_specific = isD & notOthers;
    uniqueD = countUnique(isD_specific,notOthers,Adj_filter);
    fprintf(1,'%u/%u %s - []\n',sum(uniqueD),sum(isD),diseases{i});
end
