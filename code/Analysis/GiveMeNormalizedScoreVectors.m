function [geneNames,data] = GiveMeNormalizedScoreVectors(whichDiseases,whatMeasurement,similarityType,whatProperty,orderByWhat)
% Can provide multiple diseases, will loop over them:
%-------------------------------------------------------------------------------
if nargin < 3
    similarityType = '';
end
if nargin < 4
    whatProperty = '';
end
if nargin < 5
    orderByWhat = 'none';
end
%-------------------------------------------------------------------------------

% Retrieve:
numDiseases = length(whichDiseases);
geneNames = cell(numDiseases,1);
geneWeightsNorm = cell(numDiseases,1);
for i = 1:numDiseases
    whatDisease = whichDiseases{i};
    [geneNames{i},geneWeightsNorm{i}] = GiveMeNormalizedScoreVector(whatDisease,whatMeasurement,similarityType,whatProperty);
end

% Agglomerate:
data = horzcat(geneWeightsNorm{:});
geneNames = geneNames{1};

% Sort:
% switch orderByWhat
% case 'none'
%     return
% case 'max'
%     maxGene = max(data,[],1);
% case 'mean'
%     maxGene = mean(data,1);
% end
% [maxGeneSort,ix] = sort(maxGene,'descend');
% ix(isnan(maxGeneSort)) = []; % remove top genes that are actually NaNs
% data = data(:,ix);
% geneNames = geneNames(ix);

end
