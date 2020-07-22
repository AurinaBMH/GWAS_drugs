function data = VisualizeWeightingVectors(whichDiseases, whatMeasurement, similarityType, whatProperty)

if nargin < 3
    similarityType = 'MAGMAdefault'; 
end
if nargin < 4
    whatProperty = 'P'; 
end
numTop = 50; % Look at the top X genes in particular
orderByWhat = 'max';
cMapMax = 0.03;

%-------------------------------------------------------------------------------
numDiseases = length(whichDiseases);
geneNames = cell(numDiseases,1);
geneWeightsNorm = cell(numDiseases,1);
for i = 1:numDiseases
    whatDisease = whichDiseases{i};
    [geneNames{i},geneWeightsNorm{i}] = GiveMeNormalizedScoreVector(whatDisease,...
                                    whatMeasurement,similarityType,whatProperty);
end

% Combine multiple diseases:
% [~,ia,ib] = intersect(whichDiseases,indicatorTable.Properties.VariableNames,'stable');
% indicatorTable = indicatorTable(:,ib);
data = horzcat(geneWeightsNorm{:})';

% Order by maximum weight:
switch orderByWhat
case 'max'
    maxGene = max(data,[],1);
case 'mean'
    maxGene = mean(data,1);
end
[maxGeneSort,ix] = sort(maxGene,'descend');
ix(isnan(maxGeneSort)) = []; % remove top genes that are actually NaNs
data = data(:,ix(1:numTop));
geneNames = geneNames{1}(ix(1:numTop));

%-------------------------------------------------------------------------------
f = figure('color','w');
imagesc(data)
colormap([1,1,1;BF_getcmap('blues',9,false)])
% colormap([1,1,1;jet(12)])
ax = gca;
ax.YTick = 1:numDiseases;
ax.YTickLabel = whichDiseases;
ax.XTick = 1:numTop;
ax.XTickLabel = geneNames;
ax.XTickLabelRotation = 90;
cB = colorbar;
cB.Label.String = 'Normalized treatment weight';
%caxis([0,cMapMax])
% Title:
extraText = '';
if strcmp(whatMeasurement,'GWAS')
    extraText = sprintf('%s-%s',similarityType,whatProperty);
end
title(sprintf('%s weightings over genes (%s)',whatMeasurement,extraText),...
                    'interpreter','none')
f.Position = [855   870   569   234];
end
