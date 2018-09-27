% Parameters:
whichDiseases = {'ADHD','MDD','BIP','SZP'}; % 'diabetes','pulmonary','cardiology','gastro'};
% whatMeasurement = 'Drug';
whatMeasurement = 'GWAS';
% similarityType = 'DNA';
% whatProperty = 'percGWAS';
similarityType = 'PPI_th0';
whatProperty = 'percPPIneighbors1';
numTop = 30; % Look at the top X genes in particular

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
% maxGene = max(data,[],1);
maxGene = mean(data,1);
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
caxis([0,0.15])
% Title:
extraText = '';
if strcmp(whatMeasurement,'GWAS')
    extraText = sprintf('%s-%s',similarityType,whatProperty);
end
title(sprintf('%s weightings over genes (%s)',whatMeasurement,extraText),...
                    'interpreter','none')
f.Position = [440   427   887   371];
