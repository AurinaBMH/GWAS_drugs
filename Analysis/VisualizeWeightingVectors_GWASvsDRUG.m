% plot drug and gene scores for each disorder
function dataALL = VisualizeWeightingVectors_GWASvsDRUG(whatDiseaseGWAS, whatDiseaseDRUG, similarityType, whatProperty, plotSeparate)

if nargin <3
    similarityType = 'MAGMAdefault';
    fprintf('Ploting MAPPING %s by DEFAULT', similarityType)
end
if nargin <4
    whatProperty = 'P';
    fprintf('Ploting MEASURE %s by DEFAULT', whatProperty)
end
if nargin <5
    plotSeparate = 1;
end
numTop = 100;

[geneNamesDrug,geneWeightsNormDrug] = GiveMeNormalizedScoreVector(whatDiseaseDRUG,'Drug');
                                
[geneNamesGWAS,geneWeightsNormGWAS] = GiveMeNormalizedScoreVector(whatDiseaseGWAS,...
                                    'GWAS',similarityType,whatProperty);
                                
[geneNames,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug,'stable');
geneWeightsNormGWAS = geneWeightsNormGWAS(ia);
geneWeightsNormDrug = geneWeightsNormDrug(ib);

data = horzcat(geneWeightsNormGWAS, geneWeightsNormDrug)';
maxGene = max(data,[],1);

[maxGeneSort,ix] = sort(maxGene,'descend');
ix(isnan(maxGeneSort)) = []; % remove top genes that are actually NaNs
dataALL = data(:,ix);
data = data(:,ix(1:numTop));

geneNames = geneNames(ix(1:numTop));
if plotSeparate
f = figure('color','w');
end
imagesc(data)
colormap([1,1,1;BF_getcmap('blues',9,false)])
ax = gca;
ax.YTick = 1:2;
ax.YTickLabel = {'GWAS'; 'Drug'};
ax.XTick = 1:numTop;
ax.XTickLabel = geneNames;
ax.XTickLabelRotation = 90;
cB = colorbar;
cB.Label.String = 'Normalized treatment weight';
caxis([0,max(data(1,:))])
% Title:
extraText = sprintf('%s-%s',similarityType,whatProperty);
title(sprintf('%s weightings over genes (%s)',whatDiseaseGWAS,extraText),...
                    'interpreter','none')
f.Position = [855   870   569   234];
end