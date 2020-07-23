% plot drug and gene scores for each disorder
function VisualizeWeightingVectors_GWASvsDRUG(whatDiseaseGWAS, whatDiseaseDRUG, similarityType, whatProperty)

if nargin <3
    similarityType = 'MAGMAdefault';
    fprintf('Ploting MAPPING %s by DEFAULT', similarityType)
end
if nargin <4
    whatProperty = 'P';
    fprintf('Ploting MEASURE %s by DEFAULT', whatProperty)
end
numTop = 100;

[geneNamesDrug,geneWeightsNormDrug] = GiveMeNormalizedScoreVector(whatDiseaseGWAS,'Drug');
                                
[geneNamesGWAS,geneWeightsNormGWAS] = GiveMeNormalizedScoreVector(whatDiseaseDRUG,...
                                    'GWAS',similarityType,whatProperty);
                                
[geneNames,ia,ib] = intersect(geneNamesGWAS,geneNamesDrug,'stable');
geneWeightsNormGWAS = geneWeightsNormGWAS(ia);
geneWeightsNormDrug = geneWeightsNormDrug(ib);

data = horzcat(geneWeightsNormGWAS, geneWeightsNormDrug)';
maxGene = max(data,[],1);

[maxGeneSort,ix] = sort(maxGene,'descend');
ix(isnan(maxGeneSort)) = []; % remove top genes that are actually NaNs
data = data(:,ix(1:numTop));
geneNames = geneNames(ix(1:numTop));

f = figure('color','w');
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
%caxis([0,cMapMax])
% Title:
extraText = sprintf('%s-%s',similarityType,whatProperty);
title(sprintf('%s weightings over genes (%s)',whatDisease,extraText),...
                    'interpreter','none')
f.Position = [855   870   569   234];
end