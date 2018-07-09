function GiveMeComparisonScores(whatDiseaseGWAS,whatProperty,whatDiseaseDrug)
%-------------------------------------------------------------------------------
% Returns similarity scores from matching
%-------------------------------------------------------------------------------

if nargin < 1
    whatDiseaseGWAS = 'SZP';
end
if nargin < 2
    whatProperty = 'percPPIneigh1LD';
end
if nargin < 3
    whatDiseaseDrug = 'SZP';
end

%-------------------------------------------------------------------------------
% Load data:
%-------------------------------------------------------------------------------
% (1) Similarity of each gene to GWAS hits for this disorder:
fileName = sprintf('resultsTable_%s-0.0.mat',whatDiseaseGWAS);
load(fileName,'resultsTable');

% (2) Treatment weights of each gene implicated in the disorder:
normalizeWithinDrugs = true; % weight genes lower if they occur in drugs with large numbers of gene targets
[indicatorTable,percIndicatorTable] = ImportTreatmentLists(normalizeWithinDrugs);

% Combine two datasets on overlap:
[~,ia,ib] = intersect(resultsTable.gene,indicatorTable.Properties.RowNames,'stable');
indicatorTable = indicatorTable(ib,:);
resultsTable = resultsTable(ia,:);

%-------------------------------------------------------------------------------
% Get the property of interest from the GWAS-hit-similarity table:
geneWeights_GWAS = resultsTable.(whatProperty);
% Get weights of each gene onto treatment targets:
geneWeights_treatment = indicatorTable.(whatDiseaseDrug);

% Keep only non-NaN values and rescale both to unit interval:
bothGood = (~isnan(geneWeights_GWAS) & ~isnan(geneWeights_treatment));
f_unit = @(x)BF_NormalizeMatrix(x,'mixedSigmoid'); %/sum(x);
geneWeights_GWAS = f_unit(geneWeights_GWAS(bothGood));
geneWeights_treatment = f_unit(geneWeights_treatment(bothGood));
geneNames = resultsTable.gene(bothGood);

[rho,p] = corr(geneWeights_GWAS,geneWeights_treatment,'type','Spearman');
title(sprintf('rho = %.2f (p = %.3g)\n',rho,p));

% Compute similarity score for both (contribution to correlation):
simScore = (geneWeights_GWAS + geneWeights_treatment);
[~,ix] = sort(simScore,'descend');

%-------------------------------------------------------------------------------
% Scatter plot:
f = figure('color','w'); hold on
plot(geneWeights_GWAS,geneWeights_treatment,'ok')
xlabel(sprintf('gene score—-%sGWAS-%s',whatDiseaseGWAS,whatProperty))
ylabel(sprintf('gene score—-%sdrugs',whatDiseaseDrug))
[rho,p] = corr(geneWeights_GWAS,geneWeights_treatment,'type','Spearman');
title(sprintf('rho = %.2f (p = %.3g)\n',rho,p));

for j = 1:20
    ep = 0.01;
    text(ep+geneWeights_GWAS(ix(j)),ep+geneWeights_treatment(ix(j)),geneNames{ix(j)});
end

end
