% Idea here is to be able to be able to zoom in on a distribution and look at the
% contributing genes
%-------------------------------------------------------------------------------

whatDiseaseGWAS = 'SZP';
whatProperty = 'percPPIneighbors1DiseaseLD';
whatDiseaseDrug = 'SZP';

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
f = figure('color','w');
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

% Compute dot-product similarity:
simScore = (geneWeights_GWAS+geneWeights_treatment);
% simScoreNonZero = simScore(simScore > 0);
% geneNamesNonZero = resultsTable.gene(simScore > 0);
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

%-------------------------------------------------------------------------------
% Plot distribution, labeled by genes:
f = figure('color','w');
[f,x] = ksdensity(simScore,'npoints',500); %,'support',[0,max(simScore)+0.1]);
plot(f,x,'k')
hold on
minAnnotate = 2;
for j = 1:length(simScore)
    % r = rand(1);
    r = rand(1)*f(find(x >= simScore(j),1));
    plot(r,simScore(j),'.k');
    if abs(simScore(j)) > minAnnotate
        text(r,simScore(j),geneNames{j})
    end
end
ylabel(sprintf('%s-GWAS (%s), %s-drug similarity score',whatDiseaseGWAS,whatProperty,whatDiseaseDrug))
title(sprintf('%s-%s (%s drugs)',whatDiseaseGWAS,whatProperty,whatDiseaseDrug))
ax = gca();
% ax.YLim = [0,ax.YLim(2)];

%-------------------------------------------------------------------------------
% What about looking at how the barcodes match up:
numTop = 500;
barcodes = [geneWeights_GWAS,geneWeights_treatment]';
[~,ix] = sort(simScore,'descend');
f = figure('color','w');
imagesc(BF_NormalizeMatrix(barcodes(:,ix(1:numTop))','mixedSigmoid')')
ax = gca;
ax.YTick = 1:2;
ax.YTickLabel = {'GWAS-score','Treatment-score'};
ax.XTick = 1:numTop;
ax.XTickLabel = geneNames(ix(1:numTop));
