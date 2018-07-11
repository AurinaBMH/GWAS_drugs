% Idea here is to be able to be able to zoom in on a distribution and look at the
% contributing genes
%-------------------------------------------------------------------------------

whatDiseaseGWAS = 'SZP';
whatProperty = 'percPPIneigh1Mapped';
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
% Get the property of interest from the GWAS-hit-similarity table:
geneWeights_GWAS = resultsTable.(whatProperty);
% Get weights of each gene onto treatment targets:
geneWeights_treatment = indicatorTable.(whatDiseaseDrug);

% Keep only non-NaN values and rescale both to unit interval:
bothGood = (~isnan(geneWeights_GWAS) & ~isnan(geneWeights_treatment));
geneWeights_GWAS = geneWeights_GWAS(bothGood);
geneWeights_treatment = geneWeights_treatment(bothGood);
geneNames = resultsTable.gene(bothGood);

%-------------------------------------------------------------------------------
% Compute combined similarity score for each gene:
f_unit = @(x)BF_NormalizeMatrix(x,'mixedSigmoid'); %/sum(x);
simScore = (f_unit(geneWeights_GWAS)+f_unit(geneWeights_treatment));
% simScoreNonZero = simScore(simScore > 0);
% geneNamesNonZero = resultsTable.gene(simScore > 0);
% [~,ix] = sort(simScore,'descend');

%-------------------------------------------------------------------------------
% Scatter plot:
f = figure('color','w'); hold on
plot(geneWeights_GWAS,geneWeights_treatment,'ok')

% Annotations:
usedInTreatment_ind = find(geneWeights_treatment > 0);
for k = 1:length(usedInTreatment_ind)
    ind_k = usedInTreatment_ind(k);
    % Annotate lines:
    plot(ones(2,1)*geneWeights_GWAS(ind_k),[0,geneWeights_treatment(ind_k)],'-k')
    % Annotate text labels:
    ep = 0.01*max(geneWeights_treatment);
    text(ep+geneWeights_GWAS(ind_k),ep+geneWeights_treatment(ind_k),geneNames{ind_k});
end

[~,highGWAS_ind] = sort(geneWeights_GWAS,'descend');
for k = 1:20
    ind_k = highGWAS_ind(k);
    % Annotate text labels:
    if geneWeights_treatment(ind_k)==0
        ep = 0.01*max(geneWeights_treatment);
        text(ep+geneWeights_GWAS(ind_k),ep+geneWeights_treatment(ind_k),geneNames{ind_k},...
                    'color','r','rotation',90);
    end
end

% Add a distribution density underneath:
[f,x] = ksdensity(geneWeights_GWAS,'npoints',500);
plot(x,-f/max(f)*max(geneWeights_treatment)*0.5);
ax = gca;
ax.YLim = [-max(geneWeights_treatment)*0.5,max(geneWeights_treatment)];
ax.XLim = [min(geneWeights_GWAS),max(geneWeights_GWAS)];

% Labels and title:
xlabel(sprintf('gene score-%sGWAS-%s',whatDiseaseGWAS,whatProperty))
ylabel(sprintf('gene score-%sdrugs',whatDiseaseDrug))
[rho,p] = corr(geneWeights_GWAS,geneWeights_treatment,'type','Kendall');
title(sprintf('Kendall-tau = %.2f (p = %.3g)\n',rho,p));

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
title(sprintf('%s-%s (%s drugs); %u genes',whatDiseaseGWAS,whatProperty,whatDiseaseDrug,sum(bothGood)))
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
