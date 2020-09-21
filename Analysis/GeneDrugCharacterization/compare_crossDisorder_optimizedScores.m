function compare_crossDisorder_optimizedScores(whatMeasures, disorderTrain, disorderTest)
%X - predictor data: all gene scores from GWAS
%y - response: drug score
if nargin < 1
    whatMeasures = 'all'; % 'reduced' or 'all';
end

whatNull = 'randomDrugP';
params = SetDefaultParams();
whatNorm = params.whatNorm;


[geneWeightsGWAS_ALL_train, drugScores_ord, similarityTypes, PPImeasures_names] = give_GWASandDRUG_scores(disorderTrain, whatMeasures);
% select other disorder scores to compare to
geneWeightsGWAS_ALL_test = give_GWASandDRUG_scores(disorderTest, whatMeasures);

% when different disorders are compared, remove all values that are NaN at least in one; 
INDnanTrain = all(isnan(geneWeightsGWAS_ALL_train),1);
INDnanTest = all(isnan(geneWeightsGWAS_ALL_test),1);

INDnan = INDnanTrain | INDnanTest; 
geneWeightsGWAS_ALL_train(:, INDnan) = [];
geneWeightsGWAS_ALL_test(:, INDnan) = [];


% plot randomNullP-based p-values for each mapping method
Dname = disorderTest(isstrprop(disorderTest,'alpha'));

[diseaseResultsR, diseaseResultsP, ~,~, measureNames] = compareGWASvsDRUGmatches({disorderTest}, whatNull, Dname, PPImeasures_names, similarityTypes);
close all;
% replace p=0 to the lowest possible with this number of nulls to avoid InF
% in plotting
diseaseResultsP(diseaseResultsP==0) = 1/params.numNull; 
Pvals_mapp = -log10(diseaseResultsP(~isnan(diseaseResultsP(:))))';
[Pvals_mapp,sIND]  = unique(Pvals_mapp);
measureNames = measureNames(sIND);

colors = zeros(length(Pvals_mapp), 3);
for tt=1:length(Pvals_mapp)
    if contains(measureNames{tt}, 'PPI')
        colors(tt,:) = [50/255,136/255,189/255]; % blue
    elseif contains(measureNames{tt}, 'Allen')
        colors(tt,:) = [.45 .45 .45]; % grey
    elseif contains(measureNames{tt}, 'eQTL')
        colors(tt,:) = [153/255,213/255,148/255]; % green
    else % chromatin and MAGMAdefault
        colors(tt,:) = [254/255,224/255,139/255]; % yellow
    end
end

figure('color','w', 'Position', [100 100 500 500]);

for ww = 1:length(Pvals_mapp)
    plot([Pvals_mapp(ww); Pvals_mapp(ww)], repmat(ylim',1,size(Pvals_mapp,2)), ...
        'LineWidth', 3, 'Color', colors(ww,:), 'LineStyle', ':')
    hold on;
end

xlabel('-log10(P)');
title(sprintf('%s trained on %s', disorderTest, disorderTrain));
xlim([0 3.5])
hold on;

% apply linear regression
mdl = fitlm(geneWeightsGWAS_ALL_train,drugScores_ord);

yTest = predict(mdl,geneWeightsGWAS_ALL_test);
yTestNorm = normalizeScoreVector(yTest, whatNorm);

Pval_comb = compare_to_null(disorderTest, yTestNorm, drugScores_ord, whatNull);
pPLOT = -log10(Pval_comb);

% plot randomNullP-based p-value for combined measure
plot([pPLOT; pPLOT], repmat(ylim',1,size(pPLOT,2)), 'LineWidth', 3, 'Color', [227/255,74/255,51/255]);

% save to file
figureName = sprintf('figures/Train%s_Test%s_compareOptimised_%s', disorderTrain, disorderTest, whatMeasures);
print(gcf,figureName,'-dpng','-r300');


end


