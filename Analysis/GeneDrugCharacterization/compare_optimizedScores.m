function [Ptable] = compare_optimizedScores(whatDiseases_GWAS, whatMeasures, doPlot)
%X - predictor data: all gene scores from GWAS
%y - response: drug score
if nargin < 1
    params = SetDefaultParams();
    whatDiseases_GWAS = params.whatGWAS;
end
if nargin < 2
    whatMeasures = 'reduced'; % 'reduced' or 'all';
end

if nargin < 3
    doPlot = true;
end

whatNull = 'randomDrugP';
params = SetDefaultParams();
whatNorm = params.whatNorm;

numGWAS = length(whatDiseases_GWAS);

% run lasso for each GWAS
for i = 1:numGWAS
    
    whatGWAS = whatDiseases_GWAS{i};
    [geneWeightsGWAS_ALL, drugScores_ord, similarityTypes, PPImeasures_names] = give_GWASandDRUG_scores(whatGWAS, whatMeasures);
    
    % remove rows that are all NaN and measure names
    INDnan = all(isnan(geneWeightsGWAS_ALL),1);
    geneWeightsGWAS_ALL(:, INDnan) = [];
    
    % plot randomNullP-based p-values for each mapping method
    Dname = whatGWAS(isstrprop(whatGWAS,'alpha'));
    
    [~, diseaseResultsP, ~,~, measureNames] = compareGWASvsDRUGmatches({whatGWAS}, whatNull, Dname, PPImeasures_names, similarityTypes);
    close all;
    
    Pvals_mapp = -log10(diseaseResultsP(~isnan(diseaseResultsP(:))))';
    %[Pvals_mapp,sIND]  = unique(Pvals_mapp);
    %measureNames = measureNames(sIND);
    
    if doPlot
        
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
        title(sprintf('%s', whatGWAS));
        xlim([0 4])
        hold on;
    end
    
    % apply linear regression
    mdl = fitlm(geneWeightsGWAS_ALL,drugScores_ord);
    ypred = predict(mdl,geneWeightsGWAS_ALL);
    ypredNorm = normalizeScoreVector(ypred, whatNorm);
    
    Pval_comb = compare_to_null(whatGWAS, ypredNorm, drugScores_ord, whatNull);
    if Pval_comb==0
        Pval_comb = 1/params.numNull;
    end
    pPLOT = -log10(Pval_comb);
    
    measureNames{length(measureNames)+1} = 'Combined';
    Pvals_mapp(length(Pvals_mapp)+1) = pPLOT;
    
    Ptable.(whatGWAS).measureName = measureNames;
    Ptable.(whatGWAS).Pvals = Pvals_mapp;
    
    clearvars Pvals_mapp measureNames
    
    % plot randomNullP-based p-value for combined measure
    if doPlot
        plot([pPLOT; pPLOT], repmat(ylim',1,size(pPLOT,2)), 'LineWidth', 3, 'Color', [227/255,74/255,51/255]);
        % save to file
        figureName = sprintf('figures/%s_compareOptimised_%s', whatGWAS, whatMeasures);
        print(gcf,figureName,'-dpng','-r300');
    end
    
end

end


